# qe-powder-raman 待办

最近审查：2026-11-20

现有主流程已可运行。优先修复会接收无效输入、混用旧结果或把不完整计算标为成功的问题，再完善物理验证、谱图和报告。当前方法范围是三维绝缘晶体的 Γ 点一级、非共振、谐波粉末拉曼，使用已验证的 NC + LDA 路线。

## 验证基线

- `python tests/run_tests.py` 在本机 Conda Python 3.11.15 环境下通过全部 165 项测试（2 项为可选真 QE 冒烟测试，未设 `QE_REAL_TESTS=1` 时跳过）。ASE 3.28.0、spglib 2.7.0。系统 `python3` 缺少 ASE，直接运行会产生依赖相关失败。集成测试基于 LiF、AlAs、MoS₂ 的日志夹具覆盖解析与失败路径；单元测试覆盖配置、结构、UPF、输入生成、执行、谱图与报告。
- 本机 QE 8.0dev 已在 6 个 MP-20 体系（LiF、MgH₂、Al₂O₃、AlAs、CdS、MoS₂）上完成 `pw.x → ph.x → dynmat.x` 实跑，日志与回归数值见 `examples/README.md` 的「预期结果与容差（回归基准）」。真实流程冒烟：LiF 与 AlAs 的 `all --execute`（`--np 1`、`kpt_density: 6.0`，即 2 × 2 × 2 网格）。LiF 得到 6 支模式、0 支拉曼活性模式；AlAs 得到 3 支拉曼活性模式、两个 LO 方向输出及 Stokes 谱。粗网格仅用于验证流程，峰位和强度仍需收敛验证。

## P0 · 输入、执行与结果完整性

- [x] **P0-01 拒收部分占位和不适用的结构。** 将 LiF 的 Li 占位改为 `0.50` 后，ASE 返回 `{'0': {'Li': 0.5}, '1': {'F': 1}}`，预检仍打印 `all sites fully occupied`。当前代码遍历字典键，并在转换异常后默认通过。统一 `preflight`、`structure`、`prepare` 的结构校验，检查占位、物种混占、三维周期性、晶胞和坐标的有效性。**验收**：完整占位通过，零/部分占位、混占、非三维周期和退化晶胞给出明确错误，失败时不生成 SCF 输入。见 [qe_struct.py](scripts/qe_struct.py#L70)。

- [x] **P0-02 按 UPF 内容验证赝势。** 复现：将 `F.pz-hgh.UPF` 保存为 `Li.fake.UPF` 后，Li 的清单记录 `element=F`，`errors` 仍为空。核对文件头中的元素、明确的 NC 类型和受支持的 XC 标识，处理同元素多个候选，并使输入生成直接使用已验证的映射。**验收**：错误元素、未知类型、US/PAW、非受支持的 XC 均被拒绝，候选歧义有可执行的解决提示，内置 10 种赝势合格。见 [qe_pp.py](scripts/qe_pp.py#L75)。

- [x] **P0-03 统一配置解析和路径语义。** 已复现 `prepare` 忽略 YAML 的 `pseudo_dir`，以及同时给出 `--outdir` 与 YAML `outdir` 时，生成目录与 QE `outdir` 值分离。`all` 的准备阶段也优先取 YAML 的 `bin_dir`，执行阶段优先取 CLI。区分项目输出目录与 QE scratch，统一 CLI、YAML、环境变量和默认值的优先级，明确相对路径的基准及旧配置兼容方式。把固定占据、Γ 点流程以及数值字段的允许值集中校验。**验收**：在不同工作目录调用，`prepare`、`all`、`run` 和生成的 `run.sh` 使用一致的有效配置。非法占据、非有限数值、无效 MPI/OMP 参数有明确诊断。见 [qe_raman.py](scripts/qe_raman.py#L357)、[qe_common.py](scripts/qe_common.py#L97)、[qe_inputs.py](scripts/qe_inputs.py#L122)。

- [x] **P0-04 校验阶段完成条件和产物。** 已复现只保留 AlAs `dynmat.out` 前两条模式、移除 `JOB DONE` 并缺少 `scf.out` 时，`parse` 仍返回 `success`。结合各阶段完成/收敛标志、当前运行记录和关键产物验证结果，依据实际计算晶胞校验 TO 及所请求方向的 `3N` 行、模式编号和必要数值。为导入不全的历史日志保留明确的诊断状态，并区分已验证的全零拉曼与张量生成状态未知。**验收**：截断日志、缺失矩阵、部分 LO 输出和失效的上游产物均不得标为成功，现有 LiF 负对照保留，无法确认拉曼张量时报告明确标明。见 [qe_run.py](scripts/qe_run.py#L24)、[qe_parse.py](scripts/qe_parse.py#L194)、[qe_report.py](scripts/qe_report.py#L113)。

- [x] **P0-05 让重复执行和中断恢复可靠。** 用返回 0 的测试替身连续三次运行 `--skip-existing` 得到 `success → skipped → success`，且指定的 OMP 线程未传给子进程。还复现了同目录从两个 LO 方向改为无 LO 后旧输入仍被保留，运行器和解析器继续通过通配符扫描。为每阶段原子保存状态，依据当前输入、上游状态及产物决定能否复用。统一 Python 和 `run.sh` 的执行语义，处理默认清理行为及同目录并发。**验收**：第三次仍能跳过，变更上游参数使相关阶段重跑，中断可从有效阶段恢复，关闭 LO 后旧方向不进入报告，子进程收到实际设置的 OMP 环境。见 [qe_run.py](scripts/qe_run.py#L24)、[qe_inputs.py](scripts/qe_inputs.py#L199)、[run.sh.tmpl](templates/run.sh.tmpl)。

## P1 · 物理判据、可复现性与文档

- [x] **P1-01 加入电子结构适用性判定。** `fixed` 和可能磁性元素的提示只能提供部分信息，现有 SCF 夹具也仅输出最高占据能级。确定在启动 `ph.x` 前获取占据与未占据能级的方式，记录带隙判据、可获得的磁性证据和无法判定的情形。**验收**：已知有隙样例、零带隙/重叠样例、缺少空带信息的旧日志分别得到明确结论，报告能追溯判断依据。见 [qe_struct.py](scripts/qe_struct.py#L97)、[qe_parse.py](scripts/qe_parse.py#L21)、[QE `ph.x` 输入说明](https://www.quantum-espresso.org/Doc/INPUT_PH.html)。

- [x] **P1-02 统一声学模、软模与谱峰筛选。** 解析器按绝对频率最小的三支判定声学模，绘谱则仅用频率和活性阈值。建立共享的模式分类，优先用振型相对平移的投影识别 Γ 点声学模，旧日志缺少振型时注明近似。区分数值噪声、小虚频和显著光学软模，报告列出虚频与可能原因。**验收**：低频光学模可保留进谱，数值偏移的声学模不入谱，CSV、图和报告使用相同的模式判据，MoS₂ 软模夹具仍发出警告。见 [qe_parse.py](scripts/qe_parse.py#L152)、[qe_spectrum.py](scripts/qe_spectrum.py#L94)、[QE PHonon 用户指南](https://www.quantum-espresso.org/Doc/user_guide_PDF/ph_user_guide.pdf)。

- [x] **P1-03 记录并验证 LO 方向。** 校验三分量笛卡尔向量的长度、有限性与非零性，记录对应的标准化晶胞、方向向量及文件映射。当前多方向计算各有 `.dat` 和日志，但沿用同名的 `dynmat.mold`、`dynmat.axsf`，会覆盖早先结果。报告按方向展示非解析修正后的模式，标明默认谱来自未加该修正的 Γ 点结果。若提供逐模 LO–TO 差值，用振型或简并子空间匹配替代按行号配对。**验收**：AlAs 立方等价方向结果在收敛后的容差内一致，CdS 的 `q∥c` 与 `q⊥c` 按实际晶胞构造并可辨识，所保留的振型文件按方向区分。见 [qe_inputs.py](scripts/qe_inputs.py#L174)、[qe_report.py](scripts/qe_report.py#L151)、[QE `dynmat.x` 输入说明](https://www.quantum-espresso.org/Doc/INPUT_DYNMAT.html)。

- [x] **P1-04 明确展宽、Stokes 与局域场约定。** 已复现洛伦兹 `FWHM=0` 生成 4,000 个 NaN。校验 FWHM、温度、激光波长、活动阈值和输出格式，对不满足 `0 < 谱移 < 激光波数` 的 Stokes 模式给出诊断。明确当前展宽按峰高缩放，若增加面积归一化需保留旧行为并记录选项。按波数写清 Bose 指数 `hc·ν̃/(kBT)`。本地 QE 源码 `LR_Modules/dynmat_sub.f90` 的 `RamanIR` 已以 `cmfac**2` 打印 `multiply Raman by` 因子，明确可选修正如何应用一次。**验收**：非法参数在写谱前报错，合成单峰验算峰位、半高宽与权重，报告记录温度、激光、展宽方式和局域场选项。见 [qe_spectrum.py](scripts/qe_spectrum.py#L64)、[理论说明](references/raman_theory.md)。

- [x] **P1-05 持久化可复现的运行信息。** 记录输入和标准化结构、空间群/容差与精化结果、QE 版本、UPF 来源及 SHA-256、有效配置、命令、k 网格、截断能、收敛阈值与谱图设置。记录 `refine_cell` 因原子数变化或错误回退的原因。**验收**：只凭运行目录即可追溯一个峰所用的计算与后处理设置，旧目录能读取，缺失字段显示为未知。见 [qe_struct.py](scripts/qe_struct.py#L228)、[qe_pp.py](scripts/qe_pp.py#L127)、[qe_report.py](scripts/qe_report.py#L64)。

- [x] **P1-06 修正技能入口、安装说明与 CLI 示例。** `skill-creator` 的 `quick_validate.py` 因 `SKILL.md` 的 `description` 含 `->` 里的 `>` 拒绝通过。README 的 `prepare --lo 0 0 1` 退出码为 2，目前主 CLI 通过 `lo_to_directions` 读取方向。精简描述并修订配置示例与命令，写明 Python 环境和 QE 的可移植发现方式，修正缺失的 `*_sym.cif` 引用。**验收**：技能通过 `quick_validate.py`，在新目录按文档完成预检与输入生成，方法限制、已验证组合和已实现检查的说明分别可定位。见 [SKILL.md](SKILL.md)、[README.md](README.md)、[config.example.yaml](config.example.yaml)、[示例说明](examples/README.md)。

## P2 · 计算验证与方法扩展

- [x] **P2-01 建立分层回归与发布基准。** P0 修复随附针对性回归。再增设可选真 QE 测试，用 LiF 负对照、AlAs 活性/LO 和 CdS 或 MoS₂ 晶胞精化覆盖 CLI 主流程。为六个示例补齐结构版本、方法参数与预期结果，对每个峰位、简并度、强度排序注明容差和来源。**验收**：无 QE 的测试覆盖失败路径和产物链，真 QE 测试至少在本机已验证构建及一个明确的稳定版组合上有可查记录，频率回归与实验对照使用各自的参考值。见 [现有集成测试](tests/integration/test_workflow.py)、[兼容性记录](references/qe_compatibility.md)。

- [x] **P2-02 增加可选的弛豫与收敛验证。** 当前 `--relax` 只生成 `relax.in`。明确固定晶胞和变胞方案，验证弛豫成功后再用其结构生成 SCF，记录残余力、应力和对称性变化。核对斜晶胞 `kpt_density` 的现有实空间公式与文档单位，引入明确的倒空间网格间距或显式网格，并评估兼容性。**验收**：同一小体系可以复现峰位与主要活性对 `ecutwfc`、k 网格和 `tr2_ph` 的变化，弛豫失败有明确状态。见 [qe_inputs.py](scripts/qe_inputs.py#L68)、[relax 模板](templates/relax.in.tmpl)。

- [x] **P2-03 评估极性粉末的方向平均。** 待 P1-03 的方向、张量及模式追踪可靠后，比较方位采样、权重与实验散射几何对谱形的影响。**验收**：选用明确的粉末散射几何，在立方材料有对称性校验，在各向异性材料有采样收敛记录，再提供可公开使用的选项。见 [兼容性说明](references/qe_compatibility.md)、[QE 单 q 点计算说明](https://www.quantum-espresso.org/Doc/ph_user_guide/node8.html)。

## P3 · Phase D 可信度收尾（已完成）

- [x] **D-01 统一 Git 源与同步方式。** 工作目录自带 `.git`（本地快照、无 remote），权威发布源是 `AI4S-skills-collection` 仓库；当前 rsync 未排除 `.git`，可能把本地快照带入发布子目录。**验收**：同步命令显式排除 `.git`、`__pycache__`、`*.pyc`、`.DS_Store`；文档或 TODO 记录“工作目录为编辑副本、远程仓库为唯一发布源”；本地快照与当前文件内容一致（无未并入库的改动）。
- [x] **D-02 弛豫结果进入解析与报告。** `relax_manifest.json` 已记录残余力/应力/对称性，但 `parse`/`report` 未读取。**验收**：`summary.json` 含 `relax` 字段；`report.md` 展示是否执行弛豫、残余力/应力、对称性；输出文件清单包含 `relax_manifest.json`；缺少该文件时显示“未弛豫”而非报错。
- [x] **D-03 修正 `--relax` 帮助文本。** `prepare` 与 `all` 的帮助仍写 “also write relax.in”，而 `run`/`all --execute` 现在会执行并重建输入。**验收**：三处帮助分别说明“仅生成 relax.in”“执行弛豫并重建 SCF/ph/dynmat”，与 `run --relax` 行为一致。
- [x] **D-04 校验 `lraman` 生效与 Raman 张量完整性。** 当前仅依据 dynmat 是否存在 Raman 列判断，无法区分“全零张量”和“张量未生成”。**验收**：`parse` 检查 `ph.out` 的 `Raman tensor`/`RAMAN COEFFICIENTS` 块；缺失且无 Raman 列时状态为 `failed`；有列但缺张量证据时同样失败；LiF 全零负对照仍为 `success`。`summary.json` 记录 `lraman_tensor`。
- [x] **D-05 CLI 支持 `--kpt-density` / `--kpt-grid`。** 默认 `kpt_density=30` 与示例基准 6×6×6 不一致，且无法从 CLI 直接复现基准。**验收**：`prepare`/`all` 接受 `--kpt-density FLOAT` 与 `--kpt-grid NX NY NZ`（显式网格优先）；生成输入与 `inputs_manifest.json` 使用覆盖后的值；有单元测试覆盖优先级。
- [x] **D-06 QE 版本门槛。** `preflight` 只显示版本，不校验最低版本。**验收**：解析 `Program PWSCF v.X.Y`；低于最低支持版本（7.0）时 `FAIL`；无法解析时保持 `WARN`；有单元测试覆盖版本解析。

## P4 · Phase E 方法能力（已确认引入 PBE NC）

- [x] **E-01 收敛扫描子命令。** 新增 `converge`：对 `ecutwfc`、`kpt_density`、`tr2_ph` 逐点建目录跑完整 `pw.x → ph.x → dynmat.x`，`fwhm` 只做后处理展宽；输出 `convergence.{json,csv,md}`，判据为峰位 ≤ 2 cm⁻¹、Raman 活性相对变化 ≤ 5%、总能 ≤ 1e-4 Ry/atom、介电常数相对变化 ≤ 1%，并给出推荐最小可用值；`report.md` 增加“收敛扫描”章节。已用 LiF 30/45 Ry 实跑复验（ε∞ 4.08→3.37，最大峰移 15.9 cm⁻¹，判为未收敛），另有单元/集成回归覆盖。
- [x] **E-02 引入 PBE 模守恒赝势。** 已核实 QE 的 `lraman` **只支持 LDA**：`PHonon/PH/phq_setup.f90` 对 GGA 直接报 `third order derivatives not implemented with GGA`；实测 MoS₂+PBE NC 的 `pw.x` 正常、`ph.x` 立即失败。因此调整为：内置 PseudoDojo ONCVPSP-PBE-FR-PDv0.4 NC 赝势（`pseudo/pbe/`，10 元素，fully relativistic，来源+SHA-256 记录），赝势校验接受 LDA/GGA NC，manifest 记录 `xc_families`/`raman_supported`/`pseudo_source`；Raman 开启时对 GGA 直接 FAIL 并提示 LDA 或 `--no-lraman`；新增 `raman: false` / `--no-lraman` 的声子/IR/介电工作流。**验收（修订）**：PBE NC 体系可用 `all --execute --no-lraman` 完成全流程（LiF 2×2×2 已复验）；PBE + 默认 Raman 被明确拦截；非 NC/未知 XC 仍被拒绝；报告中记录 XC/来源。
- [x] **E-03 LO/IR 谱与逐方向谱。** `spectrum` 现按每个 LO 表输出 `spectrum_lo[_N].*`（标题含 `q=(x,y,z)`，方向取自 `inputs_manifest.json`），并生成与 Γ 点 TO 叠加的 `spectrum_to_vs_lo[_N].png`；新增 IR 活性谱 `spectrum_ir.*` 与双轴联合图 `spectrum_raman_ir.png`（LiF 无拉曼但有 IR 时也有效）。`--no-ir` / `--no-lo-spectra` 与配置 `spectrum_ir` / `spectrum_lo` 可关闭；AlAs（LO 365 cm⁻¹ vs TO 332 cm⁻¹）与 LiF（拉曼空白、IR 有峰）已回归。
- [x] **E-04 模式指认与实验对照。** 已实现 `assign` 子命令：按 Raman/IR 活性与频率简并给出模式指认（`mode_assignment.csv` / `assignment.md`），支持实验峰表或完整谱（>50 行自动检峰并用 Lorentzian 拟合，scipy 可选、否则矩估计），输出实验匹配偏差、RMS、误差带与 `assignment_overlay.png`；`report.md` 增加第 9 节。未计算不可约表示（需特征向量与空间群特征标表）。
- [x] **E-05 层状/低维提示。** 已实现 `analyze_layers`：沿三个晶向测量最大间隙（真空层 ≥ 8 Å 判为 slab/准二维；≥ 2.5 Å 且两侧同元素判为层状/范德华间隙），估计层数并写入 `structure_manifest.json`；预检给出 `slab_vacuum` / `layered_gap` 警告与 vc-relax/色散修正建议，`report.md` 增加「层状/低维提示」（含软模提示）。

## P5 · Phase F 工程交付

- [ ] **F-01 打包与元数据。** 增加 `pyproject.toml`、console entry point、`__version__`、`LICENSE`、`CHANGELOG`、`CITATION`。
- [ ] **F-02 CI 与质量门禁。** GitHub Actions 运行测试；增加 coverage、ruff/flake8；为 `report.md` 与谱图增加黄金文件测试。
- [ ] **F-03 CLI 自动化与并发。** 增加 `--json`、`--quiet`、日志文件；`state.json` 并发锁；多 LO 方向并行。
- [ ] **F-04 国际化与展示。** 报告/图支持英文；字体回退跨平台；报告嵌入谱图。
- [ ] **F-05 文档收尾。** 修正 `README.md:19` “templates（Step 2 填充）”、Step 0 “本步”等过期表述；补充一份真实运行产物样例（`report.md` + 谱图）与故障排查案例。

## 维护方式

- 修改输入校验、计算流程或谱图时，同步更新相关条目的复现、状态与验收证据。新增条目沿用唯一编号，记录影响和最小复现。发现新证据时可调整优先级，保留编号便于追踪。
- `[ ]` 表示验收尚未完成。相关回归通过后勾选，移至下方并记录日期、改动位置和验证命令或结果。无法继续时在原条目记录具体阻塞及恢复条件。
- 每次复审更新顶部日期与验证基线。真 QE 测试注明版本、体系及关键参数，夹具测试注明覆盖路径。

## 完成记录

### 2026-11-20 · QE 子程序参数与 Raman 知识库

- 参考 ORCA skill 的参数速查和知识库结构，新增 [QE 子程序与调用参数](references/qe-command-line.md)，扩充 [结构化 Raman 知识库](data/raman_knowledge.json) 与 [理论及结果判读](references/raman_theory.md)。记录 QE 7.5 官方手册、原始论文和本机 8.0dev 源码的出处。
- 本机 LiF 验证 `pw.x -nk 1 -i`、`ph.x -npool 1 -in`、`dynmat.x -input` 和 `-inp`。LiH 的真实 Γ 动力学矩阵验证 `q = 0, 0, 1` 与 `q(3) = 1` 的振型输出逐字节一致。源码核对了 Raman 路线的 ASR 和局域场因子。
- 知识库通过 JSON、来源索引和内部引用校验。含真实 QE 用例的 `QE_REAL_TESTS=1 python tests/run_tests.py` 共 171 项通过。


### 2026-11-20 · MP-20 第五轮谱图与报告回归

- 沿用前几轮筛选规则，排除已抽取的 12 个结构和 6 个内置示例，剩 37 条。按 CSV 原顺序使用 `random.Random(20261124).sample(rows, 3)`。实跑采用 QE 8.0dev、HGH NC/LDA、3 MPI ranks、40/160 Ry、`tr2_ph=1e-10`、`kpt_density: 6.0`。
- 用 LiH 复现先生成 CSV/PNG、再用 13 cm⁻¹ 高斯展宽、250 K、785 nm 仅生成 Stokes CSV 时，旧报告仍显示默认参数并列出上一次留下的图。`spectrum` 现在记录 `spectrum_manifest.json`，绑定实际参数、当次文件清单和模式表哈希。报告只展示可核验记录，旧目录可重绘补齐记录。解析失败的模式表不能生成谱图。
- 在 LiH 和 Li₂SO₄ 上验证带覆盖参数的重绘、旧文件保留且从报告移除，以及模式表改变后提示重绘。三个样本均验证了 `3N` 模式和重新生成的 CSV、PNG、报告与清单的一致性。含真 QE 回归的 `QE_REAL_TESTS=1 python tests/run_tests.py` 共 171 项通过。

| MP ID | 结构 | SCF 带隙 (eV) | TO 模式数 | TO 光学虚频数 |
|---|---|---:|---:|---:|
| `mp-558382` | Li₂SO₄ | 5.6599 | 42 | 1 |
| `mp-625700` | MoH₄O₅ | 2.0747 | 60 | 16 |
| `mp-23703` | LiH | 2.0734 | 6 | 3 |

LiH 保持原始晶胞，用 6 × 6 × 6 k 网格和 60/240 Ry 复验，带隙为 2.3975 eV，三支光学模约 481 cm⁻¹，均为正频。粗参数的三支虚频约 −1030 cm⁻¹。此对照同时改变了网格和截断能，尚未分离两者的影响。


### 2026-11-20 · MP-20 第四轮弛豫与对称性回归

- 沿用前几轮筛选规则，排除已抽取的 9 个结构和 6 个示例，剩 40 条。按 CSV 原顺序用 `random.Random(20261123).sample(rows, 3)` 抽取。QE 8.0dev 使用 HGH NC/LDA、3 MPI ranks、40/160 Ry、`tr2_ph=1e-10`、`kpt_density: 6.0`。
- MgAs₄ 默认 `ph.x` 在 `set_irr_sym_new` 报 `wrong representation`。在新目录关闭声子模式的对称性搜索后走通，验证了新增 YAML 键 `ph_search_sym: false` 生成同等 `ph.in`。默认继续启用对称性搜索，遇到该错误时给出复验提示。
- MoS₂ 的 `nstep=1` 真 QE 复现中，旧流程将弛豫阶段记为成功后才因缺失末态退出。修复后记录步数上限导致的未收敛，并保留原 SCF 输入。正常弛豫在 10 个 BFGS 步后收敛，重建的 SCF 晶胞和分数坐标与 QE 末态一致。补充了 `(alat=…)` 的末态晶胞和坐标换算测试。
- 三个抽样样本均校验了 `3N` 模式、两类拉曼谱和报告。表格列出原几何粗网格结果，MgAs₄ 使用声子对称性复验设置。MoS₂ 变胞弛豫对照得到 18 支模式、0 支光学虚频，原结构有 1 支。虚频与收敛仍需独立复核。
- 基于最新远端版本的发布分支与本地编辑副本均通过 169 项测试，含真 QE 用例。新增 [输入预检与回复指南](references/input-triage.md)。

| MP ID | 结构 | SCF 带隙 (eV) | TO 模式数 | TO 光学虚频数 |
|---|---|---:|---:|---:|
| `mp-7623` | MgAs₄ | 0.5894 | 60 | 8 |
| `mp-625316` | Al(OH)₃ | 4.3169 | 42 | 5 |
| `mp-1018809` | MoS₂ | 1.7001 | 18 | 1 |


2026-11-20 Phase D（D-01…D-06）完成：同步排除本地 `.git` 并记录发布方式；弛豫信息写入 `summary.json` 与 `report.md`；修正 `--relax` 帮助；校验 `ph.out` 的 Raman 张量生成（缺失即失败）；新增 `--kpt-density`/`--kpt-grid`；preflight 增加 QE 最低版本 7.0 检查；版本探测改在临时目录执行，避免 `CRASH`/`input_tmp.in` 污染工作目录。全套房具测试 111 项通过（真 QE 默认跳过）。
### 2026-11-20 · MP-20 第三轮续跑回归

- 沿用前两轮筛选规则，排除已抽取的 6 个结构和 6 个内置示例，剩 43 条。按 CSV 原顺序用 `random.Random(20261122).sample(rows, 3)` 抽取。真实 QE 8.0dev 使用 HGH NC/LDA、3 MPI ranks、40/160 Ry、`tr2_ph=1e-10`、`kpt_density: 6.0`，输出两类谱的 CSV 和 PNG。
- 用 AsF₅ 复现旧逻辑在截断能从 40 改为 45 Ry 后依然跳过 SCF，日志停在 40 Ry。新逻辑核对输入和日志哈希、命令、线程配置与关键产物，确认 45 Ry 和 42 Ry 的输入变更各自触发重算，42 Ry 保持不变时可跳过，移走 SCF 保存目录的 XML 后会重算。上游重启时，下游旧记录立即标记为过期，解析会拒绝旧日志。
- AsF₅ 使用相同配置再次运行 `all --execute --skip-existing`，SCF、声子和 `dynmat` 均跳过。单元测试在 `ph.x` 启动时注入中断，确认状态文件保留已完成的 SCF 和正在运行的 `ph.x`，续跑跳过 SCF 并重启 `ph.x`。同时覆盖 QE 零退出但缺少完成标志。含真实 QE LiF 用例的 `QE_REAL_TESTS=1 python tests/run_tests.py` 共 101 项通过。

| MP ID | 结构 | SCF 带隙 (eV) | TO 模式数 | TO 光学虚频数 |
|---|---|---:|---:|---:|
| `mp-8723` | AsF₅ | 4.5100 | 36 | 0 |
| `mp-625731` | MoH₄O₅ | 2.2443 | 60 | 19 |
| `mp-690733` | H₂SO₄ | 5.6220 | 42 | 8 |

三例均生成报告和两类谱，模式数为 `3N`。虚频保留在报告中。原始几何和粗网格仅用于工作流验证。

### 2026-11-20 · MP-20 第二轮 LO 回归

- 沿用上一轮筛选规则，排除前一轮的 3 个结构及 6 个内置示例，剩 46 条。按 CSV 原顺序使用 `random.Random(20261121).sample(rows, 3)` 抽取。三个体系均在本机 QE 8.0dev 跑通，使用 HGH NC/LDA、3 MPI ranks、40/160 Ry、`tr2_ph=1e-10`、`kpt_density: 6.0` 和笛卡尔方向 `q=(0,0,1)`。
- 用真实 MgO 输出复现了缺失或截断的 LO 日志仍通过解析，以及 LO 覆盖 TO 振型文件。解析现在根据输入清单核对所请求方向的日志、结束标志和模式编号，忽略清单外旧 LO 日志。新生成的 LO 输入使用独立 `.mold` 和 `.axsf` 文件名。
- 用三个样本的实际 `.dyn` 重跑 TO 和 LO 后处理，逐一确认缺失日志或移除完成标志时退出码为 1，恢复完整日志后能生成谱和报告，TO 振型文件在 LO 前后逐字节一致。含真实 QE LiF 用例的 `QE_REAL_TESTS=1 python tests/run_tests.py` 共 97 项通过。

| MP ID | 结构 | SCF 带隙 (eV) | TO / LO 模式数 | TO 光学虚频数 |
|---|---|---:|---:|---:|
| `mp-1065577` | MgO | 2.8545 | 12 / 12 | 0 |
| `mp-7744` | LiSO₃F | 5.0055 | 36 / 36 | 6 |
| `mp-675020` | Li₄MoO₅ | 3.2510 | 30 / 30 | 2 |

本轮沿用 CSV 几何，光学虚频保留在报告中，频率和活性仍需弛豫与收敛验证。

### 2026-11-20 · MP-20 随机抽样复验

- 从 MatterGen `data-release/mp-20/mp_20/test.csv` 的 9,046 行中，保留元素均有内置赝势、`dft_band_gap >= 0.5` eV、`abs(dft_mag_density) <= 1e-4` 的 55 行，按 CSV 原顺序用 Python `random.Random(20261120).sample(rows, 3)` 抽取。三份 CIF 来自原行的 `cif` 字段。
- 三个初次 SCF 日志都缺少空带能级。SCF 输入现在依据 UPF 价电子数为中性闭壳层晶胞加入一个空带，并完整收敛空带能级。`run` 与 `all --execute` 在已测带隙 ≤ 0.05 eV 时阻止 `ph.x` 启动，旧日志缺少空带信息时记录警告。见 [QE `pw.x` 的 `nbnd` 说明](https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm277)。
- `QE_REAL_TESTS=1 python tests/run_tests.py` 通过 94 项，包含真实 QE LiF 负对照和重叠能带阻止 `ph.x` 的回归。

| MP ID | 结构 | SCF 带隙 (eV) | 总模式数 | 拉曼活性模式数¹ | 光学虚频数 |
|---|---|---:|---:|---:|---:|
| `mp-27281` | LiH₃O₂ | 4.4370 | 36 | 18 | 1 |
| `mp-625318` | Al(OH)₃ | 4.6098 | 42 | 39 | 1 |
| `mp-715584` | MoO₃ | 1.3299 | 48 | 45 | 16 |

本机 QE 8.0dev 使用内置 HGH NC/LDA、`kpt_density: 6.0`。前两例用 60/240 Ry、`tr2_ph=1e-14`、1 MPI rank。MoO₃ 用 40/160 Ry、`tr2_ph=1e-10`、4 MPI ranks 完成快速端到端复验，其 60/240 Ry、1 rank 的 SCF 另测得 1.4366 eV。三个报告均保留光学虚频警告。检查了模式数 `3N`、两类谱的 PNG 和各 4,000 行有限数值 CSV。此处采用原始几何与粗网格，材料稳定性、峰位及活性仍需弛豫与收敛验证。

¹ 模表在活性阈值以上的光学模式数，含有活性的虚频；绘谱只纳入符合条件的正频模式。

### 前期记录

2026-11-20 建立清单，完成代码审查、缺陷复现与基线验证。
2026-11-20 Phase A（P0-01…P0-05）完成：结构/UPF 内容校验、配置与路径语义、阶段完成/产物检查、可靠重跑与 LO 清理；新增回归测试，全套房具测试通过（提交 `20337e7`）。
2026-11-20 Phase B（P1-01…P1-06）完成：电子结构适用性（gapped/metallic/undeterminable）、共享声学/软模/谱峰分类器、LO 方向校验与逐方向留档、展宽/Stokes/局域场约定、可复现信息持久化、SKILL/README/配置/示例修订；新增回归测试。
2026-11-20 Phase C（P2-01…P2-03）完成：分层回归（含可选真 QE 冒烟 `tests/integration/test_real_qe.py`，`QE_REAL_TESTS=1` 时运行）、六示例「预期结果与容差」基准、`--relax` 实际执行并解析末态结构重生成 SCF（残余力/应力/对称性写入 `relax_manifest.json`）、显式 `kpt_grid` 与 `ecutrho` 提示、极性粉末方向平均评估记录。**遗留**：`ecutwfc`/k 网格/`tr2_ph` 的系统收敛扫描与多方向 LO 收敛记录列为后续工作（P2-02/P2-03 文档已注明）。全套房具测试 92 项通过（真 QE 测试默认跳过）。
2026-11-20 Phase E（E-02）完成：核实 QE `lraman` 仅支持 LDA（`phq_setup.f90` 的 GGA 拦截，MoS₂+PBE 实测复现），将 PBE NC 定位为声子/IR 工作流——内置 `pseudo/pbe/`（PseudoDojo ONCVPSP-PBE-FR-PDv0.4，10 元素，fully relativistic），manifest 记录 XC 家族/来源/SHA-256，Raman 开启时显式拦截 GGA，新增 `raman: false` 与 `--no-lraman`；LiF PBE `all --execute` 全流程复验通过；全套房具测试 122 项通过（2 项真 QE 默认跳过）。
2026-11-20 兼容性修复：`pw.x` HDF5 构建（QE 7.5，`.hdf5` 电荷密度/波函数）下 scf→ph 断点校验误判已修复，`_artifacts()` 同时接受 `.dat` 与 `.hdf5`，并新增伪 QE HDF5 回归测试；测试 165 项通过。
2026-11-20 Phase E（E-04/E-05）完成：新增 `assign` 模式指认与实验对照（活性分类、简并分组、实验匹配、谱峰拟合、叠图、报告第 9 节）以及层状/低维检测（真空层、同元素范德华间隙、层数、结构性警告）；全套房具测试 165 项通过（2 项真 QE 默认跳过）。
2026-11-20 Phase E（E-01）完成：新增 `converge` 收敛扫描（ecutwfc/kpt_density/tr2_ph/FWHM，逐变体目录 + 统一判据 + 推荐值 + `convergence.{json,csv,md}`），`report.md` 汇总收敛结果；LiF 30/45 Ry 真 QE 复验判为未收敛，验证判据与链路。全套房具测试 165 项通过（2 项真 QE 默认跳过）。
2026-11-20 Phase E（E-03）完成：`spectrum` 增加 IR 活性谱与 Raman/IR 联合图、逐 LO 方向谱及 TO/LO 叠加对比图，方向标签来自 `inputs_manifest.json`；新增单元/集成回归（AlAs 的 TO/LO、LiF 的“拉曼空白但 IR 有峰”），全套房具测试 132 项通过（2 项真 QE 默认跳过）。
