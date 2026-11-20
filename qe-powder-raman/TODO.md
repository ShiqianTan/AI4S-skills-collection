# qe-powder-raman 待办

最近审查：2026-11-20

现有主流程已可运行。优先修复会接收无效输入、混用旧结果或把不完整计算标为成功的问题，再完善物理验证、谱图和报告。当前方法范围是三维绝缘晶体的 Γ 点一级、非共振、谐波粉末拉曼，使用已验证的 NC + LDA 路线。

## 验证基线

- `python tests/run_tests.py` 在本机 Conda Python 3.11.15 环境下通过全部 46 项测试。ASE 3.28.0、spglib 2.7.0。系统 `python3` 缺少 ASE，直接运行会产生依赖相关失败。现有集成测试基于 LiF、AlAs、MoS₂ 的日志夹具，主要覆盖解析。
- 本机 QE 8.0dev 实跑 LiF 和 AlAs 的 `all --execute`，使用 `--np 1`、`kpt_density: 6.0`（2 × 2 × 2 网格）和 CSV 输出。LiF 得到 6 支模式、0 支拉曼活性模式。AlAs 得到 3 支拉曼活性模式、两个 LO 方向的输出及 Stokes 谱。粗网格用于验证流程，峰位和强度仍需收敛验证。

## P0 · 输入、执行与结果完整性

- [ ] **P0-01 拒收部分占位和不适用的结构。** 将 LiF 的 Li 占位改为 `0.50` 后，ASE 返回 `{'0': {'Li': 0.5}, '1': {'F': 1}}`，预检仍打印 `all sites fully occupied`。当前代码遍历字典键，并在转换异常后默认通过。统一 `preflight`、`structure`、`prepare` 的结构校验，检查占位、物种混占、三维周期性、晶胞和坐标的有效性。**验收**：完整占位通过，零/部分占位、混占、非三维周期和退化晶胞给出明确错误，失败时不生成 SCF 输入。见 [qe_struct.py](scripts/qe_struct.py#L70)。

- [ ] **P0-02 按 UPF 内容验证赝势。** 复现：将 `F.pz-hgh.UPF` 保存为 `Li.fake.UPF` 后，Li 的清单记录 `element=F`，`errors` 仍为空。核对文件头中的元素、明确的 NC 类型和受支持的 XC 标识，处理同元素多个候选，并使输入生成直接使用已验证的映射。**验收**：错误元素、未知类型、US/PAW、非受支持的 XC 均被拒绝，候选歧义有可执行的解决提示，内置 10 种赝势合格。见 [qe_pp.py](scripts/qe_pp.py#L75)。

- [ ] **P0-03 统一配置解析和路径语义。** 已复现 `prepare` 忽略 YAML 的 `pseudo_dir`，以及同时给出 `--outdir` 与 YAML `outdir` 时，生成目录与 QE `outdir` 值分离。`all` 的准备阶段也优先取 YAML 的 `bin_dir`，执行阶段优先取 CLI。区分项目输出目录与 QE scratch，统一 CLI、YAML、环境变量和默认值的优先级，明确相对路径的基准及旧配置兼容方式。把固定占据、Γ 点流程以及数值字段的允许值集中校验。**验收**：在不同工作目录调用，`prepare`、`all`、`run` 和生成的 `run.sh` 使用一致的有效配置。非法占据、非有限数值、无效 MPI/OMP 参数有明确诊断。见 [qe_raman.py](scripts/qe_raman.py#L357)、[qe_common.py](scripts/qe_common.py#L97)、[qe_inputs.py](scripts/qe_inputs.py#L122)。

- [ ] **P0-04 校验阶段完成条件和产物。** 已复现只保留 AlAs `dynmat.out` 前两条模式、移除 `JOB DONE` 并缺少 `scf.out` 时，`parse` 仍返回 `success`。结合各阶段完成/收敛标志、当前运行记录和关键产物验证结果，依据实际计算晶胞校验 TO 及所请求方向的 `3N` 行、模式编号和必要数值。为导入不全的历史日志保留明确的诊断状态，并区分已验证的全零拉曼与张量生成状态未知。**验收**：截断日志、缺失矩阵、部分 LO 输出和失效的上游产物均不得标为成功，现有 LiF 负对照保留，无法确认拉曼张量时报告明确标明。见 [qe_run.py](scripts/qe_run.py#L24)、[qe_parse.py](scripts/qe_parse.py#L194)、[qe_report.py](scripts/qe_report.py#L113)。

- [ ] **P0-05 让重复执行和中断恢复可靠。** 用返回 0 的测试替身连续三次运行 `--skip-existing` 得到 `success → skipped → success`，且指定的 OMP 线程未传给子进程。还复现了同目录从两个 LO 方向改为无 LO 后旧输入仍被保留，运行器和解析器继续通过通配符扫描。为每阶段原子保存状态，依据当前输入、上游状态及产物决定能否复用。统一 Python 和 `run.sh` 的执行语义，处理默认清理行为及同目录并发。**验收**：第三次仍能跳过，变更上游参数使相关阶段重跑，中断可从有效阶段恢复，关闭 LO 后旧方向不进入报告，子进程收到实际设置的 OMP 环境。见 [qe_run.py](scripts/qe_run.py#L24)、[qe_inputs.py](scripts/qe_inputs.py#L199)、[run.sh.tmpl](templates/run.sh.tmpl)。

## P1 · 物理判据、可复现性与文档

- [ ] **P1-01 加入电子结构适用性判定。** `fixed` 和可能磁性元素的提示只能提供部分信息，现有 SCF 夹具也仅输出最高占据能级。确定在启动 `ph.x` 前获取占据与未占据能级的方式，记录带隙判据、可获得的磁性证据和无法判定的情形。**验收**：已知有隙样例、零带隙/重叠样例、缺少空带信息的旧日志分别得到明确结论，报告能追溯判断依据。见 [qe_struct.py](scripts/qe_struct.py#L97)、[qe_parse.py](scripts/qe_parse.py#L21)、[QE `ph.x` 输入说明](https://www.quantum-espresso.org/Doc/INPUT_PH.html)。

- [ ] **P1-02 统一声学模、软模与谱峰筛选。** 解析器按绝对频率最小的三支判定声学模，绘谱则仅用频率和活性阈值。建立共享的模式分类，优先用振型相对平移的投影识别 Γ 点声学模，旧日志缺少振型时注明近似。区分数值噪声、小虚频和显著光学软模，报告列出虚频与可能原因。**验收**：低频光学模可保留进谱，数值偏移的声学模不入谱，CSV、图和报告使用相同的模式判据，MoS₂ 软模夹具仍发出警告。见 [qe_parse.py](scripts/qe_parse.py#L152)、[qe_spectrum.py](scripts/qe_spectrum.py#L94)、[QE PHonon 用户指南](https://www.quantum-espresso.org/Doc/user_guide_PDF/ph_user_guide.pdf)。

- [ ] **P1-03 记录并验证 LO 方向。** 校验三分量笛卡尔向量的长度、有限性与非零性，记录对应的标准化晶胞、方向向量及文件映射。当前多方向计算各有 `.dat` 和日志，但沿用同名的 `dynmat.mold`、`dynmat.axsf`，会覆盖早先结果。报告按方向展示非解析修正后的模式，标明默认谱来自未加该修正的 Γ 点结果。若提供逐模 LO–TO 差值，用振型或简并子空间匹配替代按行号配对。**验收**：AlAs 立方等价方向结果在收敛后的容差内一致，CdS 的 `q∥c` 与 `q⊥c` 按实际晶胞构造并可辨识，所保留的振型文件按方向区分。见 [qe_inputs.py](scripts/qe_inputs.py#L174)、[qe_report.py](scripts/qe_report.py#L151)、[QE `dynmat.x` 输入说明](https://www.quantum-espresso.org/Doc/INPUT_DYNMAT.html)。

- [ ] **P1-04 明确展宽、Stokes 与局域场约定。** 已复现洛伦兹 `FWHM=0` 生成 4,000 个 NaN。校验 FWHM、温度、激光波长、活动阈值和输出格式，对不满足 `0 < 谱移 < 激光波数` 的 Stokes 模式给出诊断。明确当前展宽按峰高缩放，若增加面积归一化需保留旧行为并记录选项。按波数写清 Bose 指数 `hc·ν̃/(kBT)`。本地 QE 源码 `LR_Modules/dynmat_sub.f90` 的 `RamanIR` 已以 `cmfac**2` 打印 `multiply Raman by` 因子，明确可选修正如何应用一次。**验收**：非法参数在写谱前报错，合成单峰验算峰位、半高宽与权重，报告记录温度、激光、展宽方式和局域场选项。见 [qe_spectrum.py](scripts/qe_spectrum.py#L64)、[理论说明](references/raman_theory.md)。

- [ ] **P1-05 持久化可复现的运行信息。** 记录输入和标准化结构、空间群/容差与精化结果、QE 版本、UPF 来源及 SHA-256、有效配置、命令、k 网格、截断能、收敛阈值与谱图设置。记录 `refine_cell` 因原子数变化或错误回退的原因。**验收**：只凭运行目录即可追溯一个峰所用的计算与后处理设置，旧目录能读取，缺失字段显示为未知。见 [qe_struct.py](scripts/qe_struct.py#L228)、[qe_pp.py](scripts/qe_pp.py#L127)、[qe_report.py](scripts/qe_report.py#L64)。

- [ ] **P1-06 修正技能入口、安装说明与 CLI 示例。** `skill-creator` 的 `quick_validate.py` 因 `SKILL.md` 的 `description` 含 `->` 里的 `>` 拒绝通过。README 的 `prepare --lo 0 0 1` 退出码为 2，目前主 CLI 通过 `lo_to_directions` 读取方向。精简描述并修订配置示例与命令，写明 Python 环境和 QE 的可移植发现方式，修正缺失的 `*_sym.cif` 引用。**验收**：技能通过 `quick_validate.py`，在新目录按文档完成预检与输入生成，方法限制、已验证组合和已实现检查的说明分别可定位。见 [SKILL.md](SKILL.md)、[README.md](README.md)、[config.example.yaml](config.example.yaml)、[示例说明](examples/README.md)。

## P2 · 计算验证与方法扩展

- [ ] **P2-01 建立分层回归与发布基准。** P0 修复随附针对性回归。再增设可选真 QE 测试，用 LiF 负对照、AlAs 活性/LO 和 CdS 或 MoS₂ 晶胞精化覆盖 CLI 主流程。为六个示例补齐结构版本、方法参数与预期结果，对每个峰位、简并度、强度排序注明容差和来源。**验收**：无 QE 的测试覆盖失败路径和产物链，真 QE 测试至少在本机已验证构建及一个明确的稳定版组合上有可查记录，频率回归与实验对照使用各自的参考值。见 [现有集成测试](tests/integration/test_workflow.py)、[兼容性记录](references/qe_compatibility.md)。

- [ ] **P2-02 增加可选的弛豫与收敛验证。** 当前 `--relax` 只生成 `relax.in`。明确固定晶胞和变胞方案，验证弛豫成功后再用其结构生成 SCF，记录残余力、应力和对称性变化。核对斜晶胞 `kpt_density` 的现有实空间公式与文档单位，引入明确的倒空间网格间距或显式网格，并评估兼容性。**验收**：同一小体系可以复现峰位与主要活性对 `ecutwfc`、k 网格和 `tr2_ph` 的变化，弛豫失败有明确状态。见 [qe_inputs.py](scripts/qe_inputs.py#L68)、[relax 模板](templates/relax.in.tmpl)。

- [ ] **P2-03 评估极性粉末的方向平均。** 待 P1-03 的方向、张量及模式追踪可靠后，比较方位采样、权重与实验散射几何对谱形的影响。**验收**：选用明确的粉末散射几何，在立方材料有对称性校验，在各向异性材料有采样收敛记录，再提供可公开使用的选项。见 [兼容性说明](references/qe_compatibility.md)、[QE 单 q 点计算说明](https://www.quantum-espresso.org/Doc/ph_user_guide/node8.html)。

## 维护方式

- 修改输入校验、计算流程或谱图时，同步更新相关条目的复现、状态与验收证据。新增条目沿用唯一编号，记录影响和最小复现。发现新证据时可调整优先级，保留编号便于追踪。
- `[ ]` 表示验收尚未完成。相关回归通过后勾选，移至下方并记录日期、改动位置和验证命令或结果。无法继续时在原条目记录具体阻塞及恢复条件。
- 每次复审更新顶部日期与验证基线。真 QE 测试注明版本、体系及关键参数，夹具测试注明覆盖路径。

## 完成记录

2026-11-20 建立清单，完成代码审查、缺陷复现与基线验证。
