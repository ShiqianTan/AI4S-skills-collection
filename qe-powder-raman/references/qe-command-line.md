# Quantum ESPRESSO 子程序与调用参数

适用范围为本 skill 的三维绝缘晶体 Γ 点 Raman。官方输入手册当前页面标识 QE 7.5，下面的本机实现说明依据 QE 8.0dev `215ff8a23`。区分 **shell 启动参数**、**QE 输入文件字段**与 **skill CLI 参数**，重复运行时记录版本、MPI/OpenMP 配置和生成的输入文件。

## 子程序与文件接口

| 程序 | 直接调用 | 读入与输出 |
|---|---|---|
| `pw.x` | `pw.x -in scf.in > scf.out 2>&1` | 做自洽场计算，写 `outdir/prefix.save/`。也可读 `relax.in` 做 `vc-relax`。先检查弛豫收敛，再从最终晶胞与结构生成 SCF。 |
| `ph.x` | `ph.x -in ph.in > ph.out 2>&1` | 读取同一 `prefix` 和 `outdir` 下的 SCF 数据。Γ 点声子、电子介电张量、Born 有效电荷和可选 Raman 响应写入 `fildyn`。 |
| `dynmat.x` | `dynmat.x -in dynmat.in > dynmat.out 2>&1` | 读取 `fildyn`，施加指定 ASR 并对角化。可对指定方向加非解析项。活性与退偏模表在 stdout，振型另写入文件。主要后处理由单进程执行。 |
| `q2r.x` | `q2r.x -in q2r.in > q2r.out 2>&1` | 将 `ph.x` 完整 q 网格的动力学矩阵（`fildyn` 前缀）变换为实空间力常数 `flfrc`。`zasr` 约束 Born 有效电荷。 |
| `matdyn.x` | `matdyn.x -in matdyn.in > matdyn.out 2>&1` | 读 `flfrc`，按 q 路径或网格插值，`flfrq` 保存频率，`dos=.true.` 可求声子态密度。 |

`q2r.x` 和 `matdyn.x` 属于独立的声子色散或态密度流程，需要 `ph.x` 的完整 q 网格结果。单 Γ 点 `.dyn` 和 `matdyn.x` 的色散频率不能提供这个 skill 所需的 Raman 张量和活性。见 [QE 单 q 点工作流](https://www.quantum-espresso.org/Doc/ph_user_guide/node8.html)、[`q2r.x` 输入](https://www.quantum-espresso.org/Doc/INPUT_Q2R.html) 与 [`matdyn.x` 输入](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html)。

## Shell 启动参数速查

| 参数 | 所属层 | 用法 |
|---|---|---|
| `-in FILE` | QE | 指定输入文件。同义选项有 `-i`、`-inp`、`-input`。MPI 调用优先显式使用 `-in`。 |
| `mpirun -np N` | MPI 启动器 | 设置 MPI 进程总数。集群用调度系统支持的启动器和资源配置。 |
| `OMP_NUM_THREADS` | 环境变量 | 设置每个 MPI 进程的 OpenMP 线程数，实际收益取决于构建和程序。核对 `MPI ranks × 线程数` 与申请的 CPU 数。 |
| `-nk N` / `-npool N` / `-npools N` | QE `pw.x`、`ph.x` | 将 k 点分给 N 个 pool。按每 image 的 MPI 进程数和实际 k 点数选择。只有一个 k 点时用一个 pool。 |
| `-nt N` / `-ntg N` | QE `pw.x` | FFT task group 数。小晶胞用默认配置。按实际 FFT 网格和计时调整。 |
| `-nd N` / `-ndiag N` / `-northo N` | QE `pw.x` | 子空间对角化的处理器组规模，可用时依赖编译的并行线性代数库。显式设置通常取满足 pool 容量的平方数。 |
| `-ni N` / `-nimage N` | QE `ph.x` 进阶 | 在 images 间分派 q 点或不可约扰动。官方指南要求先按每 image 的进程数跑 SCF。分派计算完成后，按该配置用 `recover=.true.` 组装动力学矩阵。普通 Γ Raman 先用单 image。 |

skill 的 `run` / `all` 暴露 `--np`、`--mpi-command`、`--omp-num-threads`、`--bin-dir`，按阶段调用 QE。细分 QE MPI 参数适用于手动调用，不能作为 skill 的同名选项传入。`--skip-existing` 根据阶段记录复用完整产物。`ph.x` 的 `recover` 从 QE scratch 恢复兼容的中间响应。比较不同输入或赝势时可用独立目录，改变并行布局时先核对恢复条件。见 [QE MPI 并行参数](https://www.quantum-espresso.org/Doc/user_guide/node20.html)、[`pw.x` 并行建议](https://www.quantum-espresso.org/Doc/pw_user_guide/node19.html) 与 [`ph.x` 并行和恢复](https://www.quantum-espresso.org/Doc/ph_user_guide/node17.html)。

## `pw.x` 初始结构与 SCF

| 输入区 | 关键项 | 在 Raman 链路中的含义 |
|---|---|---|
| `&CONTROL` | `calculation='scf'`、`prefix`、`outdir`、`pseudo_dir` | `ph.x` 要匹配 `prefix` 和 `outdir`，保留波函数和密度。skill 的 `--outdir` 是项目目录，配置 `qe_outdir` 是 QE scratch。 |
| `&SYSTEM` | `ibrav=0`、`nat`、`ntyp` | 显式写晶胞和原子。核对 `CELL_PARAMETERS angstrom`、`ATOMIC_POSITIONS crystal` 和 `ATOMIC_SPECIES`。 |
| `&SYSTEM` | `ecutwfc`、`ecutrho`，单位 Ry | 波函数和密度/势的截断能。skill 对默认 HGH NC 将 `ecutrho` 下限设为 `4 × ecutwfc`。峰位、活性及介电响应仍需收敛检查。 |
| `&SYSTEM` / `&ELECTRONS` | `occupations='fixed'`、`nbnd`、`diago_full_acc` | 固定占据要求有隙。skill 为中性闭壳层体系加一个空带并提高其对角化精度，以读取最高占据和最低空带能级。报告带隙取自实际采样 k 网格。 |
| `&ELECTRONS` / `K_POINTS automatic` | `conv_thr`、`mixing_beta`、k 网格与偏移 | 分别控制 SCF 误差阈值、混合和 Brillouin 区采样。k 网格与下面 `ph.x` 的声子 q 点是不同量。 |
| `&CONTROL` / `&CELL`（变胞弛豫） | `nstep`、`forc_conv_thr`、`press_conv_thr` | 分别控制步数上限、力阈值（Ry/Bohr）和压力阈值（kbar）。到达步数上限时保留日志，检查 BFGS 是否收敛后再使用末态。 |

见 [`pw.x` 输入手册](https://www.quantum-espresso.org/Doc/INPUT_PW.html)。

## `ph.x` Γ 点响应

| 项 | 默认或 skill 设置 | 用途 |
|---|---|---|
| `prefix`、`outdir`、`fildyn` | 与 SCF 配对，`fildyn='<prefix>.dyn'` | 将动力学矩阵及可用的长程响应交给 `dynmat.x`。 |
| `trans`、`epsil`、`lraman` | Raman 路线均为 `.true.` | 分别启用振动、电场响应和非共振 Raman 系数。`epsil` 与 `trans` 在 Γ 用于 ε∞ 与 Born 有效电荷。`lraman` 也会启用电场响应。 |
| `ldisp`、结尾 q 三元组 | `ldisp=.false.`，namelist `/` 之后为 `0.0 0.0 0.0` | 单 Γ 点。有限 q 或 `ldisp=.true.` 的色散是独立任务。 |
| `tr2_ph`、`niter_ph`、`alpha_mix(niter)` | skill 写 `tr2_ph=1e-14`，其余取 QE 默认 | 分别控制响应收敛阈值、最大迭代数和混合。遇到不收敛先检查 SCF、结构和残差，再调整。 |
| `eth_rps`、`eth_ns`、`dek` | QE 7.5 手册为 `1e-9`、`1e-12`、`1e-3` | 专门控制 Raman 响应中的投影方程阈值、非自洽波函数阈值和波函数对 k 求导的差分步长。skill 当前不自动扫描这些项。 |
| `amass(i)` | skill 按 `ATOMIC_SPECIES` 顺序写质量 | 指定物种质量。研究同位素时核对 `ph.x` 与后处理读入的质量。 |
| `search_sym` | QE 默认 `.true.`。skill 可用 YAML `ph_search_sym: false` | 关闭声子模式对称性分析的可选复验。适用于已核对结构后仍出现 `set_irr_sym_new: wrong representation` 的诊断场景。 |
| `asr` | `ph.x` 的逻辑值属于 Γ 快速算法 | 本机 QE 的 `lraman` 会关闭该快速算法及 `ph.x` 中的这个 ASR 开关。依赖 `dynmat.x` 的 `asr` 控制后处理。 |

见 [`ph.x` 输入手册](https://www.quantum-espresso.org/Doc/INPUT_PH.html)。`ph.x` 的 ASR 与 Raman 兼容细节还核对了本机 `PHonon/PH/phq_readin.f90`。

## `dynmat.x` 方向、ASR 和导出

| `&INPUT` 项 | 用途 |
|---|---|
| `fildyn` | 指向刚生成的 Γ 动力学矩阵。含 ε∞ 与 Born 电荷才能加入三维极性非解析项，含 Raman 张量才能求 Raman 活性。 |
| `q` | 三个笛卡尔分量，用于 Γ 附近传播方向。省略或全零得到未加非解析项的结果。`q = 0.0, 0.0, 1.0` 与 `q(3) = 1.0` 在其余分量为零时等价。晶轴方向需用实际晶胞换算为笛卡尔向量。 |
| `asr` | 这里是字符串。QE 默认为 `'no'`，skill 使用 `'crystal'` 的三平移投影。`'simple'` 是另一种校正方式。对于残余光学虚频仍需排查数值收敛和结构。 |
| `filout` / `fileig` | `filout` 写频率及归一化位移，`fileig` 可另写动力学矩阵的正交本征矢。比较简并或混合模式时保留本征矢。避免把 `filout` 和重定向的 `.out` 命名成同一个文件。 |
| `filmol` / `filxsf` | 导出 Molden 与 AXSF。各 LO 方向使用独立文件名，便于与 TO 振型对照。 |
| `amass(i)` / `lperm` / `lplasma` | 分别覆盖质量、计算模对介电响应的贡献、计算模有效等离子体频率（会启用 `lperm`）。属于手动高级选项。 |

见 [`dynmat.x` 输入手册](https://www.quantum-espresso.org/Doc/INPUT_DYNMAT.html)。在本机 QE 8.0dev 用同一 LiH Γ 动力学矩阵实测，整向量与第三分量写法的 `filout` 逐字节一致，省略 `q` 给出另一组频率。

## 从现成工作目录手动复核

以下命令在 `prepare` 生成的目录中执行。先核对 SCF 和声子是否收敛、日志是否完成、响应与矩阵是否生成，再运行下一行。

```bash
export OMP_NUM_THREADS=1
mpirun -np 4 pw.x -nk 1 -in scf.in > scf.out 2>&1
mpirun -np 4 ph.x -nk 1 -in ph.in > ph.out 2>&1
dynmat.x -in dynmat.in > dynmat.out 2>&1
# 如果 prepare 生成了单个 LO 方向输入
dynmat.x -in dynmat_lo.in > dynmat_lo.out 2>&1
```

优先用 skill 的 `run` 保留阶段签名和续跑状态。手工复验时使用独立目录，保留原始日志、`.dyn` 和振型。解释输出时见 [Raman 理论与判读](raman_theory.md)，遇到异常输入时查 [输入预检与回复指南](input-triage.md)。
