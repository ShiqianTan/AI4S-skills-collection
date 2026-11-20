# 常见失败模式与处理

按出现频率排序，来自 6 个 MP-20 体系（LiF/MgH2/Al2O3/AlAs/CdS/MoS2）的实测。

## 1. 六方/三方角噪声 → `divide_class` 错误 / 段错误

**症状**：`ph.x` 读入波函数后立即报
`Error in routine divide_class (1): Wrong classes for C_6v`（CdS），
或 `SIGSEGV/SIGBUS` 崩溃（MoS2）。

**原因**：MP-20 等数据库的 CIF 六方角 γ 含 ~1e-5° 噪声（119.99998586°），
旋转矩阵带 ~2e-7 偏差，超过 `divide_class` 的共轭容差 1e-7。

**处理**：`spglib.refine_cell(symprec=1e-3)` 将 γ 精确到 120.0°。本 skill 的
`prepare` 已自动执行并写出 `validated.cif`；手动时使用 `*_sym.cif`。

## 2. P1 CIF 破缺对称 → 拉曼支数膨胀 / ph.x 代价暴涨

**症状**：Al2O3 刚玉 CIF 为 P1，spglib 只识别 2 个对称操作，30 个模式产生
13 支“拉曼”/11 支“红外”（理想 R-3c 应为 7 拉曼 + 4 红外），且 ph.x 需算 30 个
微扰（~20 分钟）。

**处理**：生产计算前先对称化（`refine_cell` 或 pymatgen 的 `get_primitive_structure`
/ `SpacegroupAnalyzer`），或做 vc-relax。清洗只修数值噪声，**不会**强行升对称。

## 3. 菱面体原胞只识别子群

**症状**：AlAs 闪锌矿用菱面体原胞（角度≈60°），QE 只报 6 个对称操作（C3v）
而非 F-43m 的 24 个。

**处理**：结果正确但非最优；生产建议改用常规立方晶胞（4 原子）以充分利用对称。

## 4. 软模 / 负频率

**症状**：MoS2 未弛豫结构出现 −17.29 cm⁻¹ 的软光学模（虚频）。

**处理**：`parse` 会把含软模的体系标记为 `warning`。先做结构弛豫
（`prepare --relax`）再算拉曼；软模通常表示该几何不稳定。

## 5. SCF 不收敛（`occupations='fixed'`）

**症状**：固定占据下 SCF 迭代数很大或达不到 `conv_thr`。

**处理**：核对体系确为绝缘且有带隙；必要时提高 `mixing_beta` 或改用 `mixing_mode='local-TF'`。
金属/窄带隙体系超出本 skill 的 MVP 范围。

## 6. 赝势不兼容

**症状**：`preflight` 报 `FAIL`：`is PAW / is ultrasoft / functional not LDA`。

**处理**：`lraman` 只支持 **NC + LDA(PZ)**。更换为 HGH LDA 模守恒赝势
（内置 10 元素，或 `--pseudo-dir` 指定外部 NC/LDA 集合）。

## 7. 磁性元素提示

**症状**：`preflight` 对含 Mo、Fe、Co、Ni 等元素报 `WARN magnetic_hint`。

**处理**：MVP 假设非磁。磁性元素提示需结合已知磁序或独立电子结构结果核对。必须显式描述磁序时，选择已验证支持该磁态的 Raman 方法。

## 8. LO-TO 缺失 / 变量名写错

**症状**：`dynmat.out` 打印 `A direction for q was not specified: TO-LO splitting will be absent`。

**处理**：给 `dynmat.x` 的 `q` 设置非零笛卡尔向量。例如 `q = 0, 0, 1`，或只设置第三分量 `q(3) = 1`。各向异性晶体可沿实际晶胞的 `q∥c` 与 `q⊥c` 分别对照，并保留各方向振型。见 [QE 子程序与调用参数](qe-command-line.md)。
