# 固体粉末 Raman 理论与结果判读

本 skill 计算三维周期绝缘晶体的 **Γ 点一级、非共振、谐波 Raman**。谐波动力学矩阵给出频率和本征矢，电子极化响应对原子位移的导数给出 Raman 张量。QE 的 `lraman` 用均匀电场的二阶电子响应求这些系数。`dynmat.x` 再把原子位移响应投影到正常模上。见 [Lazzeri 与 Mauri 原始论文](https://doi.org/10.1103/PhysRevLett.90.036401)、[QE `ph.x` 输入手册](https://www.quantum-espresso.org/Doc/INPUT_PH.html)。

## 模式、选择定则与单位

设实际计算晶胞有 N 个原子，Γ 点有 3N 支模式，包括三个平移声学模。Raman 活性来自模坐标改变的电子极化响应。IR 活性来自 Born 有效电荷在振型上的投影。对具有反演中心的晶体，一个具有确定宇称的电偶极模式不能同时满足一级 Raman 和 IR 的活性条件。完整的空间群和实验几何还会限制具体模式的可见性。QE 的 IR 与 Raman 活性按模输出，简并模式仍可占多行。见 [Porezag 与 Pederson 的活性公式](https://doi.org/10.1103/PhysRevB.54.7830)、本机 QE `LR_Modules/dynmat_sub.f90` 的 `RamanIR`。

| 量 | 单位与解释 |
|---|---|
| 模式频率 | `dynmat.x` 输出 cm⁻¹ 与 THz。负号表示动力学矩阵本征值为负，即谐波虚频。 |
| IR 活性 | (D/Å)²/amu，属于 IR 自身的模活性。 |
| Raman 活性 | Å⁴/amu，来自粉末取向平均所用的张量不变量。 |
| 退偏比 | 无量纲。对非零活性，理论范围为 0–0.75。活性为零时分母也为零，退偏比没有物理意义。 |
| 电子介电张量 ε∞ | 无量纲 3 × 3 张量。用于非解析修正和 QE 的局域场约定。 |

对对称的实 Raman 张量 R，各向同性项 α、各向异性项 β²、活性 A 和退偏比 ρ 为

```text
α  = (Rxx + Ryy + Rzz) / 3
β² = [(Rxx−Ryy)² + (Ryy−Rzz)² + (Rzz−Rxx)² + 6(Rxy² + Ryz² + Rzx²)] / 2
A  = 45α² + 7β²
ρ  = 3β² / (45α² + 4β²)
```

近各向同性的张量给出较低的退偏比。全对称模式也可有各向异性，因此仅凭 ρ 不能确定不可约表示。简并簇的谱强度涉及各成员的贡献，单一成员的张量还可随简并基的选择而改变。见 [Porezag 与 Pederson](https://doi.org/10.1103/PhysRevB.54.7830)、本机 QE `RamanIR`。

### Raman 列缺失或空谱

1. 核对 `ph.in` 的 `lraman`、`ph.out` 的收敛与张量计算证据，以及 `.dyn` 中的 Raman 响应。
2. 核对 `dynmat.out` 是否从正确矩阵读取响应。本机 QE 在所有模式的 Raman 张量分量低于内部阈值时会省略 Raman 列。已验证的岩盐 LiF 是一级 Raman 禁戒的回归对照。
3. 若模表已给出活性，检查它是否因输出舍入、阈值、近零频、虚频或 Stokes 窗口被筛掉。skill 默认绘谱阈值为 `1e-4`，且要求频率大于 `1 cm⁻¹`。

`lraman=.false.` 属于声子/IR 路线。该设置或失败日志所产生的空谱不能用来判断选择定则。见 [QE `dynmat.x` 输入手册](https://www.quantum-espresso.org/Doc/INPUT_DYNMAT.html)、本机 `LR_Modules/dynmat_sub.f90` 与 [示例验证](../examples/README.md)。

## Γ 点 TO、方向性 LO 与声学求和规则

- 未设 `dynmat.x` 的 `q` 或使用全零向量时，skill 输出未加极性非解析项的 Γ 结果，简称 TO。对具有长程偶极响应的晶体，ε∞ 与 Born 有效电荷共同决定沿给定笛卡尔传播方向的非解析项。
- 非零 `q` 后仍输出完整的 3N 模表。各向异性晶体中，非解析项可以改变多支模式的频率与本征矢。按模序号相减容易在排序变化或简并处分错。可用 `fileig` 导出正交本征矢，比较模式重叠或简并子空间。逐方向 LO 谱也没有完成对极性粉末的传播方向积分。
- 振型的三条平移声学支理想上在 Γ 为零。skill 使用 `dynmat.x asr='crystal'` 的三平移投影。本机 QE Raman 路线关闭 `ph.x` 的 Γ 快速算法及其布尔 `asr`，最终需检查 `dynmat.x` 的声学残差。当前解析器以三个最小绝对频率近似识别声学模。近零的光学支或明显失真的平移模需要结合振型人工核实。

见 [QE 单 q 点声子指南](https://www.quantum-espresso.org/Doc/ph_user_guide/node8.html)、[QE `dynmat.x` 输入手册](https://www.quantum-espresso.org/Doc/INPUT_DYNMAT.html)、[命令与导出参数](qe-command-line.md)。

## 活性谱、Stokes 谱和局域场

默认 `spectrum_activity` 以原始 Raman 活性为各模式的峰高。`spectrum_stokes` 再乘散射波数和热布居因子。用 ν̃ 表示 **cm⁻¹ 波数**，λ 用 nm 表示激光波长。

```text
ν̃L = 10^7 / λ
n(T) + 1 = 1 / [1 − exp(−hc ν̃ / (kB T))]
I_Stokes,m ∝ A_m (ν̃L − ν̃m)^4 (n(T) + 1) / ν̃m
hc / kB ≈ 1.438776877 cm·K
```

- `T=0` 时 `n+1=1`。保留 `0 < ν̃m < ν̃L` 的 Stokes 模式。波数形式用 `hc`，角频率形式用 `ℏω`。
- Gaussian 和 Lorentzian 都以峰高缩放，FWHM 用 cm⁻¹。改变 FWHM 会改变峰面积与多峰重叠。这个参数用于模拟线形。激光波长和温度用于 Stokes 强度权重。默认流程没有绝对强度、共振增强、非谐寿命或热频移模型。
- `spectrum_ir` 直接展宽 IR 活性，IR 与 Raman 的双轴图可用于对照峰位。比较实验相对强度时记录激光、温度、归一化、样品几何和展宽。`spectrum_manifest.json` 保存实际采用的设置和产物。

QE 可输出两种 Clausius–Mossotti 修正因子。令 `f = 3 / (2 + trace(ε∞)/3)`。`multiply by` 行给出极化率因子 f，`multiply Raman by` 行已经给出活性因子 f²。本 skill 的输出使用 QE 未修正活性。需要该迹平均局域场近似时，对活性乘 `multiply Raman by` 一行的数值一次。该近似应按所比较的模型一致采用，不能代替各向异性散射几何与实验校准。见本机 QE `LR_Modules/dynmat_sub.f90` 的 `RamanIR`。

## 计算验收与异常诊断

| 检查 | 判读 |
|---|---|
| 结构和 SCF | 初始晶相、原子占位与赝势有据可查。启用弛豫时须核对 BFGS 收敛和最终晶胞。检查固定占据、SCF 收敛和最低空带。skill 的 0.05 eV 阈值是采样 k 网格上的停算判据。 |
| DFPT | 核对 Γ 点、电子介电张量、Born 有效电荷、Raman 响应、收敛和 `JOB DONE.`。旧日志缺少关键信息时标明未验证。 |
| 后处理 | TO 和各请求 LO 方向应符合实际晶胞的 3N 支模式。核查连续模式编号、声学残差、张量证据、频率与单位。保留逐方向矩阵和振型。 |
| 光学虚频 | 先复核平移模式和 ASR，再检查几何、k 网格、截断能和响应阈值。显著虚频可由结构或数值问题引起，应保留原始频率和模式编号。第 5 轮 LiH 粗参数出现约 −1030 cm⁻¹，固定晶胞加密至 6 × 6 × 6 并用 60 Ry 后，三支光学模约 +481 cm⁻¹。该对照同时改变了两项设置。 |
| 收敛和模式比对 | 独立变化 `ecutwfc`、电子 k 网格和 `tr2_ph`。比较主要活性峰、简并簇和 ε∞。只有模式排序稳定时，按编号比较才可靠。活性阈值与简并容差属于分析参数。 |

技术边界见 [QE 兼容性](qe_compatibility.md)，按用户输入选择回复见 [输入预检与回复指南](input-triage.md)，供 skill 查询的结构化数据见 [Raman 知识库](../data/raman_knowledge.json)。
