# QE 兼容性与约束

本 skill 针对 **Quantum ESPRESSO 8.x（CPU 构建）** 的固体粉末一级拉曼流程，
核心命令为 `pw.x`、`ph.x`（`&inputph` 带 `lraman=.true.`）与 `dynmat.x`。

## `lraman`（ph.x 拉曼张量）的硬性要求

QE 的 `lraman` 只在满足以下条件时才可用（源码 `PHonon/PH/phq_readin.f90`）：

| 要求 | 说明 |
|---|---|
| **模守恒赝势（NC）** | `lraman` 仅支持 NC 赝势；**超软（US/PAW）一律不支持** |
| **Γ 点（q=0）** | 电场/拉曼响应只在 Γ 点定义，`q = 0 0 0` |
| **绝缘、非磁** | 金属、非共线磁性、Hubbard U、meta-GGA、杂化泛函、TS-vdW 均不支持 |
| **LDA（PZ 族）+ 整数占据** | `occupations='fixed'`，不接受分数占据 |
| **CPU 构建** | GPU 拉曼未实现；本 skill 假设 CPU 版本 |

本 skill 的 `preflight` 会逐项检查，对不满足的体系给出 `FAIL` 并拒绝生成输入。

## dynmat.x 的 LO-TO 方向变量名是 `q`

`dynmat.x` 读取 `&input` namelist 中的 `q`（一个 3 分量实数向量），**不是** `q(3)`：

```fortran
&input
  fildyn = 'AlAs.dyn'
  asr = 'crystal'
  q = 1.0 0.0 0.0     ! LO-TO 沿该方向的声子
  filout = 'dynmat_lo.dat'
/
```

- 不设置 `q` → 输出 TO 声子，并打印 `A direction for q was not specified: TO-LO splitting will be absent`。
- 对各向异性晶体（如纤锌矿 CdS、层状 MoS2），应分别沿 `q∥c` 和 `q⊥c` 各跑一次，**不要**用 x/y/z 三个笛卡尔方向的朴素平均替代粉末 LO-TO。

## 数值噪声与对称性

- MP-20 等数据库导出的六方/三方 CIF 常带 ~1e-5° 的角度噪声（如 γ=119.99998586°），
  会触发 `ph.x` 的 `divide_class` 错误或段错误。`prepare` 内部用
  `spglib.refine_cell`（symprec=1e-3）自动清洗，并写出 `validated.cif`。
- 闪锌矿/岩盐的 **菱面体原胞**（角度≈60°）会被 QE 识别为 C3v 子群（如 AlAs
  仅 6 个对称操作而非 24 个）；结果正确但对称性/成本不最优，生产建议改用
  常规立方晶胞。

## 已验证组合

- QE v8.0dev（commit 215ff8a23），macOS arm64，gfortran 16.2 + OpenMPI + OpenBLAS + FFTW + ScaLAPACK。
- HGH LDA 模守恒赝势（`X.pz-hgh.UPF`，`SLA-PZ-NOGX-NOGC`），已内置 10 个元素：
  Li, F, Mg, H, Al, O, As, Cd, S, Mo。
- 可选真 QE 回归：`tests/integration/test_real_qe.py`（需设置 `QE_REAL_TESTS=1`，
  默认跳过）。以 LiF 负对照（岩盐，无一级拉曼）走通 `all --execute` 主流程，
  用于在本机已验证构建上做端到端冒烟。

## 极性粉末方向平均评估（P2-03）

本节记录“粉末拉曼谱是否需要/如何做取向平均”的评估，**当前不声称已实现
生产级的方向平均**，仅说明现状与可行路线。

**现状**：默认谱来自 Γ 点 TO（`dynmat.x` 未设 `q`，无解析修正）；LO-TO 需由
用户显式给出若干 `q` 方向，分别计算并在报告中按方向并列展示，不做自动合并。
这一取舍避免了在缺乏对称性权重时对 LO-TO 做朴素 x/y/z 平均引入的偏差。

**三种方向平均路线（评估）**：

1. **方位角采样**：在球面上按立体角均匀采样 `q`，逐方向算 LO-TO 后按权重
   （立体角/四面体）求和。最通用，但成本 ≈ 每个方向一次计算；对 `ph.x`
   昂贵的体系不现实。
2. **点群不可约方向 + 多重度权重**：利用晶体点群把 `q` 折叠到不可约楔形，
   按每个不可约方向的多重度加权。立方各向同性 ⇒ 单方向即可；六方/三方需
   `q∥c` 与 `q⊥c` 两类并配权重。这是生产上更现实的做法。
3. **实验散射几何**：粉末无择优取向，理论上等价于各向同性平均；若需与特定
   偏振配置对比，应显式写出散射张量收缩，而非套用单晶几何。

**立方对称性校验**：立方晶体（如 AlAs 闪锌矿）的 LO-TO 应各向同性，任意
等价方向结果应在数值容差内一致（可作方向平均实现的健全性检查）。

**各向异性收敛**：CdS（纤锌矿）、MoS2（2H 层状）等需检验方向数增加时谱形
是否收敛。目前只演示 `q∥c` 与 `q⊥c` 两组，未做多方向收敛记录——若未来实现
自动方向平均，必须补充该收敛证据后才能对外宣称可用。

**结论**：默认活性谱与报告的定位是“沿给定方向的 LO-TO + Γ 点 TO”，
极性粉末的自动方向平均列为后续工作（P2-03 未做生产声明）。
