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
