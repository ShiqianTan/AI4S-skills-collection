# 粉末拉曼理论要点

本 skill 计算的是 **晶体粉末的一级非共振谐波拉曼散射**，即 q≈0 的 Γ 点声子
引起的极化率涨落。`ph.x` 的 `lraman` 输出每个模式的二阶拉曼张量 R，
`dynmat.x` 将其归并为各向同性不变量并给出**粉末平均**的拉曼活性与退偏比。

## 粉末平均不变量

对每个简正模式，拉曼张量的各向同性部分 α 与各向异性部分 β² 定义为：

```
α  = (Rxx + Ryy + Rzz) / 3
β² = (|Rxx−Ryy|² + |Ryy−Rzz|² + |Rzz−Rxx|²
      + 6(|Rxy|² + |Ryz|² + |Rzx|²)) / 2
```

粉末平均的拉曼活性与退偏比（`dynmat.x` 直接输出）：

```
Raman 活性   = 45 α² + 7 β²          (单位 Å⁴/amu)
退偏比 ρ     = 3 β² / (45 α² + 4 β²)   (取值 0 ~ 0.75)
```

- 全对称（A1g 等）模 β²≈0 → ρ≈0（高度偏振）；ρ=0.75 对应纯各向异性（如 E 类简并模）。
- 本 skill 默认光谱 = **活性谱**（直接以 Å⁴/amu 为纵坐标，仅作相对强度比较）。

## Stokes 强度（可选）

物理 Stokes 散射强度需要额外乘上热布居与散射频率因子：

```
I(ω) ∝ A(ω) · (n(ω)+1) · (ω_L − ω)⁴ / ω
```

其中 `n(ω)+1 = 1/(1 − exp(−ℏc ω/kB T))` 为玻色占据，`ω_L = 1e7/λ`（λ=激光波长，nm）。
按波数写作 `ℏc ν̃ / (kB T)`，其中 `ℏc/kB = 1.438776877 cm·K`，ν̃ 单位 cm⁻¹。
`qe_raman.py spectrum --type stokes` 提供该选项；其绝对量纲为任意单位，
只用于与实验谱线的相对形状对比，**不可作绝对强度声明**。

- **展宽约定**：默认按**峰高**缩放（每个模式峰高 = 其活性），Lorentzian/Gaussian 的半高宽
  由 `--fwhm` 给定；FWHM 必须 > 0，否则在写谱前报错。
- **Stokes 窗口**：只保留满足 `0 < ω_L − ω < ω_L` 的模式；超出该窗口的峰会给出诊断并被排除。
- **Bose 指数**：`hc·ν̃/(kB·T)`，温度必须 ≥ 0 K。

## Clausius–Mossotti（局域场）修正

`dynmat.x` 打印两个因子：

- `multiply by <f> for Clausius-Mossotti correction`（极化率/介电）
- `multiply Raman by <f> for Clausius-Mossotti correction`（拉曼）

比较不同体系的拉曼强度时，应将每个模式的活性乘以该因子（或 `f²`）后再比较。
注意：本地 QE 的 `RamanIR` 以 `cmfac**2` 打印 `multiply Raman by` 因子，因此只需
**应用一次**；本 skill 默认输出的是 QE 原始 Raman 列（未乘该因子），报告会给出因子
并在备注中说明如何仅在需要时应用一次。

## 单位速查

| 量 | 单位 | 来源 |
|---|---|---|
| 频率 | cm⁻¹（另给 THz） | dynmat.x |
| IR 活性 | (D/Å)²/amu | dynmat.x |
| 拉曼活性 | Å⁴/amu | dynmat.x |
| 退偏比 | 无量纲 0–0.75 | dynmat.x |
| 介电常数 | 无量纲（3×3 张量） | ph.x |

## 参考

- D. Porezag, M. R. Pederson, Phys. Rev. B 54, 7830 (1996)（lraman 实现）。
- 粉末平均公式见 `LR_Modules/dynmat_sub.f90` 的 `RamanIR` 子程序。
