# 示例结构（MP-20 验证体系）

本目录包含 6 个从 Materials Project MP-20 `test.csv` 抽取的小型绝缘体系，
用于 `qe-powder-raman` 的流程回归测试。它们均已在本地 QE 8.0dev（LDA + HGH
模守恒赝势）下完成 `pw.x → ph.x → dynmat.x` 全链路验证。

| 文件 | 体系 | 原子数 | 空间群 | 说明 |
|---|---|---|---|---|
| `LiF_mp-1138.cif` | LiF | 2 | Fm-3m (225) | 岩盐；**负对照**：TO 声学模仅红外活性，无一阶拉曼 |
| `MgH2_mp-23710.cif` | MgH2 | 6 | P4₂/mnm (136) | 金红石；5 支拉曼活性支（最强 1230.80 cm⁻¹） |
| `Al2O3_mp-1143.cif` | Al2O3 | 10 | R-3c (167) | 刚玉；P1 CIF 破缺对称，拉曼支数膨胀（生产需先对称化） |
| `AlAs_mp-2172.cif` | AlAs | 2 | F-43m (216) | 闪锌矿；TO 332.3 / LO 365.0 cm⁻¹，LO-TO 分裂 32.6 cm⁻¹ |
| `CdS_mp-672.cif` | CdS | 4 | P6₃mc (186) | 纤锌矿；E2/A1/E1/E2 齐全，LO-TO 各向异性（q∥c / q⊥c） |
| `MoS2_mp-2815.cif` | MoS2 | 6 | P6₃/mmc (194) | 2H 层状；E2g1 369.2 / A1g 391.2 cm⁻¹，介电强各向异性 |

> **注意**：`CdS_mp-672.cif` 与 `MoS2_mp-2815.cif` 是 MP-20 的**原始** CIF，
> 其六方角 γ 含 ~1e-5° 数值噪声（如 119.99998586°），会触发 `ph.x` 的
> `divide_class` 对称性分析错误或崩溃。`prepare` 会调用 `spglib.refine_cell`
> 自动清洗并写出 `validated.cif`，请始终通过 `prepare` 生成计算输入。

## 一键生成全部输入

```bash
# 在本 skill 目录下
bash examples/run_examples.sh /path/to/output_dir
```

等价于对每个 CIF 执行：

```bash
python3 scripts/qe_raman.py prepare examples/LiF_mp-1138.cif \
    --outdir out/LiF_mp-1138 --prefix LiF
```

需要 LO-TO（闪锌矿/纤锌矿等）时，用 `--lo` 指定方向（每 3 个数一个方向），
或在 config 的 `lo_to_directions` 中列出：

```bash
python3 scripts/qe_raman.py prepare examples/CdS_mp-672.cif \
    --outdir out/CdS --prefix CdS --lo 0 0 1 1 0 0
```

## 端到端计算（可选）

```bash
python3 scripts/qe_raman.py all examples/LiF_mp-1138.cif \
    --outdir out/LiF_mp-1138 --prefix LiF --np 4 --execute
```

`--execute` 会依次执行 `prepare → run → parse → spectrum → report`；
省略 `--execute` 则只生成输入文件。

## 预期结果与容差（回归基准）

来源：本地 QE v8.0dev（commit 215ff8a23，LDA HGH NC，`ecutwfc=60 Ry`，
`ecutrho=240 Ry`，`occupations='fixed'`），完整日志见
`raman_validation_mp20/<system>/`。下表对每个体系给出主峰频率、Raman 活性与
介电（用于回归/对照）。

| 体系 | k 网格 | 主拉曼峰 (cm⁻¹, 活性 Å⁴/amu) | 其它关键量 |
|---|---|---|---|
| LiF | 6×6×6 | 无一级拉曼（负对照） | TO 280.38/280.73/280.73（仅 IR 5.0361）；ε∞=2.2247 |
| MgH2 | 4×4×4 | 1230.80 (A1g, 72.30)；922.90 (Eg, 32.33) | ΣRaman=214.33；ε∞≈4.262 |
| Al2O3 | 4×4×4 | 383.62 (2.11)；571.09 (1.77) | P1 CIF 破缺对称 → 13 支；Σ=9.02；ε∞≈3.446 |
| AlAs | 6×6×6 | 332.33/332.38/332.38 (TO, 354.91/354.77/354.77) | LO 365.00（q=(1,0,0)，分裂 32.62）；ε∞=12.2518 |
| CdS | 6×6×4 | 143.32 (A1, 17.62)；153.31 (E2, 8.78) | LO：q‖c A1→230.27，q⊥c E1→230.61；ε∞≈5.818 |
| MoS2 | 4×4×2 | 369.18 (E2g1, 706.05)；391.16 (A1g, 590.21) | 软模 -17.29（未弛豫）；ε∞=10.4051/10.4051/3.3349；Σ=2004.06 |

容差与注意事项：

- **频率**：在**相同 k 网格与参数**下复现，峰位容差约 ±2 cm⁻¹；改变 k 网格
  （如默认 `kpt_density=30`）会使峰位变化数 cm⁻¹，属正常收敛行为，不应按
  上表逐位比对。
- **Raman 活性**：与 k 网格/展宽设置相关，作**相对**比较（±5% 量级），
  不可当作绝对强度或实验强度。
- **简并**：AlAs TO、CdS E2/E1 等为简并支，表中列出等值多个峰属正常。
- **Al2O3**：示例为 MP-20 的 P1 结构，对称性破缺导致拉曼支数多于理想 R-3c
  （7 支）；生产应先对称化/弛豫。
- **CdS/MoS2**：原始 CIF 有六方角噪声，务必经 `prepare`（内部 `spglib.refine_cell`）。
- **MoS2** 含 -17.29 cm⁻¹ 软模（结构未弛豫），`parse` 会给出 `warning` 状态。
- 实验对照请使用各自文献值，本表仅作**计算回归基准**，不与实验数值混用。

## 收敛扫描示例

```bash
# 只准备输入（不跑 QE）
python scripts/qe_raman.py converge examples/LiF_mp-1138.cif \
    --outdir runs/LiF/conv --groups ecutwfc --ecutwfc 50 60 70

# 真正执行并生成 convergence.{json,csv,md}
python scripts/qe_raman.py converge examples/LiF_mp-1138.cif \
    --outdir runs/LiF/conv --groups ecutwfc,kpt_density,tr2_ph \
    --np 4 --execute
```

判据见 `references/qe_compatibility.md` 的“收敛扫描与判据”。粗网格（如
`ecutwfc 30`）与 60 Ry 基线之间出现数 cm⁻¹ 到十几 cm⁻¹ 的峰位差、以及
ε∞ 的百分比级变化，属于正常的未收敛表现；用该工具给出推荐值，不要直接
把粗网格结果当作生产数据。

## 模式指认与实验对照

```bash
# 先完成 parse；实验峰表列：freq_cm1,intensity[,error]
python scripts/qe_raman.py assign --outdir runs/AlAs --experimental exp_peaks.csv --tol 20

# 也可以给完整实验谱（>50 行）：自动检峰 + Lorentzian 拟合
python scripts/qe_raman.py assign --outdir runs/AlAs --experimental exp_spectrum.csv
python scripts/qe_raman.py report --outdir runs/AlAs   # 第 9 节汇总模式指认
```

MoS₂ 等层状体系在 `prepare`/`preflight` 会给出 `layered_gap` 警告（本例
约 4.3 Å 的 S–S 间隙、2 层），提示先做 vc-relax 并检查层间软模。
