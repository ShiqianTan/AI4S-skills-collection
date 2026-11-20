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
