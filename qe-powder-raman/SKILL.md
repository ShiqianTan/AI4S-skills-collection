---
name: qe-powder-raman
description: Compute solid-state (crystal) powder Raman spectra with Quantum ESPRESSO. Run the pw.x SCF -> ph.x (epsil + lraman, Gamma) -> dynmat.x workflow, then broaden modes into a powder Raman spectrum and a report. Use when the user wants a first-order, non-resonant, harmonic powder Raman spectrum for a periodic insulating crystal. Handles structure reading (CIF/POSCAR), pseudopotential validation (NC + LDA only), input generation, execution, output parsing, spectrum plotting, and reporting.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# Quantum ESPRESSO 固体粉末拉曼（qe-powder-raman）

本 skill 提供晶体粉末拉曼（一级、非共振、谐波）的端到端工作流：

```
pw.x (SCF) -> ph.x (epsil + trans + lraman, q=Γ) -> dynmat.x (IR / Raman / 去极化) -> 展宽成粉末谱 -> 报告
```

支持：结构读取（CIF/POSCAR/QE）、赝势校验（只允许 NC + LDA）、输入生成、
可选执行、输出解析、光谱绘制（Raman 活性谱 + 可选 Stokes 谱）、报告汇总。

## 适用与限制（MVP）

- 仅适用于 **三维周期性、绝缘（有带隙）、非磁性、整数占据** 的晶体。
- `lraman=.true.` **只支持模守恒（NC）赝势**，不支持 US-PP 和 PAW。
- 只支持 **LDA（PZ）** 泛函 + `occupations='fixed'`。
- 不支持：非共线磁性、Hubbard U、meta-GGA/杂化泛函、TS-vdW。
- 只算 **Γ 点** 一级非共振谐波拉曼；无温度/峰宽/绝对强度标定（除非显式开启 Stokes）。

## 快速开始

```bash
# 0) 环境与结构/赝势预检
python scripts/qe_raman.py preflight examples/LiF_mp-1138.cif --pseudo-dir pseudo

# 1) 结构精化 + 赝势清单 + 生成输入（scf/ph/dynmat + 可选 LO + run.sh）
python scripts/qe_raman.py prepare examples/LiF_mp-1138.cif --config config.yaml --outdir runs/LiF

# 2) 运行（pw.x -> ph.x -> dynmat.x -> LO） / 解析 / 画谱 / 报告
python scripts/qe_raman.py run      --outdir runs/LiF --stage all --np 4
python scripts/qe_raman.py parse    --outdir runs/LiF
python scripts/qe_raman.py spectrum --outdir runs/LiF --fwhm 8
python scripts/qe_raman.py report   --outdir runs/LiF

# 3) 一步到位（可选 --execute 真正计算，否则只生成输入）
python scripts/qe_raman.py all examples/LiF_mp-1138.cif --outdir runs/LiF --np 4 --execute
```

## 目录结构

```
qe-powder-raman/
├── SKILL.md                 # 本文件
├── README.md
├── requirements.txt
├── config.example.yaml      # 可编辑配置模板
├── data/                    # defaults.yaml / pseudopotential_map.json / raman_knowledge.json
├── references/              # qe_compatibility.md / raman_theory.md / failure_modes.md
├── pseudo/                  # 内置 HGH LDA NC 赝势（10 个元素）
├── templates/               # QE 输入模板（relax/scf/ph/dynmat）
├── scripts/                 # qe_raman.py 及模块
└── examples/                # 6 个 MP-20 示例结构
```

## 关键约定

- `dynmat.x` 的 LO-TO 方向变量名是 **`q`（3 分量）**，不是 `q(3)`。
  默认 `q=(0,0,0)` 只出 TO；设置 `q = x y z` 才加入非解析项（需 `ph.x` 已输出 Z* 和 ε∞）。
- 六方/三方 CIF 若角度有 ~1e-5° 数值噪声，会导致 `ph.x` 对称性分析报错或崩溃，
  必须先 `spglib.refine_cell` 精确化后再生成输入。
- 本 skill 的默认参数复现已验证的 LDA+HGH 流程；生产计算请自行收敛 `ecutwfc`、k 网格并（建议）先弛豫结构。

## 状态

六个构建步骤已全部完成。全部 10 个子命令可用：

| 子命令 | 作用 | 落地步骤 |
|---|---|---|
| `preflight` | 环境/结构/赝势预检 | Step 0 |
| `structure` | 对称精化 + 赝势清单（写 validated.cif / manifest） | Step 1 |
| `prepare` | 生成 scf/ph/dynmat(±LO) 输入 + run.sh | Step 2 |
| `run` | 执行 pw.x/ph.x/dynmat.x（state.json） | Step 3 |
| `parse` | 解析输出（raman_modes.csv / summary.json） | Step 3 |
| `spectrum` | 活性谱 + 可选 Stokes 谱（CSV/PNG/SVG） | Step 4 |
| `report` | 汇总 report.md | Step 4 |
| `all` | 串联以上全流程（`--execute` 触发计算） | Step 5 |
| `validate-input` | 输入文件轻量校验 | Step 2 |
| `list-templates` | 列出模板 | Step 0 |

单元/集成测试见 `tests/`（`python tests/run_tests.py`，无需 QE/MPI）；
示例结构见 `examples/`；兼容性/理论/失败模式见 `references/`。
