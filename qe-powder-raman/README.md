# qe-powder-raman

基于 Quantum ESPRESSO 的固体粉末拉曼（一级非共振谐波）计算 skill。

## 目录布局

```
qe-powder-raman/
├── SKILL.md                     # skill 定义（frontmatter + 中文说明）
├── README.md                    # 本文件
├── requirements.txt             # numpy / PyYAML / ase / spglib / matplotlib
├── config.example.yaml          # 可编辑配置模板
├── data/
│   ├── defaults.yaml            # 默认参数（合并进 qe_common.DEFAULT_CONFIG）
│   ├── pseudopotential_map.json # 内置赝势清单（信息性）
│   └── raman_knowledge.json     # 兼容性/单位/拒绝规则知识库
├── references/                  # 兼容性 / 理论 / 常见失败
├── pseudo/                      # 内置 HGH LDA NC 赝势（10 个元素）
├── templates/                   # QE 输入模板（Step 2 填充）
├── scripts/                     # 主程序与模块
│   ├── qe_raman.py              # CLI：preflight / prepare / run / parse / spectrum / report / all
│   ├── qe_common.py             # 常量与默认值
│   ├── qe_struct.py             # 结构读取与预检（ASE + spglib）
│   ├── qe_pp.py                 # UPF 头解析与 NC/LDA 校验
│   ├── qe_inputs.py             # 输入生成（Step 2）
│   ├── qe_run.py                # 运行封装（Step 3）
│   ├── qe_parse.py              # 输出解析（Step 3）
│   ├── qe_spectrum.py           # 光谱绘制（Step 4）
│   ├── qe_report.py             # 报告（Step 4）
│   ├── list_templates.py
│   └── validate_input.py
├── examples/                    # 6 个 MP-20 示例结构（Step 5）
└── tests/                       # 单元/夹具/集成（Step 5）
```

## 安装

```bash
pip install -r requirements.txt
```

需要本机可用的 QE 可执行文件 `pw.x` / `ph.x` / `dynmat.x`（默认读取
`config.example.yaml` 中的 `bin_dir`，可用 `--bin-dir` 或环境变量 `QE_BIN` 覆盖）。

## 使用

```bash
# 预检环境/结构/赝势
python scripts/qe_raman.py preflight <structure> [--pseudo-dir pseudo] [--bin-dir PATH]

# 结构精化 + 生成输入（写 validated.cif / pseudopotential_manifest.json / scf.in / ph.in / dynmat.in / run.sh）
python scripts/qe_raman.py prepare  <structure> --config config.yaml --outdir runs/MySys --prefix MySys [--lo 0 0 1]

# 执行计算（可只跑某阶段 --stage {scf,ph,dynmat,lo}）
python scripts/qe_raman.py run      --outdir runs/MySys --np 4

# 解析（写 raman_modes.csv / summary.json）
python scripts/qe_raman.py parse    --outdir runs/MySys

# 画谱（活性谱默认，--type stokes 给 Stokes；--broadening lorentzian|gaussian）
python scripts/qe_raman.py spectrum --outdir runs/MySys --fwhm 8 [--type both]

# 报告（写 report.md）
python scripts/qe_raman.py report   --outdir runs/MySys

# 一步到位（--execute 触发计算，否则只生成输入）
python scripts/qe_raman.py all <structure> --outdir runs/MySys --np 4 [--execute]
```

输入文件轻量校验：`python scripts/qe_raman.py validate-input runs/MySys/scf.in`。
列出模板：`python scripts/qe_raman.py list-templates`。

## 测试

```bash
python tests/run_tests.py     # 单元 + 集成（无需 QE/MPI）
```

示例结构见 `examples/`（6 个 MP-20 体系，含负对照 LiF 与层状/纤锌矿/闪锌矿等
不同拉曼场景）；`bash examples/run_examples.sh` 可批量生成全部输入。

## 构建路线（6 步，全部完成）

- Step 0 脚手架 + preflight ✅
- Step 1 结构读取 + 赝势映射/校验（`spglib.refine_cell`、`validated.cif`、manifest）✅
- Step 2 输入生成（`relax/scf/ph/dynmat` + LO-TO 可选）✅
- Step 3 运行 + 解析（`raman_modes.csv`、`summary.json`、声学/软模/ASR/LO）✅
- Step 4 光谱 + 报告（Raman 活性谱 + 可选 Stokes，`report.md`）✅
- Step 5 测试 + 夹具 + 示例 + 完整文档 ✅
