---
name: qe-matpes-benchmark
description: Validate Quantum ESPRESSO forces, stress and same-material relative energies against MatPES-PBE. Streams the MatPES JSON, selects 10 materials x 4 frames spanning a force range, writes static SCF inputs, runs QE, compares forces/stress/energies with MatPES, and plots parity. MatPES stress is in kbar and needs the kbar-to-eV/A^3 conversion plus the QE sign flip. Use as the primary near-equilibrium force/stress reference.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# MatPES-PBE benchmark (qe-matpes-benchmark)

Validates QE against MatPES-PBE, the best-agreeing force/stress reference in
the study.

```
MatPES-PBE JSON -> select 10 materials x 4 frames (force targets 0.01-2.5 eV/A)
                -> QE static SCF -> force/stress/relative-energy parity
```

Reference results: `references/matpes-results.md`.

## Scope

- 40 frames = 10 materials × 4 frames: InN, CaO, AgF, GeSe, TlCuPd2, TePb,
  LiC6, C, KH3, CaC4.
- Stress: MatPES stores **kbar**; multiply by `6.241509e-4` to eV/Å³ and flip
  the QE sign.
- Same-material relative energies only.

## Quick start

```bash
python ../pseudos/fetch_pseudos.py

python scripts/matpes_scan.py --config config.yaml          # optional group summary
python scripts/matpes_select.py --config config.yaml        # writes jobs.json + scf.in
python scripts/run_mptraj.py --config config.yaml --mptraj-dir runs/matpes \
  --tag matpes --gpu-mode gpu --devices 0,1,2,3,4,5,6,7
python scripts/plot_results.py --config config.yaml --tag matpes
```

Bundled `examples/` has three selected materials (12 frames) with structure,
MatPES reference and SCF input:

```bash
python scripts/run_mptraj.py --config config.yaml --mptraj-dir examples \
  --tag matpes_examples --gpu-mode cpu
python scripts/plot_results.py --config config.yaml --tag matpes_examples
```

## Results (from the study)

| Metric | Value |
|---|---:|
| Force MAE | **0.0123 eV/Å** |
| Force RMSE | 0.0204 eV/Å |
| Max force error | 0.111 eV/Å (InN) |
| Stress MAE | **7.40e-4 eV/Å³** |
| Relative energy | **0.28-3.99 meV/atom** |

MatPES-PBE is the most consistent of the three force/stress references
(MPtrj, sAlex, MatPES).

## Layout

```
qe-matpes-benchmark/
├── SKILL.md
├── README.md
├── requirements.txt
├── config.example.yaml
├── scripts/      # matpes_scan / matpes_select / run_mptraj / plot_results
├── examples/     # 3 materials x 4 frames with MatPES references
└── references/   # results and stress conversion
```
