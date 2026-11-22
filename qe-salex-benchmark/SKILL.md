---
name: qe-salex-benchmark
description: Validate Quantum ESPRESSO forces, stress and same-material relative energies against the sAlex validation set. Scans sAlex LMDB files, selects small-cell non-magnetic materials and frames, writes static SCF inputs, runs QE, compares with sAlex/ASE references, and plots force/stress/relative-energy parity. sAlex stress already follows the ASE convention, so no kbar conversion or sign flip is applied.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# sAlex validation (qe-salex-benchmark)

Validates QE against the sAlex validation set (LMDB, ASE-style energies, forces
and stress).

```
sAlex val/*.aselmdb -> scan materials -> select frames
                    -> QE static SCF -> force/stress/relative-energy parity
```

Reference results: `references/salex-results.md`.

## Scope

- Small cells (≤8 atoms), SSSP-covered, non-magnetic materials.
- Stress uses the **ASE convention already**; do not apply the MPtrj kbar
  conversion or the QE sign flip.
- Same-material relative energies only; absolute energies differ by
  pseudopotential offsets.

## Quick start

```bash
python ../pseudos/fetch_pseudos.py

# 1) scan LMDBs -> groups.json
python scripts/salex_scan.py --config config.yaml

# 2) select 20 materials x 4 frames, write jobs.json + scf.in
python scripts/salex_select.py --config config.yaml

# 3) run QE (GPU or CPU)
python scripts/run_mptraj.py --config config.yaml --mptraj-dir runs/salex \
  --tag salex --gpu-mode gpu --devices 0,1,2,3,4,5,6,7

# 4) plots
python scripts/plot_results.py --config config.yaml --tag salex
```

Bundled `examples/` contains three already-selected materials (12 frames) with
structure, sAlex reference and SCF input, so the run/plot path can be tested
without the LMDB files:

```bash
python scripts/run_mptraj.py --config config.yaml --mptraj-dir examples \
  --tag salex_examples --gpu-mode cpu
python scripts/plot_results.py --config config.yaml --tag salex_examples
```

## Results (from the study)

80 frames = 20 materials × 4 frames.

| Metric | 40 frames | 80 frames |
|---|---:|---:|
| Force MAE | 0.0234 eV/Å | **0.0223 eV/Å** |
| Force RMSE | 0.0460 eV/Å | 0.0462 eV/Å |
| Max force error | 0.309 eV/Å | 0.645 eV/Å (K2PS5) |
| Stress MAE | 8.93e-4 eV/Å³ | **8.94e-4 eV/Å³** |
| Relative energy | 1.1-14.9 meV/atom | **0.79-19.71 meV/atom** |

sAlex force consistency on distorted frames is better than the MPTraj-20
high-force group.

## Layout

```
qe-salex-benchmark/
├── SKILL.md
├── README.md
├── requirements.txt
├── config.example.yaml
├── scripts/      # scan / select / run_mptraj / plot_results
├── examples/     # 3 materials x 4 frames with sAlex references
└── references/   # results, workflow, stress conventions
```
