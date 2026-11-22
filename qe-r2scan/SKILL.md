---
name: qe-r2scan
description: Run the MatPES-R2SCAN benchmark with Quantum ESPRESSO using libxc and norm-conserving PseudoDojo pseudopotentials on the CPU build. Extracts the same 40 MatPES frames, prepares R2SCAN inputs (input_dft='r2scan', 60/240 Ry), runs SCF, compares forces/stress/relative energies with the MatPES-R2SCAN references, benchmarks MPI x OpenMP/cutoff/k-grid choices, and documents why GPU R2SCAN is unavailable. Use when a meta-GGA cross-check is needed and to quantify how far R2SCAN is from quantitative agreement.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# R2SCAN benchmark (qe-r2scan)

QE's R2SCAN path needs **libxc** and **norm-conserving** pseudopotentials. The
NVHPC GPU build cannot load gfortran libxc modules, so R2SCAN runs on the CPU
libxc build.

```
MatPES-R2SCAN JSON (same 40 frames as MatPES-PBE)
  -> QE SCF: input_dft='r2scan', PseudoDojo NC-SR, 60/240 Ry
  -> force/stress/relative-energy parity vs MatPES-R2SCAN
  -> CPU configuration benchmark
```

Reference results: `references/r2scan-results.md`.

## Quick start

```bash
# one-time: build the libxc CPU install (see qe-build) and unpack pseudos
python ../pseudos/fetch_pseudos.py

# 1) select the same MatPES frames and write R2SCAN inputs
python scripts/matpes_r2scan_select.py --config config.yaml \
  --pbe-jobs ../qe-matpes-benchmark/examples/jobs.json --out runs/matpes_r2scan

# 2) run with the libxc build (CPU, 4 MPI x 4 OpenMP recommended)
python scripts/run_r2scan.py --config config.yaml --r2scan-dir runs/matpes_r2scan \
  --np 4 --concurrent 16

# 3) plots
python scripts/plot_results.py --config config.yaml --tag matpes_r2scan

# optional: benchmark CPU configurations on one frame
python scripts/r2scan_bench.py --config config.yaml \
  --base-input runs/matpes_r2scan/<group>/<sid>/scf.in
```

Bundled `examples/` has three materials (12 frames) with R2SCAN structures,
references and inputs:

```bash
python scripts/run_r2scan.py --config config.yaml --r2scan-dir examples --np 1 --concurrent 1
python scripts/plot_results.py --config config.yaml --tag matpes_r2scan
```

## Results (from the study)

Same 40 `matpes_id` frames, PBE vs R2SCAN:

| Metric | PBE | R2SCAN |
|---|---:|---:|
| Force MAE | **0.0123 eV/Å** | 0.0489 eV/Å |
| Force RMSE | 0.0204 eV/Å | 0.0916 eV/Å |
| Max force error | 0.111 eV/Å | 0.726 eV/Å |
| Stress MAE | **7.40e-4 eV/Å³** | 3.59e-2 eV/Å³ |
| Relative energy | 0.28-3.99 meV/atom | 1.5-48 meV/atom |

R2SCAN is significantly worse than PBE because of PBE-generated NC
pseudopotentials, different VASP PAW R2SCAN settings, unconverged 60 Ry/7×7×7,
and meta-GGA implementation differences. Treat R2SCAN as a qualitative
cross-check, not a quantitative reference.

## Why GPU R2SCAN is unavailable

`nvfortran` cannot read gfortran-compiled libxc Fortran modules; the GPU binary
reports `libxc needed for this functional`. R2SCAN therefore uses the CPU libxc
build: 16 concurrent jobs x 4 MPI x 4 OpenMP = 256 threads.

## Layout

```
qe-r2scan/
├── SKILL.md
├── README.md
├── requirements.txt
├── config.example.yaml
├── scripts/      # matpes_r2scan_select / run_r2scan / r2scan_bench / plot_results
├── examples/     # 3 materials x 4 frames with R2SCAN references
└── references/   # results, limitations, workflow
```
