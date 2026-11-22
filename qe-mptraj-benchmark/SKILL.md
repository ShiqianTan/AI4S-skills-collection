---
name: qe-mptraj-benchmark
description: Validate Quantum ESPRESSO forces, stress and same-material relative energies against MPtrj VASP trajectories. Streams MPtrj to select lowest-energy, highest-energy and largest-force frames, prepares static SCF inputs (SSSP or MPRelaxSet-like), runs on GPU/CPU, compares forces/stress/energies with VASP, and plots force parity and force-error-vs-distortion. Use for near-equilibrium force/stress validation of MLIPs or generated structures, and to see where high-distortion frames break down.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# MPtrj force/stress benchmark (qe-mptraj-benchmark)

Runs QE static SCF at exact MPtrj VASP frame geometries and compares:

```
MPtrj JSON -> select 3 frames/material (lowest E, highest E, max |F|)
           -> QE static SCF (tprnfor, tstress)
           -> force / stress / relative-energy parity
```

Reference results: `references/mptraj-results.md`.

## Scope and known limits

- Forces and stress are offset-free and compared directly (after sign/unit
  conversion for stress).
- Absolute total energies differ by pseudopotential zeros; only **same-material
  relative energies** are compared.
- Near-equilibrium frames are reliable. High-distortion frames can fail badly;
  see `references/force-error-vs-distortion.md`.

## Quick start

```bash
python ../pseudos/fetch_pseudos.py

# 1) select frames (streams MPtrj_2022.9_full.json)
python scripts/mptraj_select.py --config config.yaml --preset 20

# 2) prepare static SCF inputs + jobs.json
python scripts/mptraj_prepare.py --config config.yaml --mptraj-dir runs/mptraj_20 --mode sssp

# 3) run on GPUs (or --gpu-mode cpu)
python scripts/run_mptraj.py --config config.yaml --mptraj-dir runs/mptraj_20 \
  --tag mptraj20 --gpu-mode gpu --devices 0,1,2,3,4,5,6,7

# 4) plots
python scripts/plot_results.py --config config.yaml --tag mptraj20
```

Bundled `examples/` contains four already-selected materials (AlAs, AgI, MgTe,
InBr; 3 frames each) so the prepare/run/plot path can be exercised without the
full MPtrj file:

```bash
python scripts/mptraj_prepare.py --config config.yaml --mptraj-dir examples --mode sssp
python scripts/run_mptraj.py --config config.yaml --mptraj-dir examples --tag mptraj_examples --gpu-mode cpu
python scripts/plot_results.py --config config.yaml --tag mptraj_examples
```

## Frame selection

For each material: lowest `energy_per_atom`, highest `energy_per_atom`, and the
frame with the largest force. A median-energy frame is used only if fewer than
three unique frames exist.

## Results (from the study)

| Set | Frames | Force MAE | Force RMSE | Max error |
|---|---:|---:|---:|---:|
| MPTraj-5 | 5×3 | **0.0075 eV/Å** | 0.0283 eV/Å | 0.184 eV/Å |
| MPTraj-20 | 20×3 | 0.1428 eV/Å | 0.5171 eV/Å | 6.37 eV/Å |

Force error vs. VASP max |F| (MPTraj-20):

| VASP max |F| | frames | MAE | RMSE | max |
|---|---:|---:|---:|---:|
| <0.2 eV/Å | 31 | **0.031** | 0.069 | 0.428 |
| 0.2-0.5 | 12 | 0.075 | 0.160 | 0.747 |
| 0.5-1.0 | 9 | 0.215 | 0.493 | 1.975 |
| >1.5 | 7 | **0.682** | 1.384 | 6.372 |

The MPRelaxSet-like configuration does not fix high-force frames (MAE 0.145 vs
0.143 eV/Å). See `qe-convergence-attribution` for the k/ecut/smearing/pseudo
factor scan.

## Layout

```
qe-mptraj-benchmark/
├── SKILL.md
├── README.md
├── requirements.txt
├── config.example.yaml
├── scripts/          # find_candidates / select / prepare / run / plot
├── examples/         # 4 materials x 3 frames with VASP references
└── references/       # results, force-vs-distortion, stress conventions
```
