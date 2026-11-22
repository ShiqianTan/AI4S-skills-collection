---
name: qe-mp20-benchmark
description: Reproduce Materials Project MP20 geometry, band gap and formation-energy values with Quantum ESPRESSO. Prepares MP CIFs with SSSP or MPRelaxSet-like settings, runs vc-relax + fixed SCF + dense nscf on GPU or CPU, collects parity metrics, applies MP2020 corrections, and plots gap/volume/formation-energy parity. Use to validate QE against MP/VASP on relaxed structures and thermodynamic screening. Includes k/ecut convergence and SSSP precision sensitivity checks.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# MP20 benchmark (qe-mp20-benchmark)

Validates QE + SSSP against the Materials Project on the MP20 subset:

```
MP CIF -> vc-relax -> (fixed SCF) -> dense nscf band gap
        -> QE elemental references -> raw formation energy -> MP2020 correction
        -> parity plots and summary
```

Reference results: `references/mp20-results.md`.

## Scope

- Periodic, non-magnetic, main-group systems (MP20 test split).
- SSSP 1.3.0 PBE efficiency (default) or MPRelaxSet-like settings.
- Band gap from a **dense uniform nscf** (0.15 Å⁻¹), not a high-symmetry path.
- Formation energy referenced to QE-relaxed elemental ground states, with and
  without MP2020 corrections.
- `e_above_hull` is **not** reproduced (needs the full competing phase diagram).

## Quick start

```bash
# 0) one-time: unpack the bundled SSSP pseudopotentials
python ../pseudos/fetch_pseudos.py

# 1) prepare 20 systems + elemental references (SSSP or MPRelaxSet-like)
python scripts/mp20_prepare.py --config config.yaml --mode sssp
python scripts/mp20_prepare.py --config config.yaml --mode mprelax

# 2) run on GPUs 0-7 (or --gpu-mode cpu)
python scripts/run_mp20.py --config config.yaml --mp20-dir runs/mp20_20 \
  --gpu-mode gpu --devices 0,1,2,3,4,5,6,7

# 3) collect, correct, plot
python scripts/collect_mp20.py --config config.yaml --mp20-dir runs/mp20_20 --tag mp20_20
python scripts/apply_mp2020.py --config config.yaml --tag mp20_20
python scripts/plot_results.py --config config.yaml --tag mp20_20

# Optional sensitivity checks
python scripts/convergence_test.py --config config.yaml --mp20-dir runs/mp20_20
python scripts/precision_test.py --config config.yaml
```

Outputs land in `output.results_dir` (default `runs/results`):
`{tag}_results.csv`, `{tag}_summary.json`, `{tag}_mp2020.csv`,
`{tag}_mp2020_summary.json`, and parity PNGs.

## Configuration

`config.yaml` (see `config.example.yaml`) provides QE binaries, MPI launcher,
dataset roots, and pseudo directories. `--mode sssp` uses the cutoffs already
stored in `data/mp20_20/manifest.json`; `--mode mprelax` recomputes them with
the MPRelaxSet mapping (ENCUT 520 eV floor, reciprocal density 64, EDIFF 5e-5
eV/atom, SIGMA 0.05 eV, NSW 99).

## Bundled data

`data/mp20_20/` contains the 20 CIFs, the elemental reference CIFs, and the MP
reference manifest (gap, formation energy, space group, volume). Pseudos are in
the repository-level `pseudos/` directory.

## Results (from the study)

| Metric | SSSP | MP-like |
|---|---:|---:|
| gap MAE | 0.1686 eV | 0.1680 eV |
| volume MAE | 0.4347 % | 0.4053 % |
| MP2020 Ef MAE | 0.0145 eV/atom | 0.0131 eV/atom |
| raw Ef MAE | 0.1724 eV/atom | 0.1730 eV/atom |

Main outliers: PbS/PbSe (~-0.5 eV), InN (-0.48 eV), AgI (-0.34 eV), GaN
(-0.32 eV). k/ecut and SSSP precision tests move the gap by <0.01 eV, so the
remaining gap disagreement is mostly the MP reference definition, not QE
numerics. See `references/band-gap-caveats.md`.

## Layout

```
qe-mp20-benchmark/
├── SKILL.md
├── README.md
├── requirements.txt
├── config.example.yaml
├── scripts/          # prepare / run / collect / mp2020 / plot / convergence / precision
├── data/mp20_20/     # CIFs + MP reference manifest
└── references/       # results, band-gap caveats, formation-energy notes, workflow
```
