---
name: qe-convergence-attribution
description: Attribute QE-vs-VASP band-gap and high-force disagreements to controllable factors. For band gaps, fixes the structure and compares high-symmetry k-path vs dense uniform grids, SSSP efficiency/precision and PseudoDojo pseudopotentials, initial vs final structures, and PBE/PBEsol/R2SCAN. For forces, runs a one-factor-at-a-time k/ecut/smearing/pseudopotential matrix on the worst MPtrj frames. Use to decide whether a disagreement is numerical, pseudopotential, structural/reference, or a real electronic/PES difference.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# Convergence and error attribution (qe-convergence-attribution)

This skill answers "why do QE and VASP disagree?" with controlled experiments,
instead of adding more systems.

```
Gap:   fixed structure -> path vs dense k, pseudo sets, initial vs final, functional
Force: worst MPtrj frames -> k / ecut / smearing / pseudopotential matrix
```

Reference results: `references/bandgap-controlled.md` and
`references/high-force-matrix.md`.

## Part 1: controlled band gap

```bash
# dense uniform grid (recommended), three pseudopotential sets
python scripts/gap_controlled.py --config config.yaml \
  --wbm-dir ../qe-wbm-basin/examples --mode dense \
  --pseudo-sets sssp_eff,sssp_prec,pdojo --tag gap_dense

# high-symmetry path (shows off-path CBM error)
python scripts/gap_controlled.py --config config.yaml \
  --wbm-dir ../qe-wbm-basin/examples --mode path \
  --pseudo-sets sssp_eff,sssp_prec,pdojo --tag gap_path

# initial vs final structure, or PBEsol/R2SCAN
python scripts/gap_controlled.py --config config.yaml --wbm-dir runs/wbm_gap \
  --mode dense --start initial --tag gap_initial
python scripts/gap_controlled.py --config config.yaml --wbm-dir runs/wbm_gap \
  --mode dense --pseudo-sets pdojo --functional r2scan --tag gap_r2scan
```

- `--nscf-density 0.08` gives the most converged dense grid.
- `r2scan` automatically uses the libxc build and requires PseudoDojo
  norm-conserving pseudos.
- Reference gaps come from each system's `meta.json` / manifest
  (`vasp_bandgap_pbe_eV`).

## Part 2: high-force matrix

```bash
# select the 5 worst-force frames from a prior MPtrj run
python scripts/force_matrix.py --config config.yaml \
  --mptraj-dir runs/mptraj_20 --top-n 5 --tag force_matrix

# or name the frames explicitly
python scripts/force_matrix.py --config config.yaml \
  --mptraj-dir runs/mptraj_20 --sids sid1,sid2 --tag force_matrix
```

Variants: baseline, k 0.15/0.20/0.25/0.30 Å⁻¹, ecut 0.9/1.2/1.5×, MV vs
Gaussian smearing, tetrahedra, SSSP precision, PseudoDojo, and a combined
prec+k0.20+ecut1.2 case.

## Results (from the study)

Controlled gap on 20 WBM semiconductors, fixed VASP-final structure:

| Pseudopotential | gap MAE | RMSE | bias | max |
|---|---:|---:|---:|---:|
| SSSP efficiency | 0.1317 eV | 0.1798 | -0.0202 | 0.4088 |
| SSSP precision | 0.1296 eV | 0.1794 | -0.0181 | 0.4089 |
| PseudoDojo NC-SR | 0.1305 eV | 0.1777 | -0.0181 | 0.4069 |
| **all** | **0.1306 eV** | 0.1790 | -0.0188 | 0.4089 |

Pseudopotential differences are <0.002 eV. Dense uniform grids halve the error:

| Method | gap MAE | RMSE | max |
|---|---:|---:|---:|
| high-symmetry path | 0.1317 eV | 0.1798 | 0.4088 |
| dense 0.10 Å⁻¹ | **0.0665 eV** | 0.0935 | 0.2617 |
| dense 0.08 Å⁻¹ | **0.0764 eV** | 0.1126 | 0.2803 |

Initial-vs-final at dense 0.08: final MAE **0.076 eV**, initial MAE **0.614 eV**
(18/20 finals closer; Al2Ga2Se4 and Al2P1S4 are exceptions).

High-force matrix (5 worst frames, 65 tasks): baseline force MAE **0.8072
eV/Å**; the best single factor changes it by only **0.003-0.021 eV/Å** (0.4-2.6
%), and the maximum error stays at **6.5-6.8 eV/Å**. Tetrahedra cannot be used
for forces.

## Interpretation rules

- If changing the pseudopotential leaves the gap unchanged, the residual is not
  a pseudopotential error.
- If the initial structure matches VASP better than the final, check which
  structure the reference actually corresponds to.
- If k/ecut/smearing/pseudo cannot reduce a high-force error, it is an
  electron-state / local-minimum / PES difference, not a convergence problem.
- Dense k-grids are the correct QE band-gap baseline; high-symmetry paths are
  for band shapes.

## Layout

```
qe-convergence-attribution/
├── SKILL.md
├── README.md
├── requirements.txt
├── config.example.yaml
├── scripts/      # gap_controlled / force_matrix / plot_results
├── examples/     # 3 WBM semiconductor structures with k-paths
└── references/   # bandgap-controlled, high-force-matrix, workflow
```
