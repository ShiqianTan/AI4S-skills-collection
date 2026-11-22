---
name: qe-wbm-basin
description: Reproduce WBM IS2RE structure optimization with Quantum ESPRESSO and quantify basin recovery. Prepares WBM initial + VASP final structures, runs vc-relax from the initial structure and a static SCF at the VASP final, matches structures with StructureMatcher, computes MP2020 formation energies, and runs a basin-recovery benchmark (VASP-final / initial / perturbed / volume-scan starts) including high-distortion rescue. Use to test whether QE lands in the same local minimum as VASP and to decide on basin-check / multi-start workflows.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# WBM structure and basin recovery (qe-wbm-basin)

WBM IS2RE is a structure-relaxation benchmark: the input is an unrelaxed
structure, the reference is the VASP relaxed structure. This skill measures both
the geometry/energy agreement and, critically, whether QE recovers the **same
basin** as VASP.

```
WBM initial + VASP final -> QE vc-relax (initial) + SCF (VASP final)
                          -> structure match / volume / formation energy
Basin benchmark           -> starts: VASP-final, initial, perturbed, volume scan
                          -> recovery rate + high-distortion rescue
```

Reference results: `references/wbm-results.md` and `references/basin-workflow.md`.

## Quick start

```bash
python ../pseudos/fetch_pseudos.py

# 1) select structures + prepare inputs (base20 / medium20 / high20)
python scripts/wbm_select.py --config config.yaml --preset base20

# 2) initial relax + VASP-final SCF
python scripts/run_wbm.py --config config.yaml --wbm-dir runs/wbm --tag wbm --gpu-mode gpu

# 3) structure/energy comparison and formation energies
python scripts/analyze_wbm.py --config config.yaml --frames runs/results/wbm_frames.json

# 4) basin recovery (initial / VASP-final / perturbed, optional volume scan)
python scripts/basin_recovery.py --config config.yaml --wbm-dir runs/wbm \
  --tag wbm_basin --cases vaspfinal,initial,perturb --gpu-mode gpu

# 5) plots
python scripts/plot_results.py --config config.yaml --tag wbm
```

Bundled `examples/` has four already-prepared systems (PtRbS2, I6Mg2, Rb4Se4,
Al2S4) with initial/VASP-final CIFs and inputs, so `run_wbm.py` and
`basin_recovery.py` can be exercised without the WBM dataset:

```bash
python scripts/run_wbm.py --config config.yaml --wbm-dir examples --tag wbm_examples --gpu-mode cpu
python scripts/basin_recovery.py --config config.yaml --wbm-dir examples \
  --out runs/wbm_basin_examples --tag wbm_basin_examples --backend cpu
```

## Results (from the study)

WBM-20 (initial start), matched 16/20:

| Metric | All 20 | Matched 16 |
|---|---:|---:|
| Volume MAE | 2.17 % | **1.57 %** |
| Volume median | 0.54 % | 0.44 % |
| Lattice abc MAE | 0.195 Å | **0.085 Å** |
| Angle MAE | 3.68° | 0.33° |
| Structure match | 80 % | 100 % |
| Formation-energy MAE | 27.9 meV/atom | **16.7 meV/atom** |

Four mismatches: PtRbS2, I6Mg2 (QE wrong basin), Rb4Se4, Al2Pb2Se4 (QE finds a
lower minimum). Basin recovery:

| Start | n | match VASP final |
|---|---:|---:|
| VASP final | 20 | **95 %** (19/20) |
| WBM initial | 20 | **80 %** (16/20) |
| perturbed | 20 | **85 %** (17/20) |
| volume scan | 16 | **100 %** |

For 40 conventional structures, initial-start recovery is 90 % (36/40). For 20
**high-distortion** structures it is only 30 %; adding VASP-final starts raises
it to 85 % (17/20). See `references/basin-workflow.md`.

## Key rules

- Judge the basin by **structure matching**, not by energy; a lower QE energy
  may be a different (real) minimum.
- Prefer the VASP final as a start when available.
- Without a VASP final, use initial + perturbed multi-start and a volume scan for
  anomalous systems.
- WBM CSE has no VASP forces/stress, so the "force" metric is the residual QE
  force at the VASP final.

## Layout

```
qe-wbm-basin/
├── SKILL.md
├── README.md
├── requirements.txt
├── config.example.yaml
├── scripts/          # wbm_select / run_wbm / analyze_wbm / basin_recovery / plot_results
├── examples/         # 4 prepared systems (initial + VASP final + inputs)
└── references/       # results and basin workflow
```
