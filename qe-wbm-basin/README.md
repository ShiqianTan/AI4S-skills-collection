# qe-wbm-basin

Reproduce WBM IS2RE structure optimization with QE and quantify basin recovery.

## Install

```bash
pip install -r requirements.txt
python ../pseudos/fetch_pseudos.py
```

Needs the matbench-discovery WBM directory and `WBM_IS2RE.aselmdb`.

## Workflow

```bash
cp config.example.yaml config.yaml   # set datasets.wbm / wbm_is2re and QE binaries
python scripts/wbm_select.py --config config.yaml --preset base20
python scripts/run_wbm.py --config config.yaml --wbm-dir runs/wbm --tag wbm --gpu-mode gpu
python scripts/analyze_wbm.py --config config.yaml --frames runs/results/wbm_frames.json
python scripts/basin_recovery.py --config config.yaml --wbm-dir runs/wbm --tag wbm_basin --gpu-mode gpu
python scripts/plot_results.py --config config.yaml --tag wbm
```

`wbm_select.py` presets: `base20` (original 20), `medium20` (40-structure
addition), `high20` (high-distortion stress test).

## Scripts

| Script | Purpose |
|---|---|
| `wbm_select.py` | read WBM CSE/IS2RE/summary; write initial/VASP-final CIFs + inputs + manifest |
| `run_wbm.py` | vc-relax from initial + SCF at VASP final; structure matching; summary |
| `analyze_wbm.py` | MP2020 formation energies vs WBM |
| `basin_recovery.py` | multi-start basin benchmark (VASP-final/initial/perturb/volume scan) |
| `plot_results.py` | volume/formation-energy parity and basin recovery bar chart |

## Bundled example

`examples/` contains four prepared systems (PtRbS2, I6Mg2, Rb4Se4, Al2S4) with
`initial.cif`, `vasp_final.cif`, `vc_relax.in`, `scf_vasp_final.in`, and a
manifest. Run with `--wbm-dir examples` and `--gpu-mode cpu`.

## Basin rules

- Use `StructureMatcher`, not energy, to decide if the basin is the same.
- Prefer the VASP final as a start.
- For distorted systems without a VASP final, use initial + perturbed starts and
  a volume scan; some systems still fail (QE vs VASP PES difference).

See `references/wbm-results.md` and `references/basin-workflow.md`.
