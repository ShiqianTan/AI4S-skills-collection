# qe-matpes-benchmark

Validate QE forces, stress and same-material relative energies against
MatPES-PBE.

## Install

```bash
pip install -r requirements.txt
python ../pseudos/fetch_pseudos.py
```

Needs `ijson` for streaming the MatPES JSON.

## Workflow

```bash
cp config.example.yaml config.yaml   # set datasets.matpes and QE binaries
python scripts/matpes_scan.py --config config.yaml
python scripts/matpes_select.py --config config.yaml
python scripts/run_mptraj.py --config config.yaml --mptraj-dir runs/matpes --tag matpes --gpu-mode gpu
python scripts/plot_results.py --config config.yaml --tag matpes
```

## Scripts

| Script | Purpose |
|---|---|
| `matpes_scan.py` | stream MatPES, summarize small non-magnetic groups |
| `matpes_select.py` | select 10×4 frames, write `scf.in`/`jobs.json`/`matpes_frame.json` |
| `run_mptraj.py` | shared static-SCF runner/parser |
| `plot_results.py` | force, force-magnitude, relative-energy and stress parity |

## Bundled example

`examples/` has three selected materials (12 frames) with `structure.cif`,
`matpes_frame.json`, `scf.in`, and a ready `jobs.json`:

```bash
python scripts/run_mptraj.py --config config.yaml --mptraj-dir examples --tag matpes_examples --gpu-mode cpu
python scripts/plot_results.py --config config.yaml --tag matpes_examples
```

## Stress convention

MatPES stress is in kbar. Jobs set `stress_unit: kbar`, `stress_sign: -1`.
`plot_results.py` applies `KBAR_TO_EV_A3 = 6.241509e-4` and negates the QE
stress. See `references/matpes-results.md`.
