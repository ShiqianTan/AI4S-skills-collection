# qe-mptraj-benchmark

Validate QE forces, stress and same-material relative energies against MPtrj
VASP trajectories.

## Install

```bash
pip install -r requirements.txt
python ../pseudos/fetch_pseudos.py
```

Requires `ijson` for streaming the full MPtrj JSON.

## Workflow

```bash
cp config.example.yaml config.yaml
python scripts/mptraj_select.py --config config.yaml --preset 20
python scripts/mptraj_prepare.py --config config.yaml --mptraj-dir runs/mptraj_20 --mode sssp
python scripts/run_mptraj.py --config config.yaml --mptraj-dir runs/mptraj_20 --tag mptraj20 --gpu-mode gpu
python scripts/plot_results.py --config config.yaml --tag mptraj20
```

Use `--preset 5` for the smaller chemistry subset, `--mode mprelax` for
MPRelaxSet-like settings, and `--gpu-mode cpu` without a GPU build.

## Bundled example

`examples/manifest.json` plus four materials × three frames (AlAs, AgI, MgTe,
InBr) with VASP forces/stress. This lets you run the pipeline without the 1.8 GB
MPtrj file:

```bash
python scripts/mptraj_prepare.py --config config.yaml --mptraj-dir examples --mode sssp
python scripts/run_mptraj.py --config config.yaml --mptraj-dir examples --tag mptraj_examples --gpu-mode cpu
python scripts/plot_results.py --config config.yaml --tag mptraj_examples
```

## Scripts

| Script | Purpose |
|---|---|
| `find_candidates.py` | scan MPtrj for small non-magnetic high-force materials |
| `mptraj_select.py` | select frames + write `manifest.json`/`structure.cif`/`vasp_frame.json` |
| `mptraj_prepare.py` | write `scf.in` and `jobs.json` |
| `run_mptraj.py` | run SCF, compare forces/stress/energy, write frames + summary |
| `plot_results.py` | force parity, error vs Fmax, relative-energy parity |

## Outputs

`{tag}_frames.json`, `{tag}_frames.csv`, `{tag}_summary.json`, and three PNGs in
`output.results_dir`.

## Caveats

- Same-material relative energies only; absolute energies are pseudopotential
  dependent.
- Stress needs the QE sign flip and kbar → eV/Å³ conversion.
- High-force frames are not a numerical-settings problem; see
  `references/force-error-vs-distortion.md`.
