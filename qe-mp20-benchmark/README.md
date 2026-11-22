# qe-mp20-benchmark

Reproduce Materials Project MP20 geometry, band gap and formation energy with
Quantum ESPRESSO.

## Install

```bash
pip install -r requirements.txt
python ../pseudos/fetch_pseudos.py
```

Requires a QE build (`qe-build`) and the shared `config.yaml`.

## Workflow

```bash
cp config.example.yaml config.yaml   # edit QE/paths
python scripts/mp20_prepare.py --config config.yaml --mode sssp
python scripts/run_mp20.py --config config.yaml --mp20-dir runs/mp20_20 --gpu-mode gpu --devices 0,1,2,3,4,5,6,7
python scripts/collect_mp20.py --config config.yaml --mp20-dir runs/mp20_20 --tag mp20_20
python scripts/apply_mp2020.py --config config.yaml --tag mp20_20
python scripts/plot_results.py --config config.yaml --tag mp20_20
```

Use `--mode mprelax` for the MPRelaxSet-like configuration, and
`--gpu-mode cpu` when no GPU build is available.

## Scripts

| Script | Purpose |
|---|---|
| `mp20_prepare.py` | write `vc_relax.in` + `manifest.json` for 5 or 20 systems and elemental references |
| `run_mp20.py` | run vc-relax, optional fixed SCF, dense nscf; write `run_status.json` |
| `collect_mp20.py` | compare with MP; raw formation energies; summary |
| `apply_mp2020.py` | apply MP2020 anion/oxide corrections |
| `plot_results.py` | gap/volume/formation-energy parity plots |
| `convergence_test.py` | k-point and cutoff gap convergence for AlAs/AgI |
| `precision_test.py` | re-run AlAs/AgI with SSSP precision pseudos |

## Notes

- The band gap must come from a dense uniform nscf, not a high-symmetry path
  (the latter misses off-path CBM and inflates MAE by about 2×).
- Raw QE formation energies must not be compared directly with MP; apply
  MP2020 and use consistent elemental references.
- `e_above_hull` is out of scope for this skill.

See `SKILL.md` and `references/` for the full method and measured results.
