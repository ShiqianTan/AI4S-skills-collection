# qe-salex-benchmark

Validate QE forces, stress and same-material relative energies against sAlex.

## Install

```bash
pip install -r requirements.txt
python ../pseudos/fetch_pseudos.py
```

Needs `ase` (LMDB reader), `pymatgen`, `numpy`, `matplotlib`.

## Workflow

```bash
cp config.example.yaml config.yaml   # set datasets.salex and QE binaries
python scripts/salex_scan.py --config config.yaml
python scripts/salex_select.py --config config.yaml
python scripts/run_mptraj.py --config config.yaml --mptraj-dir runs/salex --tag salex --gpu-mode gpu
python scripts/plot_results.py --config config.yaml --tag salex
```

## Scripts

| Script | Purpose |
|---|---|
| `salex_scan.py` | scan `*.aselmdb`, group small non-magnetic materials |
| `salex_select.py` | select 20×4 frames, write `scf.in`/`jobs.json`/`salex_frame.json` |
| `run_mptraj.py` | shared static-SCF runner/parser (forces/stress/energy) |
| `plot_results.py` | force, force-magnitude, relative-energy and stress parity |

## Bundled example

`examples/` has three selected materials (12 frames) with `structure.cif`,
`salex_frame.json`, `scf.in`, and a ready `jobs.json`:

```bash
python scripts/run_mptraj.py --config config.yaml --mptraj-dir examples --tag salex_examples --gpu-mode cpu
python scripts/plot_results.py --config config.yaml --tag salex_examples
```

## Conventions

- sAlex stress is ASE convention (eV/Å³, no sign flip); jobs set
  `stress_unit: eV_A3` and `stress_sign: 1`.
- Only same-material relative energies are meaningful.
- See `references/salex-results.md` for measured accuracy.
