# _shared — `qe_bench_core`

Shared Python utilities used by every `qe-*` benchmark skill in this collection.
They are extracted from the original `quantum_espresso/repro/scripts` workflow so
that skills stay small and consistent.

## Modules

| Module | Purpose |
|---|---|
| `paths` | config/env resolution (`config.example.yaml`), QE/MPI/dataset paths, subprocess env |
| `structure` | CIF/POSCAR reading, pymatgen↔ASE conversion |
| `qeio` | cutoffs, k-grids, `pw.x` input generation (scf/relax/nscf/bands) |
| `runner` | `pw.x` execution, JMP DONE detection, scratch handling |
| `parse` | energy/force/stress/gap parsing, band blocks |
| `metrics` | MAE/RMSE/bias/relative-energy helpers |
| `mp2020` | Materials Project 2020 compatibility corrections |
| `plotting` | parity plots |
| `pseudos` | SSSP manifests, PseudoDojo maps, valence parsing, tarball extraction |

## Using the core from a skill

Each skill's `scripts/_bootstrap.py` inserts the repository root on
`sys.path` and loads the user's `config.yaml`:

```python
from _bootstrap import cfg, ROOT, RESULTS, WORK
from qe_bench_core import qeio, runner, parse, metrics
```

Path resolution priority: `--config` > `QE_BENCH_CONFIG` > `QE_*` env vars >
`config.yaml` > defaults. See `config.example.yaml`.
