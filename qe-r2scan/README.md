# qe-r2scan

Run the MatPES-R2SCAN benchmark with QE libxc + PseudoDojo NC pseudopotentials.

## Install

```bash
pip install -r requirements.txt
python ../pseudos/fetch_pseudos.py
```

Needs `ijson` for streaming the MatPES JSON and a CPU libxc QE build
(`qe-build/scripts/build_libxc.sh`).

## Workflow

```bash
cp config.example.yaml config.yaml   # set qe.libxc_bin_dir, datasets.matpes, pseudos
python scripts/matpes_r2scan_select.py --config config.yaml --pbe-jobs ../qe-matpes-benchmark/examples/jobs.json
python scripts/run_r2scan.py --config config.yaml --r2scan-dir runs/matpes_r2scan --np 4 --concurrent 16
python scripts/plot_results.py --config config.yaml --tag matpes_r2scan
```

## Scripts

| Script | Purpose |
|---|---|
| `matpes_r2scan_select.py` | extract same MatPES IDs from R2SCAN, write `scf.in`/`jobs.json` |
| `run_r2scan.py` | run SCF on the libxc build, compare forces/stress/energies |
| `r2scan_bench.py` | benchmark MPI×OpenMP/cutoff/k-grid configurations |
| `plot_results.py` | force/stress parity and PBE-vs-R2SCAN bars |

## Bundled example

`examples/` has three materials (12 frames) with `structure.cif`,
`r2scan_frame.json`, `scf.in`, and a ready `jobs.json`:

```bash
python scripts/run_r2scan.py --config config.yaml --r2scan-dir examples --np 1 --concurrent 1
```

## Notes

- R2SCAN requires `input_dft='r2scan'` and norm-conserving pseudopotentials.
- Use the CPU libxc build; GPU NVHPC cannot read gfortran libxc modules.
- R2SCAN force MAE is about 4× PBE on MatPES; it is a qualitative cross-check.
- See `references/r2scan-results.md` and `references/r2scan-limitations.md`.
