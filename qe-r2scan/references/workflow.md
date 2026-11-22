# R2SCAN workflow details

## 1. Select

`matpes_r2scan_select.py` reads the MatPES-PBE `jobs.json` to get the same 40
`matpes_id` values, streams `MatPES-R2SCAN-2025.1.json`, and writes:

```
runs/matpes_r2scan/<group>/<sid>/
  structure.cif
  r2scan_frame.json   # R2SCAN energy/forces/stress
  scf.in              # input_dft='r2scan', PseudoDojo, 60/240 Ry
runs/matpes_r2scan/jobs.json
```

If the PBE jobs file is unavailable, use a MatPES selection from
`qe-matpes-benchmark` first.

## 2. Run

`run_r2scan.py` runs each SCF with `qe.libxc_bin_dir/pw.x` via `qe.mpirun`.
Production: `--np 4 --concurrent 16`, `OMP_NUM_THREADS=4`.

The script parses ASE `espresso-out`, compares forces, stress (kbar → eV/Å³,
QE sign flip) and within-material relative energies, and writes
`matpes_r2scan_frames.{json,csv}` and `matpes_r2scan_summary.json`.

## 3. Benchmark

`r2scan_bench.py --base-input .../scf.in` rewrites ecut/k and runs:

| tag | ecut | k | MPI | OpenMP |
|---|---:|---:|---:|---:|
| e60_k9_4x4 | 60 | 9 | 4 | 4 |
| e60_k9_8x4 | 60 | 9 | 8 | 4 |
| e60_k9_4x8 | 60 | 9 | 4 | 8 |
| e60_k9_2x8 | 60 | 9 | 2 | 8 |
| e60_k9_1x16 | 60 | 9 | 1 | 16 |
| e45_k7_4x4 | 45 | 7 | 4 | 4 |
| e45_k7_8x4 | 45 | 7 | 8 | 4 |

It reports wall time, energy, Fmax, and differences vs the 60 Ry/9×9×9 baseline.

## Failure modes

- `libxc needed for this functional` → wrong build (GPU or no libxc).
- Pseudopotential error → non-NC pseudo; use PseudoDojo NC.
- Missing frames → the R2SCAN dataset does not contain all PBE `matpes_id`s;
  the script warns and continues with the intersection.
