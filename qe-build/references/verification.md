# Build verification

## Si diamond smoke test

Input: `examples/si/si.in` (PBE) or `examples/si/si_r2scan.in` (R2SCAN).
Pseudo: `pseudo/Si_r.upf`.

| Quantity | CPU | GPU |
|---|---|---|
| SCF iterations | 6 | 6 |
| Total energy (Ry) | -16.92302587 | -16.92302587 |
| Wall time | 0.82 s | 2.25 s (tiny system, launch overhead) |
| Banner | — | `GPU acceleration is ACTIVE` |
| Status | `JOB DONE` | `JOB DONE` |

Energy agreement is to the printed precision, validating that the GPU build is
not changing the physics.

## R2SCAN check

Si with `input_dft='r2scan'` must print `Exchange-correlation = R2SCAN` and
reach `JOB DONE`. Do not compare its absolute energy with the PBE value; the
functional differs.

## Checklist before benchmarking

- [ ] CPU Si energy equals `-16.92302587 Ry`.
- [ ] GPU Si energy equals the CPU value and the GPU banner is present.
- [ ] libxc Si R2SCAN converges with the R2SCAN banner.
- [ ] `pw.x`, `ph.x`, `pp.x` exist under each prefix's `bin/`.
- [ ] Shared `config.yaml` points at the correct build for each task.
