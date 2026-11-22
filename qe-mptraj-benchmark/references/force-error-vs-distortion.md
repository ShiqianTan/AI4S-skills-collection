# Force error grows with distortion

MPTraj-20, QE vs VASP static forces, grouped by the VASP maximum force:

| VASP max |F| (eV/Å) | frames | force MAE | force RMSE | max error |
|---|---:|---:|---:|---:|
| <0.2 | 31 | **0.031** | 0.069 | 0.428 |
| 0.2-0.5 | 12 | 0.075 | 0.160 | 0.747 |
| 0.5-1.0 | 9 | 0.215 | 0.493 | 1.975 |
| >1.5 | 7 | **0.682** | 1.384 | 6.372 |

## Why

High-force frames are highly sensitive to:

- pseudopotential choice,
- k-points and smearing,
- the electronic state reached by the SCF,
- the optimization path / local structure.

A one-factor-at-a-time scan of k (0.15-0.30 Å⁻¹), ecut (1.0-1.5×), smearing
(Gaussian/MV), and pseudopotential (SSSP efficiency/precision, PseudoDojo) changes
the force MAE by only 0.003-0.021 eV/Å on the five worst frames — about
0.4-2.6 % — while the maximum error stays at 6.5-6.8 eV/Å. The disagreement is
therefore not a numerical-convergence problem.

## Practical rules

- Validate forces on near-equilibrium structures; report the |F| distribution.
- For distorted frames, use same-structure comparison and flag `high-force`.
- Do not use tetrahedra for forces: QE warns, some frames print no forces, and
  SCF can fail.
- For problematic SCF convergence (e.g. TiMoPd), use smearing + `local-TF` +
  CG.

## Related

The high-force factor scan is packaged in `qe-convergence-attribution`.
