# sAlex results

Dataset: `/path/to/sAlex/val` (348 LMDB files, 553,218 frames in the study).
Selection: ≤8 atoms, SSSP-covered, non-magnetic → **20 materials × 4 frames =
80 frames**.

QE 7.6 + SSSP 1.3.0 PBE efficiency, static SCF. sAlex stress uses the ASE
convention, so no kbar conversion or sign flip.

| Metric | 40 frames | 80 frames |
|---|---:|---:|
| Force MAE | 0.0234 eV/Å | **0.0223 eV/Å** |
| Force RMSE | 0.0460 eV/Å | 0.0462 eV/Å |
| Max force error | 0.309 eV/Å | 0.645 eV/Å (K2PS5) |
| Stress MAE | 8.93e-4 eV/Å³ | **8.94e-4 eV/Å³** |
| Relative energy | 1.1-14.9 meV/atom | **0.79-19.71 meV/atom** |

## Interpretation

- sAlex is one of the best force/stress references in this study: force MAE
  0.022 eV/Å and stress MAE 8.9e-4 eV/Å³.
- Its distorted frames are much more consistent than the MPTraj-20 high-force
  group (MAE 0.68 eV/Å there).
- Relative energies are sub-20 meV/atom.
- No kbar conversion: sAlex stores stress in eV/Å³ with the ASE sign already.

## Selection recipe

1. Scan every `*.aselmdb` in `val/`, group by `sid` prefix.
2. Keep groups with ≥4 frames and atoms ≤8.
3. Drop elements outside SSSP and the excluded heavy/radioactive/magnetic set.
4. Rank groups by force spread and keep 20; take 4 frames each (low-energy,
   high-energy, high-force, and a fourth representative frame).

The exact selected IDs are embedded in `scripts/salex_select.py` (`SELECTED`).
