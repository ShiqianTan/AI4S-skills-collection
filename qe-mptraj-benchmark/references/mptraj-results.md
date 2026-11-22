# MPtrj results

QE 7.6 + SSSP 1.3.0 PBE efficiency static SCF at exact MPtrj VASP frame
geometries. Forces and stress compared directly; energies only within a
material.

## MPTraj-5

Materials: AgI, MgTe, InBr, CdTe (wurtzite), AlAs; 3 frames each (lowest energy,
highest energy, largest force).

| Metric | Value |
|---|---:|
| Force MAE | **0.0075 eV/Å** |
| Force RMSE | 0.0283 eV/Å |
| Max force-component error | 0.184 eV/Å (InBr) |
| Stress MAE (after sign/unit conversion) | **2.57e-4 eV/Å³** |
| Relative-energy MAE | **0.39-2.04 meV/atom** per material |

Most frames are near equilibrium (several have zero VASP forces); the statistics
are dominated by InBr, AgI and MgTe.

## MPTraj-20 (60 frames)

Materials include AgI, MgTe, InBr, CdTe, AlAs, TiMoPd, CuWSe4, KVSe4, NaIn2Bi,
SiPt6 and others.

| Metric | SSSP | MP-like |
|---|---:|---:|
| Force MAE | 0.1428 eV/Å | 0.1454 eV/Å |
| Force RMSE | 0.5171 eV/Å | 0.5290 eV/Å |
| Max force error | 6.37 eV/Å | 6.59 eV/Å |
| Stress MAE | 2.467e-3 eV/Å³ | 2.481e-3 eV/Å³ |
| Relative energy | 0.1-90 meV/atom | 0.1-95 meV/atom |

The two configurations are nearly identical for forces/stress. Relative
energies are more sensitive to the coarse MP k-grid (e.g. YGaPt2 16 → 72
meV/atom).

## Error vs distortion

See `force-error-vs-distortion.md`. The key result: near-equilibrium frames are
excellent (MAE 0.031 eV/Å for |F|<0.2 eV/Å), while frames above 1.5 eV/Å have
MAE 0.68 eV/Å and a 6.4 eV/Å maximum. This is not repaired by k/ecut/smearing/
pseudopotential changes (`qe-convergence-attribution`).

## Interpretation

- QE static SCF is reliable for near-equilibrium force/stress validation.
- Same-material relative energies are good (sub-meV to a few meV/atom for
  well-behaved sets), but coarse MP k-grids can inflate them.
- Do not use a handful of high-force frames to conclude that QE is broken; the
  disagreement is dominated by electron-state / local-minimum differences.
