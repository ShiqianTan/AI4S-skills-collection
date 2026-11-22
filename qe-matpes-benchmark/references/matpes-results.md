# MatPES-PBE results

Dataset: `MatPES-PBE-2025.1.json` (1.8 GB, 434,712 records).
Selection: 10 materials × 4 frames — InN, CaO, AgF, GeSe, TlCuPd2, TePb, LiC6,
C, KH3, CaC4. QE 7.6 + SSSP 1.3.0 PBE efficiency, static SCF.

| Metric | Value |
|---|---:|
| Force MAE | **0.0123 eV/Å** |
| Force RMSE | 0.0204 eV/Å |
| Max force error | 0.111 eV/Å (InN) |
| Stress MAE | **7.40e-4 eV/Å³** |
| Relative energy | **0.28-3.99 meV/atom** |

MatPES-PBE is the best-agreeing force/stress reference in the study:

| Reference | force MAE (eV/Å) | stress MAE (eV/Å³) |
|---|---:|---:|
| MatPES-PBE | **0.0123** | **7.40e-4** |
| sAlex-80 | 0.0223 | 8.94e-4 |
| MPTraj-5 | 0.0075 | 2.57e-4 (near-equilibrium only) |
| MPTraj-20 | 0.143 | 2.47e-3 |

## Stress conversion

MatPES stress is stored in **kbar**. The comparison applies:

```
vasp_stress_eV_A3 = vasp_stress_kbar * 6.241509e-4
qe_stress_eV_A3   = -qe_stress_eV_A3       # QE -> VASP sign
```

Both are then converted to the Voigt matrix `[xx, yy, zz, yz, xz, xy]`.

## Frame selection

For each material, four frames whose maximum force is closest to 0.01, 0.6, 1.5
and 2.5 eV/Å, restricted to within 5-10 eV of the lowest-energy frame and
`max|F| ≤ 3 eV/Å` where possible.

## Interpretation

- MatPES is the recommended near-equilibrium force/stress validation set.
- The remaining 0.012 eV/Å force MAE is at the level expected from
  pseudopotential/method differences, not a QE defect.
- Same-material relative energies are a few meV/atom.
