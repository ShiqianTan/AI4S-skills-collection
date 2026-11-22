# MatPES-R2SCAN results

Same 40 `matpes_id` frames as MatPES-PBE. QE R2SCAN uses PseudoDojo NC-SR PBE
standard pseudopotentials and the CPU libxc build.

## PBE vs R2SCAN

| Metric | PBE | R2SCAN |
|---|---:|---:|
| Force MAE | **0.0123 eV/Å** | 0.0489 eV/Å |
| Force RMSE | 0.0204 eV/Å | 0.0916 eV/Å |
| Max force error | 0.111 eV/Å | 0.726 eV/Å |
| Stress MAE | **7.40e-4 eV/Å³** | 3.59e-2 eV/Å³ |
| Relative energy | 0.28-3.99 meV/atom | 1.5-48 meV/atom |

R2SCAN force MAE is about 4× PBE; stress MAE is about 50× larger. This is far
from quantitative agreement.

## Causes

1. **Pseudopotentials**: PseudoDojo NC-SR PBE pseudopotentials, not an R2SCAN-
   specific set.
2. **Different VASP setup**: MP R2SCAN PAW settings differ from the QE NC setup.
3. **Convergence**: 60 Ry / 7×7×7 is not fully converged for meta-GGA.
4. **Implementation**: libxc R2SCAN vs VASP PAW R2SCAN differ.

## Production setting

After a CPU benchmark, 60 Ry / 7×7×7 (0.25 Å⁻¹) was chosen: about 2× faster than
60 Ry / 9×9×9 with a force difference of about 0.01 eV/Å. Parallelization:
16 concurrent jobs × 4 MPI × 4 OpenMP = 256 threads.

## Functional sensitivity in the gap test

In the controlled WBM gap test (dense 0.10, PseudoDojo), R2SCAN opened the gap
for Al2Ga2Se4 (0.197 → 0.459 vs VASP 0.45) and Al2CaTe4 (1.422 → 1.687 vs VASP
1.59), but made As2Pb2Sr worse (0.156 → 0.109 vs VASP 0.28). No single
functional fixes all residuals; confirm the reference's functional first.

## Interpretation

R2SCAN is useful for a **qualitative** cross-check (trends, whether a gap opens),
but should not be used as the sole quantitative reference until R2SCAN-compatible
pseudopotentials and convergence are established.
