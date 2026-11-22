# WBM results

Data: WBM CSE (VASP final structure/energy), WBM summary (MP2020 formation
energy, hull, band gap, volume), and `WBM_IS2RE.aselmdb` (initial structures).
QE 7.6 + SSSP, MPRelaxSet-like vc-relax from the same WBM initial structure;
static SCF at the VASP final.

## WBM-20: initial-start geometry and energy

| Metric | All 20 | Matched 16 | Mismatched 4 |
|---|---:|---:|---:|
| Volume MAE | 2.17 % | **1.57 %** | 4.56 % |
| Volume median | 0.54 % | **0.44 %** | 2.15 % |
| Volume max | 13.43 % | 5.94 % | 13.43 % |
| Lattice abc MAE | 0.195 Å | **0.085 Å** | 0.633 Å |
| Angle MAE | 3.68° | **0.33°** | 17.1° |
| Structure match rate | 80 % | 100 % | 0 % |
| Structure RMS median | 0.0030 | **0.0030** | — |
| Residual QE force at VASP final (mean) | 0.0328 eV/Å | 0.0367 | 0.0173 |
| Formation-energy MAE | 27.9 meV/atom | **16.7 meV/atom** | 72.9 meV/atom |
| Formation-energy max | 112.3 meV/atom | 44.0 | 112.3 |
| QE(VASP final) - QE(QE final) MAE | 21.3 meV/atom | **7.05** | 78.2 |

The four mismatches:

| System | Original QE | From VASP final | From initial + improved settings | Diagnosis |
|---|---|---|---|---|
| PtRbS2 | -13.4 %, wrong basin | matches, dV -0.6 %, RMS 0.002 | still wrong basin, +150.5 meV/atom | basin selection |
| I6Mg2 | -1.6 %, wrong basin | matches, dV -0.4 %, RMS 0.0007 | still wrong basin, +40.6 meV/atom | basin selection |
| Rb4Se4 | wrong basin but lower energy | matches, dV +1.0 %, RMS 0.0004 | lower than VASP by 95.2 meV/atom | QE lower minimum |
| Al2Pb2Se4 | wrong basin but lower energy | matches, dV +1.7 %, RMS 0.0074 | QE lower by 21.3 meV/atom | QE lower minimum |

## What the numbers mean

- For structures that stay in the VASP basin, QE geometry is excellent: volume
  MAE 1.57 %, abc MAE 0.085 Å, angle MAE 0.33°, RMS median 0.003.
- Residual QE forces at the VASP final are small (~0.03 eV/Å), so the VASP final
  is a QE stationary point.
- 20 % of the original set lands in a different basin; this is the main WBM
  risk, not a numerical failure.

## WBM CSE limitation

The WBM computed-structure entries do not contain VASP forces or stress. The
"force" column above is the residual QE force at the VASP final structure, not a
VASP-vs-QE force difference.
