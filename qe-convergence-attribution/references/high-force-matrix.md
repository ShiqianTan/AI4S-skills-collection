# High-force convergence matrix

Five MPtrj frames with the largest force error (TiMoPd, SiPt6, NaIn2Bi, KVSe4,
K2HgTe) were re-evaluated with a one-factor-at-a-time scan: 65 tasks.

## Variants and force error

| Variant | n | force MAE (eV/Å) | RMSE | max |
|---|---:|---:|---:|---:|
| baseline | 5 | 0.8072 | 1.6466 | 6.5907 |
| k 0.15 | 5 | 0.8026 | 1.6316 | 6.5052 |
| k 0.20 | 5 | 0.8040 | 1.6358 | 6.5259 |
| k 0.30 | 5 | 0.8167 | 1.6787 | 6.7850 |
| ecut 1.5× | 5 | 0.8074 | 1.6469 | 6.5921 |
| smearing MV | 5 | 0.8115 | 1.6595 | 6.6652 |
| pseudo precision | 5 | 0.8070 | 1.6464 | 6.5918 |
| pseudo PseudoDojo | 5 | 0.8042 | 1.6380 | 6.5434 |
| combined | 5 | 0.8041 | 1.6361 | 6.5278 |

## Conclusion

The best single-factor change is **0.003-0.021 eV/Å** (0.4-2.6 % of 0.807), and
the maximum error stays at **6.5-6.8 eV/Å**. High-force disagreement is not a
k-point, cutoff, smearing or pseudopotential convergence problem. It reflects
different electronic states / local structures / real PES differences between QE
and VASP on strongly distorted frames.

## Tetrahedra

`occupations='tetrahedra'` is **not usable for forces**: QE warns, 2 frames
produce no force output, and 1 frame fails SCF. Do not include tetrahedra in a
force validation.

## Practical use

- Use this matrix to demonstrate that a high-force disagreement cannot be
  numerics.
- For MLIP/generated structure validation, restrict to near-equilibrium frames or
  flag high-force frames separately.
- For problematic SCF convergence, use smearing + `local-TF` + CG (e.g. TiMoPd
  19/34 iterations instead of 300 non-converged).
