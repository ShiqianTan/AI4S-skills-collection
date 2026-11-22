# Basin-recovery workflow

## Why basin checks are required

A lower total energy does **not** mean the same structure. QE and VASP can reach
different local minima from the same distorted initial structure. Always judge
the basin with `StructureMatcher` (this skill uses
`primitive_cell=False, scale=False, attempt_supercell=False, ltol=0.2,
stol=0.3, angle_tol=5`).

## Recommended workflow

1. Start from the VASP final if available.
2. Otherwise start from the WBM initial **and** a perturbed final
   (3 % volume + 0.08 Å random displacement).
3. For anomalous systems add a second perturbation and a volume scan
   (0.95, 0.97, 1.03, 1.05 ×).
4. Classify each system by structure match, not energy.
5. Flag mismatches and verify with a third method or a VASP single point.

## Measured recovery rates

### WBM-20 (original)

| Start | n | match VASP final | mean volume difference |
|---|---:|---:|---:|
| VASP final | 20 | **95 %** (19/20) | +0.38 % |
| WBM initial | 20 | **80 %** (16/20) | -0.67 % |
| perturbed | 20 | **85 %** (17/20) | +1.35 % |
| volume scan | 16 | **100 %** (16/16) | +0.43 to +0.93 % |

- Any start recovers 95 % (19/20); only Al2S4 fails all three conventional
  starts.
- The lowest-QE-energy start matches the VASP basin only 85 % of the time — do
  not use energy to pick the basin.
- PtRbS2: initial start lands in another basin (133 meV/at higher, 14.2 %
  smaller); VASP-final start and volume scan recover it.
- I6Mg2: all starts match; volume scan 1.03× gives the lowest energy.
- Al2Pb2Se4/Rb4Se4: perturbations find a lower non-VASP basin; volume scan
  recovers the VASP basin at higher energy.

### 40 conventional structures

| Set | n | initial-start match | mean volume diff | mean RMS |
|---|---:|---:|---:|---:|
| Original WBM-20 | 20 | 80 % | -0.67 % | 0.0144 |
| Additional 20 | 20 | **100 %** | +0.20 % | 0.0031 |
| Combined 40 | 40 | **90 %** | -0.24 % | 0.0081 |

### 20 high-distortion structures

| Setting | match VASP final | mean volume diff |
|---|---:|---:|
| fast vc-relax | 26.3 % (5/19 valid) | +23.7 % |
| MPRelaxSet-like | 30.0 % (6/20) | +20.3 % |
| MPRelaxSet + VASP-final multi-start | **85 %** (17/20) | — |

- High-distortion failures come with huge volume changes (+50 to +114 %).
- k-points, cutoffs and convergence thresholds do not fix them.
- VASP-final starts rescue 11/14 of the MPRelaxSet failures; the remaining
  KF3/SrF3/InS fail every start — likely a real QE/VASP PES difference.

## Production rules

- With a VASP final: start there.
- Without one, and for distorted structures: initial + perturbed multi-start and
  a volume scan; do not trust a single initial-start run.
- Use QE tetrahedra **never** for relaxations/forces; it cannot provide reliable
  forces.
- Report "all systems" and "matched systems" metrics separately and list
  outliers.
