# Controlled band-gap attribution

A controlled experiment fixes the **structure** and varies one factor at a time.
This separates QE numerical/pseudopotential errors from reference-structure and
functional effects.

## Setup

20 WBM semiconductors (0.05-4.0 eV, five per gap quartile), `n_sites <= 8`,
main-group/non-magnetic. The structure is the **VASP final** from WBM CSE unless
`--start initial` is used. Reference: WBM `bandgap_pbe` (eV).

## Pseudopotentials are not the main error

| Pseudopotential | n | MAE | RMSE | bias | max |
|---|---:|---:|---:|---:|---:|
| SSSP efficiency | 20 | 0.1317 | 0.1798 | -0.0202 | 0.4088 |
| SSSP precision | 20 | 0.1296 | 0.1794 | -0.0181 | 0.4089 |
| PseudoDojo NC-SR | 20 | 0.1305 | 0.1777 | -0.0181 | 0.4069 |
| **all** | **60** | **0.1306** | **0.1790** | **-0.0188** | **0.4089** |

MAE spread <0.002 eV. The outliers (Al2Ge2Te2, Al2In2Te4, Al2Ga2Se4,
Al2CaTe4, As2Pb2Sr) are nearly identical across libraries. Do not expect a
pseudopotential swap to fix a systematic gap offset.

## Path vs dense grid

| Method | gap MAE | RMSE | max |
|---|---:|---:|---:|
| high-symmetry path | 0.1317 | 0.1798 | 0.4088 |
| dense 0.10 Å⁻¹ | **0.0665** | 0.0935 | 0.2617 |
| dense 0.08 Å⁻¹ | **0.0764** | 0.1126 | 0.2803 |

The path misses off-path conduction-band minima and roughly doubles the error.
Al2Ge2Te2: path 1.039 → dense 0.08 0.612 (VASP 0.63). Al2In2Te4: 0.580 → 0.075
(VASP 0.21). As2Se2Te: 0.621 → 0.393 (VASP 0.40).

## Initial vs final structure

Dense 0.08, all 20 systems:

| Start | gap MAE |
|---|---:|
| VASP final | **0.076 eV** |
| WBM initial | **0.614 eV** |

18/20 finals are closer; only Al2Ga2Se4 and Al2P1S4 are closer from the initial
structure.

## Residual attribution (3 outliers)

Final vs initial × two pseudopotentials, dense 0.08:

| System | VASP | final SSSP-eff | final PDOJO | initial SSSP-eff | initial PDOJO |
|---|---:|---:|---:|---:|---:|
| Al2Ga2Se4 | 0.45 | 0.177 | 0.186 | **0.449** | **0.456** |
| Al2CaTe4 | 1.59 | 1.310 | 1.331 | 0.860 | 0.874 |
| As2Pb2Sr | 0.28 | 0.103 | 0.102 | 0.027 | 0.025 |

- Pseudopotential effect ≤0.022 eV.
- **Al2Ga2Se4 is a reference-structure issue**: the initial structure gives
  0.449-0.456 eV (VASP 0.45), the final only 0.18 eV.
- Al2CaTe4 and As2Pb2Sr are not fixed by structure or pseudopotential
  (0.18-0.28 eV residual) — likely reference/functional or a real QE/VASP
  electronic-structure difference.

## Functional dependence (dense 0.10, PseudoDojo)

| System | PBE | PBEsol | R2SCAN | VASP |
|---|---:|---:|---:|---:|
| Al2Ga2Se4 | 0.197 | ~0.20 | 0.459 | 0.45 |
| Al2CaTe4 | 1.422 | ~1.42 | 1.687 | 1.59 |
| As2Pb2Sr | 0.156 | ~0.15 | 0.109 | 0.28 |

PBEsol tracks PBE; R2SCAN opens the gap for the first two but is worse for
As2Pb2Sr. No single functional fixes all residuals.

## Recommended report

Always state the structure, functional, k treatment and pseudopotential family of
both sides. The QE side should use dense k (0.08-0.10 Å⁻¹). MP `band_gap` CSV
values are not a controlled PBE reference.
