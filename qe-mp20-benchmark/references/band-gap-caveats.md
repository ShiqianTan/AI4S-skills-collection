# Band-gap caveats

## Use a dense uniform grid, not a high-symmetry path

A high-symmetry k-path can miss the true conduction-band minimum (off-path CBM)
and overestimate the gap. On 20 WBM semiconductors fixed at the VASP final
structure:

| Method | gap MAE | RMSE | max |
|---|---:|---:|---:|
| high-symmetry path | 0.1317 eV | 0.1798 eV | 0.4088 eV |
| dense 0.10 Å⁻¹ | **0.0665 eV** | 0.0935 eV | 0.2617 eV |
| dense 0.08 Å⁻¹ | **0.0764 eV** | 0.1126 eV | 0.2803 eV |

The path method roughly doubles the error. Use at least 0.08-0.10 Å⁻¹.

## Pseudopotentials are not the main error source

The same 20 systems with SSSP efficiency, SSSP precision and PseudoDojo NC-SR
give gap MAEs of 0.1317, 0.1296 and 0.1305 eV (spread <0.002 eV). Outliers are
identical across libraries. Do not expect a pseudopotential swap to fix a
systematic gap offset.

## The MP `band_gap` value is not a controlled PBE reference

MP gap values may mix functionals, structures and workflows. Examples:

- `CdTe mp-406`: MP gap 1.2285 eV contradicts other CdTe polymorphs (~0.62 eV).
- `Al2Ga2Se4`: the WBM *initial* structure gives 0.449-0.456 eV (VASP 0.45),
  while the VASP final gives 0.18 eV. The reference likely corresponds to a
  different structure.

For a rigorous test, compare on the **same structure**, with the **same
functional** (PBE), the **same k treatment** and a **matched pseudopotential
family**. This skill's dense nscf is the right QE side; the VASP side must be
controlled separately.

## Post-relaxation gap

For MPRelaxSet-like runs, the vc-relax uses smearing. Run a fixed-occupation
SCF, then the nscf, so the gap is not computed from a smeared density.

## GaN

GaN needs `nosym=.true.` and `diagonalization='cg'` for the nscf; otherwise the
MP k-grid can become a large P1 mesh.
