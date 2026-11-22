# Formation energy and MP2020

## Definition

Raw QE formation energy per atom:

```
Ef_raw = (E_compound - Σ n_el · E_el_QE) / N_atoms
```

where `E_el_QE` are QE-relaxed **elemental ground states** (references in
`data/mp20_20/ref_*.cif`), not VASP/MP energies. Absolute QE and VASP total
energies differ by pseudopotential zeros, so only formation energies and
same-material relative energies are comparable.

## MP2020 correction

The Materials Project applies the MP2020 anion/oxide correction (and fitted
elemental references). Apply the same correction to QE formation energies:

```python
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
# scripts/apply_mp2020.py wraps this for the collected CSV
```

The correction removes almost all of the raw disagreement for the systems where
it applies:

| Stage | MAE (SSSP) |
|---|---:|
| Raw QE vs MP | 0.1724 eV/atom |
| QE + MP2020 vs MP | **0.0145 eV/atom** |

MP-like settings give 0.0131 eV/atom.

## Rules

- Always report raw and corrected values separately.
- Use consistent elemental references; O/N gas-phase references are a common
  source of apparent disagreement.
- Do not compare absolute total energies across codes.
- `e_above_hull` requires the full competing phase diagram and consistent
  references; it is not reproduced by this skill.
- WBM formation-energy MAE was 27.9 meV/atom (max 112 meV/atom), dominated by
  the basin-mismatched PtRbS2 system.
