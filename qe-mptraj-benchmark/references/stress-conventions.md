# Stress and energy conventions

## Units

- VASP/MPtrj stress is reported in **kbar**.
- QE `tstress` reports **Ry/bohr³**; ASE's `get_stress()` returns eV/Å³.
- 1 kbar = 6.241509×10⁻⁴ eV/Å³ (`KBAR_TO_EV_A3 = 0.0006241509`).

## Sign

QE stress has the **opposite sign** to the VASP/MPtrj convention. In the
benchmark the QE Voigt matrix is negated before comparison:

```python
qe_mat = voigt_to_matrix(atoms.get_stress())
qe_mat = -qe_mat                 # QE -> VASP convention
vasp_mat = vasp_stress * KBAR_TO_EV_A3
```

sAlex uses the ASE convention already, so no kbar conversion or sign flip is
applied there.

## Voigt order

The helper uses `[xx, yy, zz, yz, xz, xy]`:

```python
[[v0, v5, v4],
 [v5, v1, v3],
 [v4, v3, v2]]
```

MatPES stores stress in kbar and needs the same conversion and sign flip as
MPtrj.

## Energies

Absolute VASP and QE total energies differ by pseudopotential offsets. Compare:

- forces (offset-free),
- stress (offset-free),
- **relative** energies within one material: subtract each method's minimum.

The summary reports a per-material `energy_offset_eV_atom` only for diagnostics;
it is not a quality metric.
