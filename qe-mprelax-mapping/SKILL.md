---
name: qe-mprelax-mapping
description: Map the Materials Project MPRelaxSet (pymatgen commit c4998d92) to Quantum ESPRESSO parameters. Converts ENCUT, EDIFF_PER_ATOM, SIGMA, KPOINTS reciprocal density, NSW/IBRION/ISIF and NELM into ecutwfc/ecutrho, conv_thr, degauss, k-grid and vc-relax settings, and records the intentional deviations (ISMEAR=-5, ISPIN=2, LDAU, ALGO=Fast). Use before any MP20/MPtrj/MatPES/WBM run so the QE settings are traceable to the VASP reference.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# MPRelaxSet → Quantum ESPRESSO mapping (qe-mprelax-mapping)

The Materials Project reference values were produced with the `MPRelaxSet`
(`pymatgen` commit `c4998d92`). This skill turns that INCAR/KPOINTS definition
into a QE configuration that is as close as the two codes allow, and documents
where they necessarily diverge.

Source definition: [`data/MPRelaxSet_c4998d92.yaml`](data/MPRelaxSet_c4998d92.yaml).

## Mapping table

| MPRelaxSet (VASP) | QE setting | Rule |
|---|---|---|
| `ENCUT = 520 eV` | `ecutwfc`, `ecutrho` | `ecutwfc = max(520/13.6057, max SSSP cutoff_wfc)`; `ecutrho = max(8·ecutwfc, max SSSP cutoff_rho)` |
| `EDIFF_PER_ATOM = 5e-5 eV` | `conv_thr` | `max(1e-8, 5e-5·natoms/13.6057)` Ry |
| `SIGMA = 0.05 eV` | `degauss` | `0.05/13.6057` Ry, Gaussian smearing |
| `KPOINTS reciprocal_density = 64` | k-grid + shift | `MPRelaxSet(structure).kpoints`; gamma-centred mesh → shift (0,0,0), else (1,1,1) |
| `NELM = 100` | `electron_maxstep = 100` | |
| `NSW = 99`, `IBRION = 2`, `ISIF = 3` | `vc-relax`, BFGS, `nstep = 99` | |
| `ALGO = Fast` | Davidson (CG for difficult systems) | |
| `PREC = Accurate` | default cutoffs | |
| `EDIFFG`/force threshold | `forc_conv_thr = 5e-4 Ry/bohr` | |

### Intentional deviations

| VASP | QE choice | Reason |
|---|---|---|
| `ISMEAR = -5` (tetrahedra) | Gaussian smearing 0.05 eV | QE tetrahedra gives unreliable forces for relaxations; documented in `references/deviations.md` |
| `ISPIN = 2` | `nspin = 1` | Selected benchmark systems are non-magnetic; enable `nspin=2` for magnetic systems |
| `LDAU = true` | not enabled | Selected systems have no Co/Cr/Fe/Mn/Mo/Ni/V/W in an F/O environment |
| `LWAVE=false`, `LORBIT=11` | not mapped | post-processing only |

### Extra benchmark rules

- After `vc-relax`, run one **fixed-occupation SCF**, then the **nscf** for the
  band gap, so the gap is not computed from a smeared density.
- GaN needs `nosym=.true.` and `diagonalization='cg'` in the nscf.

## Quick start

```bash
# Print the mapping for a CIF (uses data/MPRelaxSet_c4998d92.yaml + SSSP JSON)
python scripts/show_mapping.py examples/mp-2172_AlAs.cif \
  --sssp-json /path/to/SSSP_1.3.0_PBE_efficiency.json

# Machine-readable output
python scripts/show_mapping.py examples/mp-2172_AlAs.cif \
  --sssp-json /path/to/SSSP_1.3.0_PBE_efficiency.json --json
```

Or from Python:

```python
from mprelax_config import mapping_for
m = mapping_for(structure, sssp)
```

## Layout

```
qe-mprelax-mapping/
├── SKILL.md
├── README.md
├── requirements.txt
├── config.example.yaml
├── scripts/mprelax_config.py   # pure functions
├── scripts/show_mapping.py     # CLI
├── data/MPRelaxSet_c4998d92.yaml
├── data/mapping_table.json
├── references/vasp-qe-mapping.md
├── references/deviations.md
└── examples/mp-2172_AlAs.cif
```
