# qe-mprelax-mapping

Convert the Materials Project `MPRelaxSet` (pymatgen `c4998d92`) into Quantum
ESPRESSO parameters, and document the intentional deviations.

## Install

```bash
pip install -r requirements.txt
```

Requires `pymatgen` (for `MPRelaxSet.kpoints`) and PyYAML.

## Usage

```bash
python scripts/show_mapping.py examples/mp-2172_AlAs.cif \
  --sssp-json /path/to/SSSP_1.3.0_PBE_efficiency.json
python scripts/show_mapping.py examples/mp-2172_AlAs.cif \
  --sssp-json /path/to/SSSP_1.3.0_PBE_efficiency.json --json
```

Output example (AlAs):

```
system: AlAs  natoms: 2
ecutwfc (Ry): 38.2     # max(520 eV, SSSP Al/As)
ecutrho (Ry): 305.8
conv_thr (Ry): 7.350e-06
smearing: gaussian  degauss (Ry): 0.00367493
k-grid: (6, 6, 6) shift=(0, 0, 0)
electron_maxstep: 100  nstep: 99
```

Exact numbers depend on the SSSP manifest and the structure's reciprocal
density.

## Library

`mprelax_config.py` exposes:

- `mp_kpoints(structure)` → `(grid, shift)`
- `mp_cutoffs(structure, sssp)` → `(ecutwfc, ecutrho)`
- `mp_electrons(natoms)` → `(conv_thr, degauss)`
- `mapping_for(structure, sssp)` → full dict

## Deviations

See `references/vasp-qe-mapping.md` and `references/deviations.md`. The two
important ones are `ISMEAR=-5` → Gaussian smearing (QE tetrahedra is not safe
for forces) and `ISPIN=2` → `nspin=1` for the non-magnetic benchmark subset.
