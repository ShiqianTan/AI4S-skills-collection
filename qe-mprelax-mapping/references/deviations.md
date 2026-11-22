# Intentional deviations from MPRelaxSet

MPRelaxSet cannot be reproduced exactly in QE. These deviations are deliberate
and should be reported alongside any QE-vs-VASP comparison.

## `ISMEAR = -5` → Gaussian smearing

VASP uses tetrahedron integration with Blöchl corrections. QE tetrahedra
(`occupations='tetrahedra'`) cannot provide reliable forces (QE warns, some
frames print no forces, and SCF can fail). The benchmark therefore uses Gaussian
smearing `degauss=0.05 eV`.

Consequence: small energy differences from the reference remain even when the
geometry matches. `SIGMA`/`degauss` is a numerical setting, not physics.

## `ISPIN = 2` → `nspin = 1`

The selected MP20/MPtrj/MatPES/WBM benchmark subset is non-magnetic (main-group,
closed-shell). Setting `nspin=1` saves time and avoids spurious magnetic
solutions. Magnetic systems must use `nspin=2` with sensible initial moments.

## `LDAU = true` → no Hubbard U

MPRelaxSet enables U for Co/Cr/Fe/Mn/Mo/Ni/V/W in F/O environments. None of the
selected systems contains such a pair, so `lda_plus_u` is not enabled. Any
extension to transition-metal oxides/fluorides must add `U` and match the VASP
`LDAUU/LDAUL/LDAUTYPE` values.

## `ALGO = Fast` → Davidson (CG fallback)

Davidson is robust for the insulating systems in the benchmark. Use CG
(`diagonalization='cg'`) with `nosym=.true.` for GaN and other difficult cases.

## Post-processing keys

`LWAVE`, `LORBIT`, `LDAUPRINT` are not mapped; they only affect files that the
benchmark does not consume.

## What is *not* a deviation

- `ENCUT=520 eV` is honoured as a lower bound.
- The `reciprocal_density=64` k-mesh is taken directly from pymatgen.
- `NSW=99`, `IBRION=2`, `ISIF=3` map to `vc-relax` + BFGS.
- `EDIFF_PER_ATOM` and `SIGMA` are converted with the correct unit factors.
