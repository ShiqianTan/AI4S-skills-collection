# VASP MPRelaxSet → QE mapping details

Reference: `data/MPRelaxSet_c4998d92.yaml`, from pymatgen `MPRelaxSet`
(`c4998d92`). These settings are what the Materials Project used for the
reference energies, so matching them as closely as QE allows is the point of
this skill.

## Cutoffs

```
ENCUT = 520 eV
ecutwfc_Ry = max(520 / 13.605693122994, max_elements SSSP.cutoff_wfc)
ecutrho_Ry = max(8 * ecutwfc_Ry, max_elements SSSP.cutoff_rho)
```

The 520 eV floor keeps QE from using a lower plane-wave cutoff than VASP, while
the SSSP recommendations prevent under-converged pseudopotentials. The `8×`
density ratio is appropriate for norm-conserving pseudopotentials; ultrasoft/PAW
SSSP entries carry their own `cutoff_rho`.

## SCF convergence and smearing

```
conv_thr_Ry = max(1e-8, 5e-5 eV/atom * natoms / 13.605693122994)
degauss_Ry  = 0.05 eV / 13.605693122994
smearing    = Gaussian
electron_maxstep = 100
mixing_beta = 0.7
```

MP uses `EDIFF_PER_ATOM=5e-5 eV`; QE's `conv_thr` is a total energy threshold
in Ry, hence the `natoms` factor. `SIGMA=0.05 eV` maps directly to `degauss`.

## k-points

`MPRelaxSet(structure).kpoints` uses `reciprocal_density=64`. Read the generated
grid and shift from pymatgen rather than re-deriving it, so the QE mesh matches
the MP mesh. Gamma-centred meshes use shift `(0,0,0)`; Monkhorst-Pack uses
`(1,1,1)`.

## Relaxation

```
calculation = 'vc-relax'     # ISIF=3
ion_dynamics = 'bfgs'        # IBRION=2
cell_dynamics = 'bfgs'
nstep = 99                   # NSW=99
forc_conv_thr = 5e-4 Ry/bohr
etot_conv_thr = 1e-5 Ry
```

## Electronic solver

MP uses `ALGO=Fast`. QE Davidson is the closest robust choice. For problematic
systems (e.g. GaN), use CG (`diagonalization='cg'`) and `nosym=.true.`.

## Post-relaxation band gap

Do **not** compute the gap from the smeared vc-relax density. Run one
fixed-occupation SCF, then a dense nscf, and take the gap from that.

## Units

- 1 Ry = 13.605693122994 eV
- 1 bohr = 0.529177210903 Å
- VASP stress is in kbar; 1 kbar = 6.241509×10⁻⁴ eV/Å³
- QE stress has the opposite sign to the VASP/MPtrj convention.
