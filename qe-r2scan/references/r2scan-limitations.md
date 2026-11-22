# R2SCAN limitations and fixes

## Requirements

- QE built with `--with-libxc` (`qe-build/scripts/build_libxc.sh`).
- **Norm-conserving** pseudopotentials; ultrasoft/PAW are rejected. Use
  PseudoDojo NC-SR PBE standard (`pseudos/pseudodojo/`).
- `input_dft = 'r2scan'` in `&SYSTEM`.
- CPU build: the NVHPC GPU build cannot load gfortran-compiled libxc modules.

## GPU is unavailable

Attempting R2SCAN with the NVHPC build fails with:

```
libxc needed for this functional
```

`nvfortran` cannot read libxc Fortran modules built by gfortran. A working GPU
R2SCAN needs nvfortran-compatible libxc bindings, which were not available.

## Recommended fixes / next steps

1. **R2SCAN-compatible pseudopotentials** (or a NC set closer to the VASP PAW
   setup).
2. **Convergence matrix** for ecut and k-points; meta-GGA is more sensitive than
   PBE.
3. **GPU support** by building libxc with nvfortran and linking it into QE.
4. **Controlled VASP reference**: match the VASP R2SCAN INCAR (PAW, cutoff,
   k-grid, smearing) before drawing quantitative conclusions.

## CPU configuration

- 60 Ry / 7×7×7 is ~2× faster than 60 Ry / 9×9×9 with ~0.01 eV/Å force
  difference.
- 16 concurrent × 4 MPI × 4 OpenMP = 256 threads was the production setting.

## When to use R2SCAN

- Qualitative trend checks.
- Testing whether a band gap opens or a relative ordering changes.
- Not as the only quantitative reference for forces/stress/formation energies.
