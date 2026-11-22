# CPU libxc build (R2SCAN)

R2SCAN is a meta-GGA. QE requires both `--with-libxc` and **norm-conserving**
pseudopotentials. The GPU NVHPC build cannot read gfortran-compiled libxc
Fortran modules, so R2SCAN uses a dedicated CPU build.

## Dependencies

```bash
conda install -n qe -c conda-forge libxc
```

Verified libxc: 7.1.2 (conda-forge).

## Build

```bash
conda activate qe
export CC=$CONDA_PREFIX/bin/mpicc
export FC=$CONDA_PREFIX/bin/mpif90
export F90=$FC
export MPIF90=$FC
export BLAS_LIBS="-L$CONDA_PREFIX/lib -lopenblas"
export LAPACK_LIBS="-L$CONDA_PREFIX/lib -lopenblas"
export SCALAPACK_LIBS="-L$CONDA_PREFIX/lib -lscalapack"
export FFT_LIBS="-L$CONDA_PREFIX/lib -lfftw3_omp -lfftw3"
./configure --prefix=.../qe-install-cpu-libxc \
  --enable-parallel --enable-openmp --with-scalapack=yes \
  --with-libxc=yes --with-libxc-prefix=$CONDA_PREFIX
python3 scripts/patch_make_inc.py --make-inc make.inc --profile cpu-libxc \
  --mpif90 $CONDA_PREFIX/bin/mpif90 --conda-prefix $CONDA_PREFIX
make -j64 all && make install
```

`patch_make_inc.py --profile cpu-libxc` adds `-D__LIBXC` to `DFLAGS` in
addition to the CPU fixes.

## Pseudopotentials

Use PseudoDojo NC-SR PBE standard (`pseudodojo/` in `qe-r2scan`). SSSP
ultrasoft/PAW pseudopotentials are not valid for this path.

## Run settings

```fortran
&SYSTEM
  input_dft = 'r2scan'
  ecutwfc = 60
  ecutrho = 240
/
```

Typical production: 60 Ry, 7x7x7, Gaussian smearing 0.02 Ry, 16 concurrent
jobs x 4 MPI x 4 OpenMP = 256 threads.

## Verification

Si R2SCAN SCF prints `Exchange-correlation = R2SCAN` and converges. GPU attempts
fail with `libxc needed for this functional` because the NVHPC build cannot load
the gfortran modules.
