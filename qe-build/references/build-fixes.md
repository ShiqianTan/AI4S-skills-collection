# CPU build: failures and fixes

Verified with QE 7.6 (`qe-7.6-1-g61569eb48`), conda GCC/gfortran 11.4.0,
OpenMPI 5.0.11, OpenBLAS 0.3.34, FFTW 3.3.11.

## 1. gfortran type mismatch in `UtilXlib/mp.f90`

```
Error: Type mismatch between actual argument ... (REAL(8)/COMPLEX(8))
```

QE 7.6 MPI calls are too loose for modern gfortran. Append to `FFLAGS`:

```
-fallow-argument-mismatch -fallow-invalid-boz
```

## 2. Wrong linker driver

Activating conda sets `LD=x86_64-conda-linux-gnu-ld`, which rejects `-Wl,...`
options. Set `LD` to the MPI compiler wrapper:

```
LD = /path/to/conda/bin/mpif90
```

## 3. OpenMP runtime not linked

Append `-fopenmp` to `LDFLAGS`.

## 4. Threaded FFTW symbols missing

Link the OpenMP FFTW library before the serial one:

```
FFT_LIBS = -L$CONDA_PREFIX/lib -lfftw3_omp -lfftw3
```

## Working configuration

```bash
export CC=$CONDA_PREFIX/bin/mpicc
export FC=$CONDA_PREFIX/bin/mpif90
export F90=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gfortran
export MPIF90=$CONDA_PREFIX/bin/mpif90
export BLAS_LIBS="-L$CONDA_PREFIX/lib -lopenblas"
export LAPACK_LIBS="-L$CONDA_PREFIX/lib -lopenblas"
export SCALAPACK_LIBS="-L$CONDA_PREFIX/lib -lscalapack"
export FFT_LIBS="-L$CONDA_PREFIX/lib -lfftw3_omp -lfftw3"
./configure --prefix=.../qe-install-cpu \
  --enable-parallel --enable-openmp --with-scalapack=yes
python3 scripts/patch_make_inc.py --make-inc make.inc --profile cpu \
  --mpif90 $CONDA_PREFIX/bin/mpif90 --conda-prefix $CONDA_PREFIX
make -j64 all && make install
```

Result: `EXIT=0`; `pw.x`, `cp.x`, `ph.x`, `pp.x`, `bands.x` installed.
