# GPU build (OpenACC/CUDA, NVHPC)

QE's GPU port needs `nvfortran` (NVHPC SDK), not only `nvcc`.

Verified: NVHPC 25.7, bundled CUDA 12.9, default target `cc120`.
The SDK was installed under the workspace, not `/opt`.

## Installation

```bash
export NVHPC_SILENT=true
export NVHPC_INSTALL_DIR=/path/to/nvhpc-sdk
export NVHPC_INSTALL_TYPE=single
# run the tarball's install script
```

## Failures and fixes

1. **Installer defaults to `/opt`** and fills a small root overlay
   (`Disk quota exceeded`). Use `NVHPC_SILENT=true` and
   `NVHPC_INSTALL_DIR=<workspace>/nvhpc-sdk`.
2. **MPI trampoline rejects the CUDA driver.** `comm_libs/mpi/bin/mpif90`
   expects the bundled CUDA 12.9, but the driver reports CUDA 13.0. Use
   `comm_libs/openmpi4/bin/mpif90`.
3. **conda compiler flags leak into `make.inc`** (`-march=nocona`,
   `LD=...-ld`). Build in a clean shell: do not `conda activate`; unset
   `CC/FC/F90/CFLAGS/FFLAGS/LDFLAGS/LD`.
4. **conda `mpi.mod` shadows NVHPC's MPI module.**
   `NVFORTRAN-F-0004-Corrupt or Old Module file .../envs/qe/include/mpi.mod`.
   Copy only the FFTW headers to an isolated dir and use it as
   `FFTW_INCLUDE`; never put the conda `include` dir on the compiler path.
5. **`make distclean` runs `git submodule deinit --all --force`.** Re-init:
   `git submodule update --init external/devxlib external/mbd external/fox`.
6. **`libgomp: TODO` at run time.** conda `libfftw3_omp` pulls in conda
   `libgomp`, conflicting with the NVHPC OpenMP runtime. Link
   `-lfftw3_threads -lfftw3` instead.
7. **No rpath for FFTW.** At run time export
   `LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`.

## Working configuration

```bash
export NV=/path/to/nvhpc-sdk/Linux_x86_64/25.7
export PATH=$NV/compilers/bin:$NV/comm_libs/openmpi4/bin:$NV/cuda/12.9/bin:$PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$NV/compilers/lib:$NV/comm_libs/openmpi4/lib:$NV/cuda/12.9/lib64:$NV/math_libs/12.9/targets/x86_64-linux/lib
export CC=$NV/comm_libs/openmpi4/bin/mpicc
export FC=$NV/comm_libs/openmpi4/bin/mpif90
export F90=$NV/compilers/bin/nvfortran
export MPIF90=$NV/comm_libs/openmpi4/bin/mpif90
export BLAS_LIBS="-L$NV/compilers/lib -lblas"
export LAPACK_LIBS="-L$NV/compilers/lib -llapack"
export FFT_LIBS="-L$CONDA_PREFIX/lib -lfftw3_threads -lfftw3"
export FFTW_INCLUDE=/path/to/qe-deps/fftw/include
./configure --prefix=.../qe-install-gpu \
  --enable-parallel --enable-openmp --with-scalapack=no \
  --with-gpu=cuda --with-cuda-cc=120 --with-cuda-runtime=12.9
make -j32 all && make install
```

## Verification

GPU Si diamond reproduces the CPU energy to printed precision
(`-16.92302587 Ry`) and prints
`GPU acceleration is ACTIVE. 1 visible GPUs per MPI rank`.

The GPU build **cannot** run R2SCAN with gfortran-built libxc modules; use the
CPU libxc build for meta-GGA.
