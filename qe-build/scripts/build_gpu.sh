#!/usr/bin/env bash
# Build Quantum ESPRESSO with OpenACC/CUDA using NVHPC (PBE, GPU).
#
# Required environment:
#   QE_SRC      path to the q-e source tree
#   PREFIX      install prefix (e.g. /path/to/qe-install-gpu)
#   NVHPC_ROOT  NVHPC SDK root (e.g. .../nvhpc-sdk/Linux_x86_64/25.7)
# Optional:
#   JOBS        parallel make jobs (default: nproc)
#   FFTW_INCLUDE  isolated FFTW include dir (default: QE_SRC/../qe-deps/fftw/include)
#   CONDA_LIB   conda lib dir holding libfftw3_threads.so
#
# IMPORTANT: run this in a clean shell. Do NOT `conda activate`, otherwise
# conda CFLAGS/LD leak into make.inc and break the NVHPC link.
set -euo pipefail

: "${QE_SRC:?set QE_SRC to the q-e source directory}"
: "${PREFIX:?set PREFIX to the install prefix}"
: "${NVHPC_ROOT:?set NVHPC_ROOT to the NVHPC SDK root}"
JOBS="${JOBS:-$(nproc)}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NV="${NVHPC_ROOT}"
CONDA_LIB="${CONDA_LIB:-${CONDA_PREFIX:-}/lib}"
FFTW_INCLUDE="${FFTW_INCLUDE:-$(cd "${QE_SRC}/.." && pwd)/qe-deps/fftw/include}"

# Remove conda compiler leakage.
unset CC FC F90 MPIF90 CFLAGS FFLAGS F90FLAGS LDFLAGS LD 2>/dev/null || true

export PATH="${NV}/compilers/bin:${NV}/comm_libs/openmpi4/bin:${NV}/cuda/12.9/bin:${PATH}"
export LD_LIBRARY_PATH="${CONDA_LIB}:${NV}/compilers/lib:${NV}/comm_libs/openmpi4/lib:${NV}/cuda/12.9/lib64:${NV}/math_libs/12.9/targets/x86_64-linux/lib:${LD_LIBRARY_PATH:-}"
export NVCOMPILERS="${NV}/compilers"
export CUDA_HOME="${NV}/cuda/12.9"

export CC="${NV}/comm_libs/openmpi4/bin/mpicc"
export FC="${NV}/comm_libs/openmpi4/bin/mpif90"
export F90="${NV}/compilers/bin/nvfortran"
export MPIF90="${NV}/comm_libs/openmpi4/bin/mpif90"
export BLAS_LIBS="-L${NV}/compilers/lib -lblas"
export LAPACK_LIBS="-L${NV}/compilers/lib -llapack"
export FFT_LIBS="-L${CONDA_LIB} -lfftw3_threads -lfftw3"

cd "${QE_SRC}"
# `make distclean` de-inits submodules; make sure they are present.
if [[ ! -d external/devxlib/src ]]; then
  git submodule update --init external/devxlib external/mbd external/fox
fi

./configure --prefix="${PREFIX}" \
  --enable-parallel --enable-openmp --with-scalapack=no \
  --with-gpu=cuda --with-cuda-cc=120 --with-cuda-runtime=12.9

python3 "${HERE}/patch_make_inc.py" \
  --make-inc make.inc --profile gpu \
  --mpif90 "${MPIF90}" --fft-include "${FFTW_INCLUDE}" --conda-prefix "${CONDA_LIB%/lib}"

make -j"${JOBS}" all
make install

echo "GPU PBE build finished: ${PREFIX}/bin/pw.x"
echo "At run time add ${CONDA_LIB} to LD_LIBRARY_PATH (FFTW is not rpath'd)."
