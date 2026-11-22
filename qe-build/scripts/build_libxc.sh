#!/usr/bin/env bash
# Build Quantum ESPRESSO with libxc for R2SCAN (CPU only).
#
# Required environment:
#   QE_SRC  path to the q-e source tree
#   PREFIX  install prefix (e.g. /path/to/qe-install-cpu-libxc)
# Optional:
#   JOBS    parallel make jobs (default: nproc)
#   CONDA_PREFIX  active conda env providing OpenMPI/libraries and libxc
set -euo pipefail

: "${QE_SRC:?set QE_SRC to the q-e source directory}"
: "${PREFIX:?set PREFIX to the install prefix}"
: "${CONDA_PREFIX:?activate the conda env that provides libxc}"
JOBS="${JOBS:-$(nproc)}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX_BIN="${CONDA_PREFIX}/bin"

export CC="${CC:-${PREFIX_BIN}/mpicc}"
export FC="${FC:-${PREFIX_BIN}/mpif90}"
export F90="${F90:-${FC}}"
export MPIF90="${MPIF90:-${FC}}"
export BLAS_LIBS="${BLAS_LIBS:--L${CONDA_PREFIX}/lib -lopenblas}"
export LAPACK_LIBS="${LAPACK_LIBS:--L${CONDA_PREFIX}/lib -lopenblas}"
export SCALAPACK_LIBS="${SCALAPACK_LIBS:--L${CONDA_PREFIX}/lib -lscalapack}"
export FFT_LIBS="${FFT_LIBS:--L${CONDA_PREFIX}/lib -lfftw3_omp -lfftw3}"
export FFTW_INCLUDE="${FFTW_INCLUDE:-${CONDA_PREFIX}/include}"

cd "${QE_SRC}"
./configure --prefix="${PREFIX}" \
  --enable-parallel --enable-openmp --with-scalapack=yes \
  --with-libxc=yes --with-libxc-prefix="${CONDA_PREFIX}"

python3 "${HERE}/patch_make_inc.py" \
  --make-inc make.inc --profile cpu-libxc \
  --mpif90 "${MPIF90}" --conda-prefix "${CONDA_PREFIX}"

make -j"${JOBS}" all
make install

echo "CPU libxc build finished: ${PREFIX}/bin/pw.x"
echo "R2SCAN requires norm-conserving pseudopotentials (e.g. PseudoDojo NC-SR)."
