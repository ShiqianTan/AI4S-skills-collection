#!/usr/bin/env bash
# Build Quantum ESPRESSO with the conda GCC/OpenMPI toolchain (PBE, CPU).
#
# Required environment:
#   QE_SRC  path to the q-e source tree (checked out at qe-7.6-1)
#   PREFIX  install prefix (e.g. /path/to/qe-install-cpu)
# Optional:
#   JOBS    parallel make jobs (default: nproc)
#   CONDA_PREFIX  active conda env providing OpenMPI/OpenBLAS/ScaLAPACK/FFTW
set -euo pipefail

: "${QE_SRC:?set QE_SRC to the q-e source directory}"
: "${PREFIX:?set PREFIX to the install prefix}"
JOBS="${JOBS:-$(nproc)}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "warning: CONDA_PREFIX is empty; relying on mpicc/mpif90 from PATH" >&2
fi
PREFIX_BIN="${CONDA_PREFIX:-/usr}/bin"

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
  --enable-parallel --enable-openmp --with-scalapack=yes

python3 "${HERE}/patch_make_inc.py" \
  --make-inc make.inc --profile cpu \
  --mpif90 "${MPIF90}" --conda-prefix "${CONDA_PREFIX:-}"

make -j"${JOBS}" all
make install

echo "CPU PBE build finished: ${PREFIX}/bin/pw.x"
