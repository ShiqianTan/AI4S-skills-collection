# qe-build

Build and verify Quantum ESPRESSO for QE-vs-VASP parity benchmarks.

- CPU PBE build (`build_cpu.sh`)
- OpenACC/CUDA GPU build with NVHPC (`build_gpu.sh`)
- CPU libxc build for R2SCAN (`build_libxc.sh`)
- Si diamond smoke tests and CPU/GPU/libxc energy cross-checks (`verify.py`)

## Install

The Python verifier needs only the standard library. The QE builds need the
compilers and libraries described in `SKILL.md`. A ready conda environment:

```bash
conda create -n qe -c conda-forge python=3.11 numpy pandas matplotlib \
  pymatgen ase spglib openmpi openblas scalapack fftw hdf5 mpi4py cmake make
conda install -n qe -c conda-forge gfortran_linux-64 gxx_linux-64
```

## Build

```bash
QE_SRC=/path/to/q-e PREFIX=/path/to/qe-install-cpu bash scripts/build_cpu.sh
NVHPC_ROOT=/path/to/nvhpc-sdk/Linux_x86_64/25.7 QE_SRC=/path/to/q-e \
  PREFIX=/path/to/qe-install-gpu bash scripts/build_gpu.sh
QE_SRC=/path/to/q-e PREFIX=/path/to/qe-install-cpu-libxc bash scripts/build_libxc.sh
```

Each script runs `configure`, applies the known `make.inc` fixes with
`scripts/patch_make_inc.py`, then `make all && make install`.

## Verify

```bash
python scripts/verify.py --pw /path/to/bin/pw.x --mpirun mpirun
python scripts/verify.py --pw /path/to/bin/pw.x --mpirun mpirun --gpu 1
python scripts/verify.py --pw /path/to/bin/pw.x --r2scan
```

The verifier writes the output next to the input (in a temporary directory by
default) and checks for `JOB DONE`, the Si total energy and, for GPU runs, the
`GPU acceleration is ACTIVE` banner.

## Verified versions

See `data/versions.json`. QE 7.6 `qe-7.6-1-g61569eb48`, GCC 11.4.0,
NVHPC 25.7 / CUDA 12.9, OpenMPI 5.0.11, OpenBLAS 0.3.34, FFTW 3.3.11.

## Notes

- Build the GPU version in a **clean shell** (do not activate conda), otherwise
  conda `-march=nocona` / `LD=...-ld` leak into `make.inc`.
- `make distclean` de-initialises QE submodules; re-run
  `git submodule update --init external/devxlib external/mbd external/fox`.
- R2SCAN needs norm-conserving pseudopotentials and `--with-libxc`; see
  `references/libxc-r2scan.md`. The GPU build cannot run R2SCAN with gfortran
  libxc modules.
