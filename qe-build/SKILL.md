---
name: qe-build
description: Build and verify Quantum ESPRESSO for QE-vs-VASP parity benchmarks. Covers the CPU PBE build, the OpenACC/CUDA GPU build with NVHPC, and the CPU libxc build required for R2SCAN. Includes the known configure/make fixes, Si diamond smoke tests, and energy cross-checks between builds. Use when preparing a reproducible QE 7.6 installation before any benchmark skill.
metadata:
  runtime: local
  execution_mode: cli
  cmd_policy: required
---

# Quantum ESPRESSO build and verification (qe-build)

This skill packages the build recipe used for the QE 7.6 vs VASP/MP study. It
produces three independent installs that later skills depend on:

| Build | Prefix | Purpose | Needs |
|---|---|---|---|
| CPU PBE | `qe-install-cpu` | MP20 / MPtrj / sAlex / MatPES / WBM PBE | conda GCC/OpenMPI/OpenBLAS/ScaLAPACK/FFTW |
| GPU PBE | `qe-install-gpu` | fast PBE runs (one GPU per MPI rank) | NVHPC SDK (nvfortran), CUDA |
| CPU libxc | `qe-install-cpu-libxc` | R2SCAN (meta-GGA) | conda `libxc`, norm-conserving pseudos |

All three must reproduce the same Si diamond energy:
`-16.92302587 Ry` at `ecutwfc=30 Ry`, `8x8x8`, PBE `Si_r.upf`.

## Environment

```bash
conda create -n qe -c conda-forge python=3.11 numpy pandas matplotlib \
  pymatgen ase spglib openmpi openblas scalapack fftw hdf5 mpi4py cmake make
conda install -n qe -c conda-forge gfortran_linux-64 gxx_linux-64
```

Verified toolchain: OpenMPI 5.0.11, OpenBLAS 0.3.34, FFTW 3.3.11, GCC 11.4,
QE 7.6 (`qe-7.6-1-g61569eb48`), NVHPC 25.7 / CUDA 12.9.

## Quick start

```bash
# CPU PBE
QE_SRC=~/qe/q-e PREFIX=~/qe/qe-install-cpu bash scripts/build_cpu.sh

# GPU PBE (clean shell: do NOT activate conda)
NVHPC_ROOT=~/nvhpc-sdk/Linux_x86_64/25.7 QE_SRC=~/qe/q-e \
  PREFIX=~/qe/qe-install-gpu bash scripts/build_gpu.sh

# CPU libxc (R2SCAN)
QE_SRC=~/qe/q-e PREFIX=~/qe/qe-install-cpu-libxc bash scripts/build_libxc.sh

# Verify (CPU / GPU / libxc)
python scripts/verify.py --pw ~/qe/qe-install-cpu/bin/pw.x --mpirun mpirun
python scripts/verify.py --pw ~/qe/qe-install-gpu/bin/pw.x --mpirun mpirun --gpu 1
python scripts/verify.py --pw ~/qe/qe-install-cpu-libxc/bin/pw.x --r2scan
```

Then point the shared config at the installs:

```yaml
qe:
  bin_dir: /path/to/qe-install-cpu/bin
  gpu_bin_dir: /path/to/qe-install-gpu/bin
  libxc_bin_dir: /path/to/qe-install-cpu-libxc/bin
```

## Known build fixes (already applied by the scripts)

1. **gfortran argument mismatch** in `UtilXlib/mp.f90` →
   `-fallow-argument-mismatch -fallow-invalid-boz`.
2. **Wrong linker driver** (`x86_64-conda-linux-gnu-ld`) → `LD=mpif90`.
3. **OpenMP runtime not linked** → `-fopenmp` in `LDFLAGS`.
4. **Threaded FFTW symbols** → `-lfftw3_omp -lfftw3` (CPU) or
   `-lfftw3_threads -lfftw3` (GPU, avoids conda `libgomp` vs NVHPC OpenMP).
5. **NVHPC installs to `/opt`** → silent install to a workspace directory.
6. **CUDA 13 driver vs bundled MPI** → use `comm_libs/openmpi4/bin/mpif90`.
7. **conda flags leak into GPU build** → build in a clean shell.
8. **conda `mpi.mod` shadows NVHPC** → isolated FFTW include dir.
9. **`make distclean` deinits submodules** → re-run
   `git submodule update --init external/devxlib external/mbd external/fox`.
10. **libxc Fortran modules** are gfortran-only → R2SCAN runs on the CPU build;
    the NVHPC GPU build cannot read them (`libxc needed for this functional`).

See `references/build-fixes.md`, `references/gpu-nvhpc.md` and
`references/libxc-r2scan.md` for the full failure/repair log.

## Acceptance

- CPU and GPU Si energies agree to printed precision (`-16.92302587 Ry`).
- GPU log contains `GPU acceleration is ACTIVE`.
- CPU libxc Si R2SCAN log contains `Exchange-correlation = R2SCAN` and converges.
- `pw.x`, `ph.x`, `pp.x` present in each install's `bin/`.

## Layout

```
qe-build/
├── SKILL.md
├── README.md
├── requirements.txt
├── config.example.yaml
├── scripts/        # build_cpu.sh / build_gpu.sh / build_libxc.sh / patch_make_inc.py / verify.py
├── templates/      # Si diamond SCF inputs (PBE, R2SCAN)
├── data/versions.json
├── references/     # build-fixes / gpu-nvhpc / libxc-r2scan / verification
├── pseudo/         # Si_r.upf (smoke test)
└── examples/si/    # ready-to-run smoke inputs
```
