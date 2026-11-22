# AI4S-skills-collection

AI for Science (AI4S) 计算 skills 集合。

## Skills

- [qe-powder-raman](./qe-powder-raman/README.md) — 基于 Quantum ESPRESSO 的固体粉末拉曼（一级非共振谐波）计算 skill

## QE VASP-parity benchmark skills

Methodology packaged from a Quantum ESPRESSO 7.6 vs VASP / Materials Project
reproduction study. The result tables and error attribution are reproduced in
each skill's `references/` directory.

**Full report (PDF):** [reports/QE-vs-VASP-parity-report.pdf](./reports/QE-vs-VASP-parity-report.pdf)

| Skill | What it does |
|---|---|
| [qe-build](./qe-build/README.md) | Build/verify QE (CPU PBE, GPU OpenACC, CPU libxc/R2SCAN) with known fixes and Si smoke tests |
| [qe-mprelax-mapping](./qe-mprelax-mapping/README.md) | Map the Materials Project MPRelaxSet (c4998d92) to QE parameters |
| [qe-mp20-benchmark](./qe-mp20-benchmark/README.md) | MP20 geometry, band gap, formation energy, MP2020 corrections, convergence checks |
| [qe-mptraj-benchmark](./qe-mptraj-benchmark/README.md) | MPtrj forces/stress/relative energies and force-error-vs-distortion |
| [qe-salex-benchmark](./qe-salex-benchmark/README.md) | sAlex forces/stress/relative energies (ASE stress convention) |
| [qe-matpes-benchmark](./qe-matpes-benchmark/README.md) | MatPES-PBE forces/stress/relative energies (kbar conversion) |
| [qe-wbm-basin](./qe-wbm-basin/README.md) | WBM IS2RE optimization, formation energies, basin recovery and high-distortion rescue |
| [qe-convergence-attribution](./qe-convergence-attribution/README.md) | Controlled band-gap and high-force k/ecut/smearing/pseudo attribution |
| [qe-r2scan](./qe-r2scan/README.md) | MatPES-R2SCAN with libxc + PseudoDojo NC, CPU configuration benchmark |

## Shared setup

Shared Python utilities live in [`_shared/qe_bench_core`](./_shared/README.md).
Each skill reads one `config.yaml` for QE binaries, datasets and
pseudopotentials (see [`_shared/config.example.yaml`](./_shared/config.example.yaml)):

```yaml
qe:
  bin_dir: /path/to/qe-install-cpu/bin
  gpu_bin_dir: /path/to/qe-install-gpu/bin
  libxc_bin_dir: /path/to/qe-install-cpu-libxc/bin
  mpirun: mpirun
  nvhpc_root: /path/to/nvhpc-sdk/Linux_x86_64/25.7
  conda_lib: /path/to/envs/qe/lib
pseudo:
  sssp_efficiency: pseudos/sssp/upf
  sssp_precision: pseudos/sssp/upf_precision
  pseudodojo: pseudos/pseudodojo
datasets:
  mp20: /path/to/MP20_dataset
  mptraj: /path/to/MPtrj_2022.9_full.json
  salex: /path/to/sAlex/val
  matpes: /path/to/matpes
  wbm: /path/to/matbench-discovery/data/wbm
  wbm_is2re: /path/to/WBM_IS2RE.aselmdb
```

Path resolution priority: `--config` > `QE_BENCH_CONFIG` > `QE_*` environment
variables > `config.yaml` > built-in defaults.

## Pseudopotentials

[`pseudos/`](./pseudos/README.md) bundles the SSSP 1.3.0 efficiency/precision
tarballs + JSON manifests and the PseudoDojo NC-SR PBE UPF files. Unpack the SSSP
archives once after cloning:

```bash
python pseudos/fetch_pseudos.py
```

## Typical workflow

```bash
# 1) build QE (CPU/GPU/libxc) and verify Si
bash qe-build/scripts/build_cpu.sh
python qe-build/scripts/verify.py --pw /path/to/bin/pw.x

# 2) near-equilibrium forces/stress
python qe-matpes-benchmark/scripts/matpes_select.py --config config.yaml
python qe-matpes-benchmark/scripts/run_mptraj.py --config config.yaml --mptraj-dir runs/matpes --tag matpes

# 3) structures/formation energies
python qe-mp20-benchmark/scripts/mp20_prepare.py --config config.yaml
python qe-mp20-benchmark/scripts/run_mp20.py --config config.yaml

# 4) structure optimization and basin checks
python qe-wbm-basin/scripts/wbm_select.py --config config.yaml
python qe-wbm-basin/scripts/run_wbm.py --config config.yaml
python qe-wbm-basin/scripts/basin_recovery.py --config config.yaml

# 5) error attribution
python qe-convergence-attribution/scripts/gap_controlled.py --config config.yaml --mode dense
python qe-convergence-attribution/scripts/force_matrix.py --config config.yaml
```
