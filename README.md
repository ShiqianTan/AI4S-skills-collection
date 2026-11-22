# AI4S-skills-collection

AI for Science (AI4S) 计算 skills 集合。

## Skills

- [qe-powder-raman](./qe-powder-raman/README.md) — 基于 Quantum ESPRESSO 的固体粉末拉曼（一级非共振谐波）计算 skill

## QE VASP-parity benchmark skills

Methodology packaged from the Quantum ESPRESSO 7.6 vs VASP/Materials Project
reproduction study. Shared code lives in [`_shared/qe_bench_core`](./_shared/README.md);
each skill reads a single `config.yaml` for QE binaries, datasets and
pseudopotentials (see [`_shared/config.example.yaml`](./_shared/config.example.yaml)).

<!-- SKILL_LIST_START -->
- [qe-build](./qe-build/README.md) — build and verify QE (CPU PBE, GPU OpenACC, CPU libxc/R2SCAN) with known fixes and Si smoke tests
- [qe-mprelax-mapping](./qe-mprelax-mapping/README.md) — map the Materials Project MPRelaxSet to QE parameters and document deviations
- [qe-mp20-benchmark](./qe-mp20-benchmark/README.md) — MP20 geometry/band-gap/formation-energy parity, MP2020 corrections, convergence checks
- [qe-mptraj-benchmark](./qe-mptraj-benchmark/README.md) — MPtrj forces/stress/relative-energy parity and force-error-vs-distortion analysis
- [qe-salex-benchmark](./qe-salex-benchmark/README.md) — sAlex forces/stress/relative-energy validation (ASE stress convention)
- [qe-matpes-benchmark](./qe-matpes-benchmark/README.md) — MatPES-PBE forces/stress/relative-energy validation with kbar conversion
- [qe-wbm-basin](./qe-wbm-basin/README.md) — WBM IS2RE structure optimization, MP2020 formation energies, and basin-recovery diagnostics
<!-- SKILL_LIST_END -->
