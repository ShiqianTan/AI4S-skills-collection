# MP20 workflow details

## Stage 1: prepare

`mp20_prepare.py` reads the bundled CIFs and the MP reference manifest and
writes one directory per system:

```
runs/mp20_20/mp-2172_AlAs/
  vc_relax.in
runs/mp20_20/references/Al/vc_relax.in
runs/mp20_20/manifest.json
```

`sssp` mode uses the cutoffs stored in the manifest (SSSP efficiency, k =
0.25/0.20/0.15 Å⁻¹). `mprelax` mode recomputes with the MPRelaxSet mapping.

## Stage 2: run

`run_mp20.py`:

1. `vc-relax` from the MP CIF (BFGS ions and cell).
2. For MPRelaxSet-like non-metals: a fixed-occupation SCF on the relaxed
   structure.
3. Dense `nscf` for the band gap.
4. Elemental references: `vc-relax` with smearing.

State is stored in `run_status.json`. Use `--force` to rerun and `--skip-refs`
to reuse existing references.

## Stage 3: collect / correct / plot

- `collect_mp20.py` reads `run_status.json`, builds `{tag}_results.csv`, and
  computes gap RMSE/MAE, volume MAE and raw formation energy.
- `apply_mp2020.py` adds MP2020 corrections using the collected total energy and
  the CIF composition.
- `plot_results.py` writes parity plots from the CSV/JSON.

## Accuracy expectations

| Property | Typical | Notes |
|---|---|---|
| volume | 0.2-0.5 % | same basin |
| band gap | 0.17 eV MAE | reference-definition limited |
| MP2020 Ef | 13 meV/atom | good |
| raw Ef | 170 meV/atom | pseudopotential/reference offset |

## Common failures

- **Gap is 0**: metallic system or smeared nscf; check `occupations` and whether
  the CBM/VBM lines are present.
- **Formation energy missing**: an element has no QE reference in
  `manifest["references"]` or the reference run failed.
- **All systems identical**: the run reused `vc_relax.out`; use `--force`.
- **GPU OOM or wrong device**: set `--devices` and check `config.yaml`
  `nvhpc_root`/`conda_lib`.
