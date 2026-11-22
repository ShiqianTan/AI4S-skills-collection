# qe-convergence-attribution

Attribute QE-vs-VASP band-gap and high-force disagreements to controllable
factors.

## Install

```bash
pip install -r requirements.txt
python ../pseudos/fetch_pseudos.py
```

Needs `seekpath` for `--mode path` and `pymatgen`, `ase`, `numpy`, `pandas`,
`matplotlib`.

## Controlled gap

```bash
cp config.example.yaml config.yaml
python scripts/gap_controlled.py --config config.yaml --wbm-dir ../qe-wbm-basin/examples \
  --mode dense --pseudo-sets sssp_eff,sssp_prec,pdojo --tag gap_dense
python scripts/gap_controlled.py --config config.yaml --wbm-dir ../qe-wbm-basin/examples \
  --mode path --pseudo-sets sssp_eff,sssp_prec,pdojo --tag gap_path
```

Options: `--start final|initial`, `--functional pbe|pbesol|r2scan`,
`--nscf-density`, `--scf-density`, `--ids`, `--prepare-only`,
`--gpu-mode gpu|cpu`.

## High-force matrix

```bash
python scripts/force_matrix.py --config config.yaml --mptraj-dir runs/mptraj_20 \
  --top-n 5 --tag force_matrix
```

Variants: baseline, k 0.15-0.30, ecut 0.9/1.2/1.5×, MV/Gaussian/tetrahedra,
SSSP precision, PseudoDojo, combined.

## Scripts

| Script | Purpose |
|---|---|
| `gap_controlled.py` | fixed-structure gap: path vs dense, pseudo sets, initial/final, functional |
| `force_matrix.py` | one-factor-at-a-time k/ecut/smearing/pseudo force matrix |
| `plot_results.py` | gap MAE by pseudo, force MAE by variant |

## Reference results

See `references/bandgap-controlled.md` and `references/high-force-matrix.md`.

- Pseudopotential choice moves the gap MAE by <0.002 eV.
- Dense k-grids roughly halve the gap error relative to a high-symmetry path.
- k/ecut/smearing/pseudo change high-force errors by only 0.4-2.6 %.
