# Workflow and design

## Controlled gap

`gap_controlled.py`:

1. Reads a fixed structure (`vasp_final.cif` or `initial.cif`).
2. Builds SCF + either a seekpath high-symmetry path (`--mode path`) or a dense
   uniform nscf (`--mode dense`).
3. Repeats for each pseudopotential set.
4. Compares to the WBM reference gap in `meta.json`.

Recommended baseline: `--mode dense --nscf-density 0.08` (or 0.10 for speed).

## Force matrix

`force_matrix.py`:

1. Reads a prior MPtrj run (`jobs.json` + `vasp_frame.json`).
2. Selects the worst `--top-n` frames by `force_max_err_eV_A` (or explicit
   `--sids`).
3. Builds baseline + k/ecut/smearing/pseudo/combined variants.
4. Runs SCF, parses forces from the text output, compares with VASP.

## Attribution logic

| Observation | Likely cause |
|---|---|
| Changing pseudopotential changes the gap by <0.01 eV | not a pseudopotential error |
| Dense grid much better than path | off-path CBM / path choice |
| Initial structure matches VASP better than final | reference-structure mismatch |
| k/ecut/smearing/pseudo cannot reduce force error | electron-state / PES / basin difference |
| R2SCAN opens some gaps but not others | functional + reference sensitivity |

## Failure modes

- `seekpath` import missing → install it (`requirements.txt`).
- R2SCAN fails with `libxc needed` → use the CPU libxc build
  (`qe-build/scripts/build_libxc.sh`) and PseudoDojo NC pseudos.
- Tetrahedra yields no forces → expected; use Gaussian/MV.
- No frames CSV for `force_matrix.py` → pass `--sids` with a valid
  `--mptraj-dir`.
