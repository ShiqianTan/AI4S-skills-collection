# sAlex workflow details

## 1. Scan

`salex_scan.py` reads every `*.aselmdb` with ASE's `connect` and builds
`groups.json`: one entry per material group with frame count, element set,
atom count and force range.

## 2. Select

`salex_select.py` uses the hard-coded `SELECTED` map (20 groups × 4 frame IDs
from the study) and the `groups.json` file list to locate the rows. For each
frame it writes:

```
runs/salex/<group>/<sid>/
  structure.cif
  salex_frame.json   # sAlex energy, forces, stress, volume
  scf.in
runs/salex/jobs.json
```

Jobs set `stress_unit: eV_A3` and `stress_sign: 1`.

## 3. Run

`run_mptraj.py` (shared with the MPtrj skill) runs one static SCF per frame,
parses ASE `espresso-out`, and compares:

- forces: component MAE/RMSE/max,
- stress: after unit/sign handling,
- relative energy within each material group (`mp_id` = group id).

## 4. Plot

`plot_results.py` writes force parity, force-magnitude parity, relative-energy
parity and stress parity.

## Failure modes

- **Missing sid**: the selected ID is not in the LMDB file list from
  `groups.json`; rerun the scan or check the dataset version.
- **No stress**: the LMDB row lacks `stress`; the frame is still usable for
  forces/energies.
- **kbar conversion applied by mistake**: sAlex is already eV/Å³; ensure the job
  has `stress_unit: eV_A3`, `stress_sign: 1`.
