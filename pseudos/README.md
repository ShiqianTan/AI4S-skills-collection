# Bundled pseudopotentials

Shared by the `qe-*` benchmark skills. Paths are configured in
`_shared/config.example.yaml`.

## Contents

```
pseudos/
├── sssp/
│   ├── SSSP_1.3.0_PBE_efficiency.tar.gz     # 103 UPF files (compressed)
│   ├── SSSP_1.3.0_PBE_precision.tar.gz      # 103 UPF files (compressed)
│   ├── SSSP_1.3.0_PBE_efficiency.json       # cutoffs + filenames
│   ├── SSSP_1.3.0_PBE_precision.json
│   ├── upf/                                  # created by fetch_pseudos.py
│   └── upf_precision/                        # created by fetch_pseudos.py
└── pseudodojo/                               # 73 NC-SR PBE standard UPF files
```

SSSP 1.3.0 is the PBE efficiency/precision library from the SSSP project.
PseudoDojo NC-SR PBE standard (`nc-sr-04_pbe_standard`) is the norm-conserving
library required for the QE R2SCAN path.

## Setup

```bash
python pseudos/fetch_pseudos.py          # extract the SSSP tarballs
python pseudos/fetch_pseudos.py --force  # re-extract
```

The extracted `upf/` directories are ignored by git; only the compressed
tarballs and JSON manifests are committed. PseudoDojo UPFs are committed
uncompressed (about 18 MB).

## Which library to use

| Task | Library |
|---|---|
| PBE MP20 / MPtrj / sAlex / MatPES / WBM | SSSP efficiency |
| PBE precision sensitivity check | SSSP precision |
| R2SCAN (meta-GGA) | PseudoDojo NC-SR PBE (norm-conserving only) |

QE's R2SCAN path rejects ultrasoft/PAW pseudopotentials.

## Provenance

- SSSP 1.3.0: `SSSP_1.3.0_PBE_efficiency` / `..._precision`, including their JSON
  recommendation manifests.
- PseudoDojo: `nc-sr-04_pbe_standard_upf` (NC-SR, PBE, standard accuracy).

If you need to refresh the libraries, download the official tarballs and replace
the files here, then rerun `fetch_pseudos.py --force`.
