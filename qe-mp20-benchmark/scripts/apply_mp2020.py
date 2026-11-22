#!/usr/bin/env python3
"""Apply MP2020 anion/oxide corrections to QE formation energies."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core.mp2020 import correction_per_atom  # noqa: E402
from qe_bench_core.parse import RY_TO_EV  # noqa: E402
from qe_bench_core.structure import read_structure  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--tag", default="mp20_20")
    ap.add_argument("--data-dir", default=None)
    args = ap.parse_args()

    p = get_paths(args.config)
    data_dir = Path(args.data_dir) if args.data_dir else p.skill_root / "data" / "mp20_20"
    results = p.results
    df = pd.read_csv(results / f"{args.tag}_results.csv")
    rows = []
    for _, r in df.iterrows():
        if r.status != "ok" or pd.isna(r.get("QE_Ef_raw_eV_atom")):
            continue
        structure = read_structure(data_dir / f"{r.mp_id}_{r.formula}.cif")
        corr, details = correction_per_atom(structure.composition, float(r.QE_E_total_Ry) * RY_TO_EV)
        ef_corr = float(r.QE_Ef_raw_eV_atom) + corr
        rows.append({
            "mp_id": r.mp_id, "formula": r.formula,
            "QE_Ef_raw_eV_atom": r.QE_Ef_raw_eV_atom,
            "MP2020_correction_eV_atom": round(corr, 4),
            "QE_Ef_MP2020_eV_atom": round(ef_corr, 4),
            "MP_Ef_eV_atom": r.MP_Ef_eV_atom,
            "Ef_diff_raw_eV_atom": round(float(r.QE_Ef_raw_eV_atom) - float(r.MP_Ef_eV_atom), 4),
            "Ef_diff_MP2020_eV_atom": round(ef_corr - float(r.MP_Ef_eV_atom), 4),
            "corrections": json.dumps(details),
        })
    out = pd.DataFrame(rows)
    out.to_csv(results / f"{args.tag}_mp2020.csv", index=False)
    summary = {
        "n_systems": int(len(out)),
        "Ef_MAE_raw_eV_atom": round(float(out.Ef_diff_raw_eV_atom.abs().mean()), 4),
        "Ef_MAE_MP2020_eV_atom": round(float(out.Ef_diff_MP2020_eV_atom.abs().mean()), 4),
        "Ef_max_abs_raw_eV_atom": round(float(out.Ef_diff_raw_eV_atom.abs().max()), 4),
        "Ef_max_abs_MP2020_eV_atom": round(float(out.Ef_diff_MP2020_eV_atom.abs().max()), 4),
    }
    (results / f"{args.tag}_mp2020_summary.json").write_text(json.dumps(summary, indent=2))
    print(out.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
