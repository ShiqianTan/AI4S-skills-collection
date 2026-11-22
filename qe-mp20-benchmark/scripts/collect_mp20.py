#!/usr/bin/env python3
"""Collect MP20 QE results and compare with the MP/VASP reference values."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core.parse import RY_TO_EV  # noqa: E402
from qe_bench_core.structure import read_structure  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--mp20-dir", default=None)
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--tag", default="mp20_20")
    args = ap.parse_args()

    p = get_paths(args.config)
    mp20 = Path(args.mp20_dir) if args.mp20_dir else p.work / "mp20_20"
    data_dir = Path(args.data_dir) if args.data_dir else p.skill_root / "data" / "mp20_20"
    results = p.results

    manifest = json.loads((mp20 / "manifest.json").read_text())
    status = json.loads((mp20 / "run_status.json").read_text())
    ref_energy = {
        el: data["energy_per_atom_Ry"]
        for el, data in status.get("references", {}).items()
        if data.get("status") == "ok"
    }

    rows = []
    for mp_id, meta in manifest["systems"].items():
        res = status["systems"].get(mp_id, {})
        base = {
            "mp_id": mp_id, "formula": meta["formula"], "nat": meta["nat"],
            "metal": meta.get("metal", False), "spacegroup_MP": meta["mp_spacegroup"],
            "MP_gap_eV": meta["mp_band_gap"], "MP_volume_A3": meta["cif_volume"],
            "MP_Ef_eV_atom": meta["mp_formation_energy_per_atom"],
            "formation_energy_available": meta.get("formation_energy_available", True),
        }
        if res.get("status") != "ok":
            rows.append({**base, "status": res.get("status")})
            continue
        structure = read_structure(data_dir / f"{mp_id}_{meta['formula']}.cif")
        counts = structure.composition.as_dict()
        e_comp = res["vc_relax_energy_Ry"]
        nat = sum(counts.values())
        ef_raw = None
        if base["formation_energy_available"] and all(el in ref_energy for el in counts):
            e_ref_sum = sum(counts[el] * ref_energy[el] for el in counts)
            ef_raw = (e_comp - e_ref_sum) / nat * RY_TO_EV
        rows.append({
            **base, "status": "ok",
            "QE_gap_eV": res["qe_band_gap_eV"],
            "gap_diff_eV": round(res["qe_band_gap_eV"] - meta["mp_band_gap"], 4),
            "QE_volume_A3": res["relaxed_volume_A3"],
            "volume_diff_pct": res["volume_error_pct"],
            "QE_Ef_raw_eV_atom": None if ef_raw is None else round(ef_raw, 4),
            "Ef_diff_eV_atom": None if ef_raw is None else round(
                ef_raw - meta["mp_formation_energy_per_atom"], 4),
            "QE_E_total_Ry": e_comp,
            "vc_relax_wall_s": res.get("vc_relax_wall_s"),
            "nscf_wall_s": res.get("nscf_wall_s"),
        })

    df = pd.DataFrame(rows)
    results.mkdir(parents=True, exist_ok=True)
    df.to_csv(results / f"{args.tag}_results.csv", index=False)
    ok = df[df.status == "ok"]
    ef_ok = ok[ok.Ef_diff_eV_atom.notna()]
    summary = {
        "n_systems": int(len(ok)),
        "n_formation_energy": int(len(ef_ok)),
        "gap_MAE_eV": round(float(ok.gap_diff_eV.abs().mean()), 4),
        "gap_RMSE_eV": round(float((ok.gap_diff_eV ** 2).mean() ** 0.5), 4),
        "gap_max_abs_eV": round(float(ok.gap_diff_eV.abs().max()), 4),
        "volume_MAE_pct": round(float(ok.volume_diff_pct.abs().mean()), 4),
        "volume_max_abs_pct": round(float(ok.volume_diff_pct.abs().max()), 4),
        "Ef_MAE_eV_atom": round(float(ef_ok.Ef_diff_eV_atom.abs().mean()), 4) if len(ef_ok) else None,
        "Ef_max_abs_eV_atom": round(float(ef_ok.Ef_diff_eV_atom.abs().max()), 4) if len(ef_ok) else None,
    }
    (results / f"{args.tag}_summary.json").write_text(json.dumps(summary, indent=2))
    cols = ["mp_id", "formula", "nat", "metal", "status", "MP_gap_eV", "QE_gap_eV",
            "gap_diff_eV", "MP_volume_A3", "QE_volume_A3", "volume_diff_pct",
            "MP_Ef_eV_atom", "QE_Ef_raw_eV_atom", "Ef_diff_eV_atom"]
    print(df[cols].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
