#!/usr/bin/env python3
"""Compute QE formation energies for the WBM validation and summarize."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core.mp2020 import correction_per_atom  # noqa: E402
from qe_bench_core.structure import read_structure  # noqa: E402

RY_TO_EV = 13.605693122994


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--frames", default=None, help="wbm_frames.json")
    ap.add_argument("--ref-run-status", default=None,
                    help="run_status.json with QE elemental reference energies")
    args = ap.parse_args()

    p = get_paths(args.config)
    frames_path = Path(args.frames) if args.frames else p.results / "wbm_frames.json"
    ref_path = Path(args.ref_run_status) if args.ref_run_status else (
        p.work / "mp20_20_mprelax" / "run_status.json"
    )
    if not ref_path.is_file():
        ref_path = p.work / "mp20_20" / "run_status.json"
    refs = json.loads(ref_path.read_text())["references"]
    ref_energy = {el: d["energy_per_atom_Ry"] * RY_TO_EV
                  for el, d in refs.items() if d.get("status") == "ok"}

    recs = [r for r in json.loads(frames_path.read_text()) if r.get("status") == "ok"]
    rows = []
    for r in recs:
        d = Path(r["dir"])
        struct = read_structure(d / "qe_relaxed.cif")
        counts = struct.composition.as_dict()
        if not all(el in ref_energy for el in counts):
            continue
        e_comp = r["qe_final_energy_eV"]
        nat = sum(counts.values())
        ef_raw = (e_comp - sum(counts[el] * ref_energy[el] for el in counts)) / nat
        corr, details = correction_per_atom(struct.composition, e_comp)
        ef_mp2020 = ef_raw + corr
        rows.append({
            "material_id": r["material_id"], "formula": r.get("formula"),
            "n_sites": r.get("n_sites"),
            "qe_Ef_raw_eV_atom": round(ef_raw, 4),
            "mp2020_correction_eV_atom": round(corr, 4),
            "qe_Ef_MP2020_eV_atom": round(ef_mp2020, 4),
            "wbm_Ef_MP2020_eV_atom": r["vasp_e_form_mp2020_eV_atom"],
            "Ef_diff_eV_atom": round(ef_mp2020 - r["vasp_e_form_mp2020_eV_atom"], 4),
            "volume_diff_pct": r["volume_diff_pct"],
            "structure_rms": r.get("structure_rms"),
            "structure_match": r.get("structure_match"),
            "qe_force_max_at_vasp_final_eV_A": r.get("qe_force_max_at_vasp_final_eV_A"),
            "qe_vaspfinal_energy_diff_meV_atom": r.get("qe_vaspfinal_energy_diff_meV_atom"),
            "qe_relaxation_energy_eV_atom": r.get("qe_relaxation_energy_eV_atom"),
            "corrections": json.dumps(details),
        })
    df = pd.DataFrame(rows)
    p.results.mkdir(parents=True, exist_ok=True)
    df.to_csv(p.results / "wbm_formation_energies.csv", index=False)
    summary = {
        "n_systems": int(len(df)),
        "Ef_MAE_eV_atom": round(float(df.Ef_diff_eV_atom.abs().mean()), 4),
        "Ef_max_abs_eV_atom": round(float(df.Ef_diff_eV_atom.abs().max()), 4),
        "Ef_RMSE_eV_atom": round(float((df.Ef_diff_eV_atom ** 2).mean() ** 0.5), 4),
    }
    (p.results / "wbm_formation_summary.json").write_text(json.dumps(summary, indent=2))
    print(df.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
