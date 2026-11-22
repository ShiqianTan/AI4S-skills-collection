#!/usr/bin/env python3
"""Scan sAlex validation LMDBs and summarize small-cell material groups."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from ase.db import connect  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402

EXCLUDE = set("""He Ne Ar Kr Xe Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr
Sc Ti V Cr Mn Fe Co Ni Y Zr Tc La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au
Po At""".split())


def group_key(sid: str) -> str:
    parts = sid.split("_")
    return "_".join(parts[:-2]) if len(parts) >= 3 else sid


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--val", default=None, help="sAlex val/ directory of *.aselmdb")
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-nat", type=int, default=8)
    args = ap.parse_args()

    p = get_paths(args.config)
    val = Path(args.val) if args.val else p.dataset("salex")
    out = Path(args.out) if args.out else p.work / "salex"
    out.mkdir(parents=True, exist_ok=True)
    sssp = set(json.loads(p.pseudo_json("efficiency").read_text()))

    files = sorted(val.glob("*.aselmdb"))
    groups: dict[str, dict] = {}
    n_rows = 0
    for fi, path in enumerate(files):
        try:
            db = connect(str(path))
        except Exception as exc:  # noqa: BLE001
            print("skip", path.name, exc, file=sys.stderr)
            continue
        for row in db.select():
            n_rows += 1
            sid = str(row.data.get("sid", ""))
            if not sid:
                continue
            gk = group_key(sid)
            atoms = row.toatoms()
            if len(atoms) > args.max_nat:
                continue
            els = {str(e) for e in atoms.get_chemical_symbols()}
            if not els.issubset(sssp) or (els & EXCLUDE):
                continue
            forces = atoms.calc.results.get("forces")
            if forces is None:
                continue
            fmax = float(np.max(np.linalg.norm(np.asarray(forces), axis=1)))
            entry = groups.setdefault(gk, {
                "group": gk, "formula": atoms.get_chemical_formula(),
                "nat": len(atoms), "elements": sorted(els),
                "frames": [], "files": set(),
            })
            entry["frames"].append({"sid": sid, "fmax": fmax,
                                    "energy": float(row.energy), "file": path.name})
            entry["files"].add(path.name)
        if (fi + 1) % 50 == 0:
            print(f"scanned {fi+1}/{len(files)} files, {n_rows} rows, {len(groups)} groups", flush=True)

    candidates = []
    for gk, g in groups.items():
        frames = g["frames"]
        if len(frames) < 4:
            continue
        fmaxs = sorted(f["fmax"] for f in frames)
        candidates.append({
            "group": gk, "formula": g["formula"], "nat": g["nat"],
            "elements": g["elements"], "n_frames": len(frames),
            "fmax_min": round(fmaxs[0], 4), "fmax_median": round(fmaxs[len(fmaxs)//2], 4),
            "fmax_max": round(fmaxs[-1], 4), "files": sorted(g["files"]),
        })
    candidates.sort(key=lambda c: (-(c["fmax_max"] - c["fmax_min"]), c["nat"], c["formula"]))
    (out / "groups.json").write_text(json.dumps(candidates, indent=2))
    print(f"rows={n_rows} groups_total={len(groups)} candidate_groups={len(candidates)}")
    for c in candidates[:60]:
        print(c["group"], c["formula"], "nat", c["nat"], "frames", c["n_frames"],
              "F", c["fmax_min"], c["fmax_median"], c["fmax_max"], c["elements"])


if __name__ == "__main__":
    main()
