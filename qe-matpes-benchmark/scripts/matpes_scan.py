#!/usr/bin/env python3
"""Scan MatPES-PBE and summarize small-cell, non-magnetic groups."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ijson  # noqa: E402
import numpy as np  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402

EXCLUDE = set("""He Ne Ar Kr Xe Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr
Sc Ti V Cr Mn Fe Co Ni Y Zr Tc La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au
Po At""".split())


def group_key(matpes_id: str) -> str:
    return matpes_id.rsplit("_", 1)[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--json", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-nat", type=int, default=8)
    ap.add_argument("--max-mag", type=float, default=0.5)
    args = ap.parse_args()

    p = get_paths(args.config)
    json_path = Path(args.json) if args.json else p.dataset("matpes") / "MatPES-PBE-2025.1.json"
    out = Path(args.out) if args.out else p.work / "matpes"
    out.mkdir(parents=True, exist_ok=True)
    sssp = set(json.loads(p.pseudo_json("efficiency").read_text()))

    groups: dict[str, dict] = {}
    n = 0
    t0 = time.time()
    with json_path.open("rb") as f:
        for rec in ijson.items(f, "item"):
            n += 1
            mid = str(rec.get("matpes_id", ""))
            if not mid:
                continue
            gk = group_key(mid)
            sites = rec.get("structure", {}).get("sites", [])
            nat = len(sites)
            if nat > args.max_nat or nat == 0:
                continue
            els = {str(s.get("species", [{}])[0].get("element")) for s in sites}
            if not els.issubset(sssp) or (els & EXCLUDE):
                continue
            mags = [abs(float(s.get("properties", {}).get("magmom", 0.0))) for s in sites]
            if max(mags) > args.max_mag:
                continue
            forces = rec.get("forces")
            if not forces:
                continue
            fmax = float(np.max(np.linalg.norm(np.asarray(forces, float), axis=1)))
            entry = groups.setdefault(gk, {
                "group": gk, "nat": nat, "elements": sorted(els),
                "n_frames": 0, "fmax_min": fmax, "fmax_max": fmax,
            })
            entry["n_frames"] += 1
            entry["fmax_min"] = min(entry["fmax_min"], fmax)
            entry["fmax_max"] = max(entry["fmax_max"], fmax)
    candidates = sorted(groups.values(), key=lambda g: (-(g["fmax_max"] - g["fmax_min"]), g["nat"]))
    (out / "groups.json").write_text(json.dumps(candidates, indent=2))
    print(f"records={n} groups={len(groups)} in {time.time()-t0:.1f}s; wrote {out/'groups.json'}")


if __name__ == "__main__":
    main()
