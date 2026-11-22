#!/usr/bin/env python3
"""Find small-cell, non-magnetic MPtrj materials with meaningful forces."""
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
La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Tc Po At""".split())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--json", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-nat", type=int, default=8)
    ap.add_argument("--min-max-force", type=float, default=0.30)
    ap.add_argument("--max-magmom", type=float, default=0.5)
    ap.add_argument("--target", type=int, default=60)
    args = ap.parse_args()

    p = get_paths(args.config)
    json_path = Path(args.json) if args.json else p.dataset("mptraj")
    out = Path(args.out) if args.out else p.work / "mptraj" / "find_candidates.json"
    sssp = set(json.loads(p.pseudo_json("efficiency").read_text()))

    candidates = []
    t0 = time.time()
    n_materials = 0
    with json_path.open("rb") as f:
        for mp_id, traj_map in ijson.kvitems(f, ""):
            n_materials += 1
            if not isinstance(traj_map, dict) or not traj_map:
                continue
            max_force = 0.0
            max_mag = 0.0
            bandgap = None
            first_structure = None
            for rec in traj_map.values():
                if not isinstance(rec, dict):
                    continue
                if first_structure is None:
                    first_structure = rec.get("structure")
                force = rec.get("force", rec.get("forces"))
                if force is not None:
                    arr = np.asarray(force, dtype=float)
                    if arr.ndim == 2 and arr.shape[1] == 3 and arr.size:
                        max_force = max(max_force, float(np.max(np.linalg.norm(arr, axis=1))))
                mag = rec.get("magmom")
                if mag is not None:
                    arr = np.asarray(mag, dtype=float)
                    if arr.size:
                        max_mag = max(max_mag, float(np.max(np.abs(arr))))
                if bandgap is None and rec.get("bandgap") is not None:
                    bandgap = float(rec["bandgap"])
            if first_structure is None:
                continue
            els = [str(site["species"][0]["element"]) for site in first_structure.get("sites", [])]
            if not els or len(els) > args.max_nat:
                continue
            if not set(els).issubset(sssp) or set(els) & EXCLUDE:
                continue
            if max_force < args.min_max_force or max_mag > args.max_magmom:
                continue
            candidates.append({
                "mp_id": mp_id, "formula": first_structure.get("formula") or "".join(els),
                "nat": len(els), "elements": sorted(set(els)),
                "max_force_eV_A": max_force, "max_magmom": max_mag, "bandgap_eV": bandgap,
            })
            if len(candidates) >= args.target:
                break
    candidates.sort(key=lambda c: c["max_force_eV_A"], reverse=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(candidates, indent=2))
    print(f"scanned {n_materials} materials in {time.time()-t0:.1f}s; wrote {len(candidates)} to {out}")


if __name__ == "__main__":
    main()
