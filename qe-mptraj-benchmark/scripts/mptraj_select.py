#!/usr/bin/env python3
"""Stream MPtrj, select representative frames, and save VASP reference data.

For each material it keeps the lowest-energy, highest-energy and largest-force
frames (plus a median-energy frame if fewer than three exist).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from decimal import Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ijson  # noqa: E402
import numpy as np  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402

TARGETS_5 = {
    "mp-2172": "AlAs", "mp-2176": "ZnTe", "mp-22925": "AgI",
    "mp-2624": "AlSb", "mp-1039": "MgTe",
    "mp-12779": "CdTe", "mp-22870": "InBr", "mp-1747": "K2Te",
    "mp-22903": "RbI", "mp-22867": "RbBr",
}
TARGETS_20 = {
    "mp-1093708": "TiMoPd", "mp-1025340": "CuWSe4", "mp-1068813": "KTe2Pd",
    "mp-1147776": "KVSe4", "mp-1097229": "NaIn2Bi", "mp-1079080": "CsCd4As3",
    "mp-1097108": "K2HgTe", "mp-1096820": "AlAgS2", "mp-1079269": "SiPt6",
    "mp-1096742": "YGaPt2", "mp-1096568": "LiCu2Au", "mp-1099282": "Mg6CuNi",
    "mp-1099109": "CaMg6Ni", "mp-1093539": "Ca2AlZn", "mp-10264": "Ni2Te2",
    "mp-1094623": "Mg4Ga2", "mp-1147775": "Ba2OsO4", "mp-11514": "NbPt2",
    "mp-1023918": "Al2Sb2", "mp-1014231": "Hf2Zn",
}


def fnum(value):
    if value is None:
        return None
    try:
        return float(value)
    except Exception:  # noqa: BLE001
        return None


def as_matrix(value):
    if value is None:
        return None
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 2 and arr.shape[1] == 3:
        return arr.tolist()
    if arr.shape == (3, 3) or arr.size == 9:
        return arr.reshape(3, 3).tolist()
    return None


def json_default(obj):
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Decimal):
        return float(obj)
    raise TypeError(type(obj))


def process_material(mp_id, formula, traj_map, out_root):
    from pymatgen.core import Structure

    frames = []
    bandgap = None
    for traj_id, rec in traj_map.items():
        if not isinstance(rec, dict):
            continue
        try:
            structure = Structure.from_dict(rec["structure"])
        except Exception:  # noqa: BLE001
            continue
        raw_energy = fnum(rec.get("uncorrected_total_energy", rec.get("corrected_total_energy")))
        corr_energy = fnum(rec.get("corrected_total_energy"))
        energy_per_atom = fnum(rec.get("energy_per_atom"))
        if energy_per_atom is None and raw_energy is not None:
            energy_per_atom = raw_energy / len(structure)
        forces = as_matrix(rec.get("force", rec.get("forces")))
        stress = as_matrix(rec.get("stress"))
        if bandgap is None and rec.get("bandgap") is not None:
            bandgap = fnum(rec.get("bandgap"))
        force_norms = np.linalg.norm(np.asarray(forces), axis=1) if forces else np.array([])
        frames.append({
            "traj_id": str(traj_id), "structure": structure, "nat": len(structure),
            "raw_energy": raw_energy, "corrected_energy": corr_energy,
            "energy_per_atom": energy_per_atom, "forces": forces, "stress": stress,
            "max_force": float(np.max(force_norms)) if force_norms.size else None,
            "mean_force": float(np.mean(force_norms)) if force_norms.size else None,
            "volume": float(structure.volume),
        })
    if not frames:
        return None

    energies = np.array([f["energy_per_atom"] for f in frames if f["energy_per_atom"] is not None])
    order = np.argsort(energies)
    candidates = [int(order[0]), int(order[-1]),
                  int(np.argmax([f["max_force"] if f["max_force"] is not None else -1 for f in frames]))]
    candidates.append(int(order[len(order) // 2]))
    selected = []
    for idx in candidates:
        if idx not in selected:
            selected.append(idx)
        if len(selected) == 3:
            break

    md = out_root / mp_id
    md.mkdir(parents=True, exist_ok=True)
    saved = []
    for rank, idx in enumerate(selected):
        f = frames[idx]
        sid = f"{mp_id}_{f['traj_id']}".replace("/", "_")
        fdir = md / sid
        fdir.mkdir(parents=True, exist_ok=True)
        f["structure"].to(filename=str(fdir / "structure.cif"))
        payload = {
            "mp_id": mp_id, "formula": formula, "traj_id": f["traj_id"], "nat": f["nat"],
            "vasp_uncorrected_total_energy_eV": f["raw_energy"],
            "vasp_corrected_total_energy_eV": f["corrected_energy"],
            "vasp_energy_per_atom_eV": f["energy_per_atom"],
            "vasp_forces_eV_A": f["forces"], "vasp_stress": f["stress"],
            "vasp_max_force_eV_A": f["max_force"],
            "vasp_mean_force_eV_A": f["mean_force"],
            "volume_A3": f["volume"], "selection_rank": rank,
        }
        (fdir / "vasp_frame.json").write_text(json.dumps(payload, indent=2, default=json_default))
        saved.append(payload)

    (md / "frames.json").write_text(json.dumps([
        {"traj_id": f["traj_id"], "nat": f["nat"], "energy_per_atom_eV": f["energy_per_atom"],
         "max_force_eV_A": f["max_force"], "mean_force_eV_A": f["mean_force"],
         "volume_A3": f["volume"]}
        for f in frames
    ], indent=2, default=json_default))
    print(f"{mp_id} {formula}: {len(frames)} frames, nat={frames[0]['nat']}, selected {len(saved)}", flush=True)
    return {
        "mp_id": mp_id, "formula": formula, "n_frames": len(frames),
        "nat": frames[0]["nat"],
        "energy_per_atom_min": float(energies.min()),
        "energy_per_atom_max": float(energies.max()),
        "max_force_min": float(min(f["max_force"] for f in frames if f["max_force"] is not None)),
        "bandgap_eV": bandgap, "selected": saved,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--json", default=None, help="MPtrj_2022.9_full.json")
    ap.add_argument("--preset", choices=["5", "20"], default="20")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    p = get_paths(args.config)
    json_path = Path(args.json) if args.json else p.dataset("mptraj")
    out_root = Path(args.out) if args.out else p.work / ("mptraj_20" if args.preset == "20" else "mptraj")
    out_root.mkdir(parents=True, exist_ok=True)
    targets = TARGETS_5 if args.preset == "5" else TARGETS_20

    found = {}
    t0 = time.time()
    with json_path.open("rb") as f:
        for mp_id, traj_map in ijson.kvitems(f, ""):
            if mp_id not in targets or not isinstance(traj_map, dict):
                continue
            result = process_material(mp_id, targets[mp_id], traj_map, out_root)
            if result is not None:
                found[mp_id] = result
            if len(found) == len(targets):
                break
    (out_root / "manifest.json").write_text(json.dumps(found, indent=2, default=json_default))
    missing = sorted(set(targets) - set(found))
    print(f"done in {time.time() - t0:.1f}s; found={len(found)}; missing={missing}", flush=True)


if __name__ == "__main__":
    main()
