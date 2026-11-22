#!/usr/bin/env python3
"""Select MatPES-PBE frames (10 materials x 4 frames) and prepare QE inputs."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ijson  # noqa: E402
import numpy as np  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core import qeio  # noqa: E402
from qe_bench_core.structure import read_structure  # noqa: E402

GROUPS = {
    "matpes-20240214_20812": "InN", "matpes-20240214_2605": "CaO",
    "matpes-20240214_7592": "AgF", "matpes-20240214_10759": "GeSe",
    "matpes-20240214_1096374": "TlCuPd2", "matpes-20240214_19717": "TePb",
    "matpes-20240214_1001581": "LiC6", "matpes-20240214_1018088": "C",
    "matpes-20240214_1184938": "KH3", "matpes-20240214_1022724": "CaC4",
}
TARGETS = [0.01, 0.6, 1.5, 2.5]
MAX_FORCE = 3.0


def group_key(mid: str) -> str:
    return mid.rsplit("_", 1)[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--json", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--spacing", type=float, default=0.20)
    args = ap.parse_args()

    from pymatgen.core import Structure

    p = get_paths(args.config)
    json_path = Path(args.json) if args.json else p.dataset("matpes") / "MatPES-PBE-2025.1.json"
    out = Path(args.out) if args.out else p.work / "matpes"
    out.mkdir(parents=True, exist_ok=True)
    sssp = json.loads(p.pseudo_json("efficiency").read_text())
    pseudo_dir = p.pseudo_dir("sssp_efficiency")

    frames_by_group: dict[str, list[dict]] = {g: [] for g in GROUPS}
    with json_path.open("rb") as f:
        for rec in ijson.items(f, "item", use_float=True):
            mid = str(rec.get("matpes_id", ""))
            gk = group_key(mid)
            if gk not in GROUPS:
                continue
            structure_dict = rec.get("structure")
            if not isinstance(structure_dict, dict):
                continue
            structure = Structure.from_dict(structure_dict)
            forces = np.asarray(rec.get("forces", []), dtype=float)
            if forces.ndim != 2 or forces.shape[1] != 3:
                continue
            fmax = float(np.max(np.linalg.norm(forces, axis=1)))
            frames_by_group[gk].append({
                "matpes_id": mid, "structure": structure,
                "energy": float(rec.get("energy")), "forces": forces,
                "stress": rec.get("stress"), "fmax": fmax,
                "mean_force": float(np.mean(np.linalg.norm(forces, axis=1))),
                "bandgap": rec.get("bandgap"), "nat": len(structure),
            })

    jobs = []
    for gk, formula in GROUPS.items():
        frames = frames_by_group[gk]
        if not frames:
            raise RuntimeError(f"no frames for {gk} {formula}")
        emin = min(fr["energy"] for fr in frames)
        chosen, used = [], set()
        for target in TARGETS:
            cand = []
            for window in (5.0, 10.0, 1e9):
                cand = [fr for fr in frames if fr["matpes_id"] not in used
                        and fr["fmax"] <= MAX_FORCE and fr["energy"] - emin <= window]
                if cand:
                    break
            if not cand:
                cand = [fr for fr in frames if fr["matpes_id"] not in used]
            best = min(cand, key=lambda fr: abs(fr["fmax"] - target))
            chosen.append(best)
            used.add(best["matpes_id"])

        for fr in chosen:
            sid = fr["matpes_id"]
            d = out / gk / sid
            d.mkdir(parents=True, exist_ok=True)
            fr["structure"].to(filename=str(d / "structure.cif"))
            (d / "matpes_frame.json").write_text(json.dumps({
                "group": gk, "formula": formula, "matpes_id": sid, "nat": fr["nat"],
                "matpes_total_energy_eV": fr["energy"],
                "matpes_energy_per_atom_eV": fr["energy"] / fr["nat"],
                "matpes_forces_eV_A": fr["forces"].tolist(),
                "matpes_stress": fr["stress"],
                "matpes_max_force_eV_A": fr["fmax"],
                "matpes_mean_force_eV_A": fr["mean_force"],
                "matpes_bandgap_eV": fr["bandgap"],
                "volume_A3": float(fr["structure"].volume),
            }, indent=2))
            ecutwfc, ecutrho = qeio.cutoff_for(fr["structure"], sssp)
            kpts = qeio.kgrid(fr["structure"], spacing=args.spacing)
            qeio.write_qe_input(d / "scf.in", fr["structure"], f"matpes_{sid}", "scf",
                                ecutwfc, ecutrho, kpts, sssp,
                                pseudo_dir=pseudo_dir, occupations="smearing")
            jobs.append({
                "mp_id": gk, "formula": formula, "traj_id": sid, "sid": sid,
                "dir": str(d.relative_to(out)), "nat": fr["nat"], "ecutwfc": ecutwfc, "ecutrho": ecutrho,
                "kpts": list(kpts), "occupations": "smearing", "mp_bandgap_eV": fr["bandgap"],
                "stress_unit": "kbar", "stress_sign": -1,
                "vasp_uncorrected_total_energy_eV": fr["energy"],
                "vasp_energy_per_atom_eV": fr["energy"] / fr["nat"],
                "vasp_forces": fr["forces"].tolist(), "vasp_stress": fr["stress"],
                "vasp_max_force_eV_A": fr["fmax"],
                "vasp_mean_force_eV_A": fr["mean_force"],
                "volume_A3": float(fr["structure"].volume),
            })
    (out / "jobs.json").write_text(json.dumps(jobs, indent=2))
    print(f"prepared {len(jobs)} MatPES frames for {len(GROUPS)} materials in {out}")


if __name__ == "__main__":
    main()
