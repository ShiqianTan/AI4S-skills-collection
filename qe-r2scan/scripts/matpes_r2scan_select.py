#!/usr/bin/env python3
"""Extract the same MatPES frames from the R2SCAN dataset and prepare QE inputs.

R2SCAN requires libxc and norm-conserving pseudopotentials; this skill uses
PseudoDojo NC-SR PBE standard, 60/240 Ry, 0.25 A^-1 k-grid.
"""
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--json", default=None)
    ap.add_argument("--pbe-jobs", default=None, help="MatPES-PBE jobs.json with the same sids")
    ap.add_argument("--out", default=None)
    ap.add_argument("--ecutwfc", type=float, default=60.0)
    ap.add_argument("--ecutrho", type=float, default=240.0)
    ap.add_argument("--spacing", type=float, default=0.25)
    args = ap.parse_args()

    from pymatgen.core import Structure

    p = get_paths(args.config)
    json_path = Path(args.json) if args.json else p.dataset("matpes") / "MatPES-R2SCAN-2025.1.json"
    out = Path(args.out) if args.out else p.work / "matpes_r2scan"
    out.mkdir(parents=True, exist_ok=True)
    pbe_jobs_path = Path(args.pbe_jobs) if args.pbe_jobs else p.work / "matpes" / "jobs.json"
    pbe_jobs = json.loads(pbe_jobs_path.read_text())
    selected = {j["sid"]: j for j in pbe_jobs}
    pseudo_dir = p.pseudo_dir("pseudodojo")

    found = {}
    with json_path.open("rb") as f:
        for rec in ijson.items(f, "item", use_float=True):
            mid = str(rec.get("matpes_id", ""))
            if mid in selected:
                found[mid] = rec
                if len(found) == len(selected):
                    break
    missing = sorted(set(selected) - set(found))
    if missing:
        print(f"warning: {len(missing)} records missing from the R2SCAN dataset")

    jobs = []
    for mid, pbe in selected.items():
        rec = found.get(mid)
        if rec is None:
            continue
        structure = Structure.from_dict(rec["structure"])
        forces = np.asarray(rec["forces"], dtype=float)
        fmax = float(np.max(np.linalg.norm(forces, axis=1)))
        elements = [str(e) for e in structure.composition.elements]
        pmap = {el: {"filename": f"{el}.upf"} for el in elements}
        gk = pbe["mp_id"]
        formula = pbe["formula"]
        d = out / gk / mid
        d.mkdir(parents=True, exist_ok=True)
        structure.to(filename=str(d / "structure.cif"))
        (d / "r2scan_frame.json").write_text(json.dumps({
            "group": gk, "formula": formula, "matpes_id": mid, "nat": len(structure),
            "matpes_r2scan_total_energy_eV": float(rec["energy"]),
            "matpes_r2scan_energy_per_atom_eV": float(rec["energy"]) / len(structure),
            "matpes_r2scan_forces_eV_A": forces.tolist(),
            "matpes_r2scan_stress": rec.get("stress"),
            "matpes_r2scan_max_force_eV_A": fmax,
            "matpes_r2scan_mean_force_eV_A": float(np.mean(np.linalg.norm(forces, axis=1))),
            "matpes_r2scan_bandgap_eV": rec.get("bandgap"),
            "volume_A3": float(structure.volume),
        }, indent=2))
        kpts = qeio.kgrid(structure, spacing=args.spacing)
        qeio.write_qe_input(d / "scf.in", structure, f"matpes_r2scan_{mid}", "scf",
                            args.ecutwfc, args.ecutrho, kpts, pmap,
                            pseudo_dir=pseudo_dir, occupations="smearing",
                            input_dft="r2scan")
        jobs.append({
            "mp_id": gk, "formula": formula, "traj_id": mid, "sid": mid,
            "dir": str(d.relative_to(out)), "nat": len(structure), "ecutwfc": args.ecutwfc,
            "ecutrho": args.ecutrho, "kpts": list(kpts), "occupations": "smearing",
            "functional": "r2scan", "mp_bandgap_eV": rec.get("bandgap"),
            "stress_unit": "kbar", "stress_sign": -1,
            "vasp_uncorrected_total_energy_eV": float(rec["energy"]),
            "vasp_energy_per_atom_eV": float(rec["energy"]) / len(structure),
            "vasp_forces": forces.tolist(), "vasp_stress": rec.get("stress"),
            "vasp_max_force_eV_A": fmax,
            "vasp_mean_force_eV_A": float(np.mean(np.linalg.norm(forces, axis=1))),
            "volume_A3": float(structure.volume),
        })
    (out / "jobs.json").write_text(json.dumps(jobs, indent=2))
    print(f"prepared {len(jobs)} R2SCAN jobs in {out}; missing={len(missing)}")


if __name__ == "__main__":
    main()
