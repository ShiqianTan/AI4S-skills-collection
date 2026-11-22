#!/usr/bin/env python3
"""Select sAlex frames, save VASP references and write QE SCF inputs."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from ase.db import connect  # noqa: E402
from ase.io import write  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core import qeio  # noqa: E402
from qe_bench_core.structure import structure_from_ase  # noqa: E402

# 20 materials x 4 frames selected in the study (see references/salex-results.md).
SELECTED = {
    "agm003604844": ["agm003604844_0_0", "agm003604844_1_27", "agm003604844_2_58", "agm003604844_2_68"],
    "agm003610536": ["agm003610536_2_0", "agm003610536_2_105", "agm003610536_2_133", "agm003610536_2_127"],
    "agm003602335": ["agm003602335_0_0", "agm003602335_2_137", "agm003602335_2_154", "agm003602335_2_162"],
    "agm003598187": ["agm003598187_0_0", "agm003598187_2_202", "agm003598187_2_215", "agm003598187_2_233"],
    "agm003611662": ["agm003611662_0_0", "agm003611662_2_85", "agm003611662_2_124", "agm003611662_2_123"],
    "agm003607271": ["agm003607271_0_0", "agm003607271_1_22", "agm003607271_1_27", "agm003607271_1_31"],
    "agm003610555": ["agm003610555_0_0", "agm003610555_2_49", "agm003610555_2_55", "agm003610555_2_58"],
    "agm003601629": ["agm003601629_0_0", "agm003601629_2_103", "agm003601629_2_109", "agm003601629_2_102"],
    "agm003598814": ["agm003598814_0_0", "agm003598814_2_51", "agm003598814_2_64", "agm003598814_2_75"],
    "agm003600957": ["agm003600957_0_0", "agm003600957_2_87", "agm003600957_2_91", "agm003600957_2_99"],
    "agm003550318": ["agm003550318_0_0", "agm003550318_2_63", "agm003550318_2_46", "agm003550318_2_42"],
    "agm003604297": ["agm003604297_2_0", "agm003604297_2_23", "agm003604297_2_26", "agm003604297_2_28"],
    "agm003600888": ["agm003600888_0_0", "agm003600888_2_32", "agm003600888_2_37", "agm003600888_2_86"],
    "agm003468171": ["agm003468171_0_0", "agm003468171_2_149", "agm003468171_2_157", "agm003468171_2_238"],
    "agm003468142": ["agm003468142_0_0", "agm003468142_2_73", "agm003468142_2_75", "agm003468142_2_77"],
    "agm003469447": ["agm003469447_0_0", "agm003469447_2_90", "agm003469447_2_202", "agm003469447_2_208"],
    "agm003468655": ["agm003468655_0_0", "agm003468655_2_127", "agm003468655_2_131", "agm003468655_2_187"],
    "agm003470612": ["agm003470612_0_0", "agm003470612_2_235", "agm003470612_2_108", "agm003470612_2_224"],
    "agm003470115": ["agm003470115_2_0", "agm003470115_2_161", "agm003470115_2_227", "agm003470115_2_222"],
    "agm003604011": ["agm003604011_2_0", "agm003604011_2_155", "agm003604011_2_174", "agm003604011_2_164"],
}
FORMULAS = {
    "agm003604844": "Ag2Hg2S2", "agm003610536": "Cu2Ge2P2", "agm003602335": "As2O2S2",
    "agm003598187": "Ge2I2Se2", "agm003611662": "F2Ge2P2", "agm003607271": "Ca2Ge2Sn2",
    "agm003610555": "As2Br2Te2", "agm003601629": "Br2S2Sb2", "agm003598814": "Ag2As2Te2",
    "agm003600957": "O2S2Sb2", "agm003550318": "Ag4Mg3Sn", "agm003604297": "H2Pd2Sr2",
    "agm003600888": "Bi2Rb2Sb2", "agm003468171": "F5GeHg2", "agm003468142": "Br5GeIn2",
    "agm003469447": "F5InK2", "agm003468655": "Bi2InO5", "agm003470612": "K2PS5",
    "agm003470115": "BNa2Se5", "agm003604011": "H2O2Zn2",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--val", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--spacing", type=float, default=0.20)
    args = ap.parse_args()

    p = get_paths(args.config)
    val = Path(args.val) if args.val else p.dataset("salex")
    out = Path(args.out) if args.out else p.work / "salex"
    out.mkdir(parents=True, exist_ok=True)
    sssp = json.loads(p.pseudo_json("efficiency").read_text())
    pseudo_dir = p.pseudo_dir("sssp_efficiency")

    groups = {g["group"]: g for g in json.loads((out / "groups.json").read_text())}
    all_sids = {sid for sids in SELECTED.values() for sid in sids}
    found = {}
    for group, sids in SELECTED.items():
        if group not in groups:
            continue
        for fname in groups[group]["files"]:
            db = connect(str(val / fname))
            for row in db.select():
                sid = str(row.data.get("sid", ""))
                if sid not in all_sids or not sid.startswith(group + "_"):
                    continue
                atoms = row.toatoms()
                forces = np.asarray(atoms.calc.results["forces"], dtype=float)
                stress = atoms.calc.results.get("stress")
                found[sid] = {
                    "group": group, "formula": FORMULAS[group], "sid": sid, "atoms": atoms,
                    "energy": float(row.energy), "forces": forces,
                    "stress": np.asarray(stress, dtype=float).tolist() if stress is not None else None,
                    "fmax": float(np.max(np.linalg.norm(forces, axis=1))),
                    "mean_force": float(np.mean(np.linalg.norm(forces, axis=1))),
                }
    missing = sorted(all_sids - set(found))
    if missing:
        raise RuntimeError(f"missing sids: {missing}")

    jobs = []
    for group, sids in SELECTED.items():
        for sid in sids:
            item = found[sid]
            atoms = item["atoms"]
            structure = structure_from_ase(atoms)
            d = out / group / sid
            d.mkdir(parents=True, exist_ok=True)
            write(d / "structure.cif", atoms)
            (d / "salex_frame.json").write_text(json.dumps({
                "group": group, "formula": item["formula"], "sid": sid, "nat": len(atoms),
                "salex_total_energy_eV": item["energy"],
                "salex_energy_per_atom_eV": item["energy"] / len(atoms),
                "salex_forces_eV_A": item["forces"].tolist(),
                "salex_stress": item["stress"],
                "salex_max_force_eV_A": item["fmax"],
                "salex_mean_force_eV_A": item["mean_force"],
                "volume_A3": float(atoms.get_volume()),
            }, indent=2))
            ecutwfc, ecutrho = qeio.cutoff_for(structure, sssp)
            kpts = qeio.kgrid(structure, spacing=args.spacing)
            qeio.write_qe_input(d / "scf.in", structure, f"salex_{sid}", "scf",
                                ecutwfc, ecutrho, kpts, sssp,
                                pseudo_dir=pseudo_dir, occupations="smearing")
            jobs.append({
                "mp_id": group, "formula": item["formula"], "traj_id": sid, "sid": sid,
                "dir": str(d), "nat": len(atoms), "ecutwfc": ecutwfc, "ecutrho": ecutrho,
                "kpts": list(kpts), "occupations": "smearing", "mp_bandgap_eV": None,
                "stress_unit": "eV_A3", "stress_sign": 1,
                "vasp_uncorrected_total_energy_eV": item["energy"],
                "vasp_energy_per_atom_eV": item["energy"] / len(atoms),
                "vasp_forces": item["forces"].tolist(), "vasp_stress": item["stress"],
                "vasp_max_force_eV_A": item["fmax"],
                "vasp_mean_force_eV_A": item["mean_force"],
                "volume_A3": float(atoms.get_volume()),
            })
    (out / "jobs.json").write_text(json.dumps(jobs, indent=2))
    print(f"prepared {len(jobs)} sAlex frames for {len(SELECTED)} materials in {out}")


if __name__ == "__main__":
    main()
