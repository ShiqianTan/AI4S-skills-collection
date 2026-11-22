#!/usr/bin/env python3
"""Generate QE static-SCF inputs for selected MPtrj reference frames."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core import mprelax, qeio  # noqa: E402
from qe_bench_core.structure import read_structure  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--mptraj-dir", default=None)
    ap.add_argument("--mode", choices=["sssp", "mprelax"], default="sssp")
    ap.add_argument("--spacing", type=float, default=0.20)
    ap.add_argument("--force-smearing", action="store_true")
    ap.add_argument("--pseudo-json", default=None)
    ap.add_argument("--pseudo-dir", default=None)
    args = ap.parse_args()

    p = get_paths(args.config)
    mptraj = Path(args.mptraj_dir) if args.mptraj_dir else p.work / (
        "mptraj_20" if (p.work / "mptraj_20" / "manifest.json").is_file() else "mptraj"
    )
    sssp_json = Path(args.pseudo_json) if args.pseudo_json else p.pseudo_json("efficiency")
    pseudo_dir = Path(args.pseudo_dir) if args.pseudo_dir else p.pseudo_dir("sssp_efficiency")
    sssp = json.loads(sssp_json.read_text())
    manifest = json.loads((mptraj / "manifest.json").read_text())

    jobs = []
    for mp_id, meta in manifest.items():
        gap = meta.get("bandgap_eV")
        for frame in meta["selected"]:
            sid = f"{mp_id}_{frame['traj_id']}"
            fdir = mptraj / mp_id / sid
            structure = read_structure(fdir / "structure.cif")
            if args.mode == "sssp":
                ecutwfc, ecutrho = qeio.cutoff_for(structure, sssp)
                kpts = qeio.kgrid(structure, spacing=args.spacing)
                k_shift = (0, 0, 0)
                occ = "smearing" if (args.force_smearing or gap is None or gap <= 0.2) else "fixed"
                settings = dict(occupations=occ)
            else:
                ecutwfc, ecutrho = mprelax.mp_cutoffs(structure, sssp)
                kpts, k_shift = mprelax.mp_kpoints(structure)
                conv_thr, degauss = mprelax.mp_electrons(len(structure))
                occ = "smearing"
                settings = dict(
                    occupations=occ, smearing="gaussian", degauss=degauss,
                    conv_thr=conv_thr, mixing_beta=mprelax.MP_SETTINGS["mixing_beta"],
                    electron_maxstep=mprelax.NELM, k_shift=k_shift,
                )
            qeio.write_qe_input(
                fdir / "scf.in", structure, f"mptraj_{sid}", "scf",
                ecutwfc, ecutrho, kpts, sssp, pseudo_dir=pseudo_dir, **settings,
            )
            jobs.append({
                "mp_id": mp_id, "formula": meta["formula"],
                "traj_id": frame["traj_id"], "sid": sid, "dir": str(fdir),
                "nat": frame["nat"], "ecutwfc": ecutwfc, "ecutrho": ecutrho,
                "kpts": list(kpts), "k_shift": list(k_shift),
                "occupations": occ, "config": args.mode,
                "mp_bandgap_eV": gap,
                "vasp_uncorrected_total_energy_eV": frame["vasp_uncorrected_total_energy_eV"],
                "vasp_energy_per_atom_eV": frame["vasp_energy_per_atom_eV"],
                "vasp_forces": frame["vasp_forces_eV_A"],
                "vasp_stress": frame["vasp_stress"],
                "vasp_max_force_eV_A": frame["vasp_max_force_eV_A"],
                "vasp_mean_force_eV_A": frame["vasp_mean_force_eV_A"],
                "volume_A3": frame["volume_A3"],
            })
    (mptraj / "jobs.json").write_text(json.dumps(jobs, indent=2))
    print(f"prepared {len(jobs)} MPtrj SCF jobs ({args.mode}) in {mptraj}")
    for job in jobs[:12]:
        print(job["sid"], job["formula"], "nat", job["nat"], "k", job["kpts"],
              job["occupations"], "Fmax(VASP)", round(job["vasp_max_force_eV_A"], 3))


if __name__ == "__main__":
    main()
