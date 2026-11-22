#!/usr/bin/env python3
"""Run QE static SCF for selected MPtrj frames and compare forces/stress."""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core import runner  # noqa: E402
from qe_bench_core.parse import KBAR_TO_EV_A3, voigt_to_matrix  # noqa: E402



def _ensure_pseudo_dir(inp: Path, pseudo_dir) -> None:
    import re

    text = inp.read_text()
    new = re.sub(r"pseudo_dir\s*=\s*'[^']*'", f"pseudo_dir = '{pseudo_dir}'", text)
    if new != text:
        inp.write_text(new)


def run_scf(job: dict, gpu, force: bool, ctx: dict) -> dict:
    d = Path(job["dir"])
    if not d.is_absolute():
        d = (Path(ctx.get("base_dir", ".")) / d).resolve()
    out = d / "scf.out"
    res = dict(job)
    res["dir"] = str(d)
    res["gpu"] = gpu
    if not force and runner.is_done(out):
        res["scf_exit"] = 0
        res["scf_wall_s"] = None
    else:
        _ensure_pseudo_dir(d / "scf.in", ctx.get("pseudo_dir"))
        runner.prepare_tmp(d)
        run = runner.run_pw(d / "scf.in", out,
                            pw=ctx["pw"], mpirun=ctx["mpirun"], np=ctx["np"],
                            env=ctx["env"], gpu=gpu, timeout=ctx["timeout"])
        res["scf_exit"] = run["exit"]
        res["scf_wall_s"] = run["wall_s"]
        if run["exit"] != 0 or not run["done"]:
            res["status"] = "scf_failed"
            return res

    from ase.io import read

    atoms = read(str(out), format="espresso-out")
    qe_forces = np.asarray(atoms.get_forces(), dtype=float)
    vasp_forces = np.asarray(job["vasp_forces"], dtype=float)
    if qe_forces.shape != vasp_forces.shape:
        res["status"] = f"shape_mismatch {qe_forces.shape} vs {vasp_forces.shape}"
        return res

    diff = qe_forces - vasp_forces
    nat = len(atoms)
    qe_energy = float(atoms.get_potential_energy())

    qe_gap = None
    if job.get("occupations", "fixed") == "fixed":
        for line in out.read_text(errors="ignore").splitlines():
            if "highest occupied, lowest unoccupied level" in line:
                vals = line.split(":")[1].split()
                qe_gap = float(vals[1]) - float(vals[0])
            elif "highest occupied level" in line and "lowest unoccupied" not in line:
                qe_gap = 0.0
    if qe_gap is None and job.get("mp_bandgap_eV") == 0:
        qe_gap = 0.0

    res.update({
        "nat": nat,
        "qe_total_energy_eV": qe_energy,
        "qe_energy_per_atom_eV": qe_energy / nat,
        "qe_forces": qe_forces.tolist(),
        "force_mae_eV_A": float(np.mean(np.abs(diff))),
        "force_rmse_eV_A": float(np.sqrt(np.mean(diff ** 2))),
        "force_max_err_eV_A": float(np.max(np.abs(diff))),
        "force_mae_norm_eV_A": float(np.mean(np.abs(
            np.linalg.norm(qe_forces, axis=1) - np.linalg.norm(vasp_forces, axis=1)))),
        "qe_max_force_eV_A": float(np.max(np.linalg.norm(qe_forces, axis=1))),
        "vasp_max_force_eV_A": float(np.max(np.linalg.norm(vasp_forces, axis=1))),
        "qe_volume_A3": float(atoms.get_volume()),
        "vasp_volume_A3": job["volume_A3"],
        "qe_band_gap_eV": qe_gap,
        "mp_bandgap_eV": job.get("mp_bandgap_eV"),
        "status": "ok",
    })
    if qe_gap is not None and job.get("mp_bandgap_eV") is not None:
        res["gap_diff_eV"] = qe_gap - job["mp_bandgap_eV"]
    try:
        qe_stress = np.asarray(atoms.get_stress(), dtype=float)
        res["qe_stress_eV_A3"] = qe_stress.tolist()
        if job.get("vasp_stress") is not None:
            qe_mat = voigt_to_matrix(qe_stress)
            if job.get("stress_sign", -1) == -1:
                # QE stress is opposite to the VASP/MPtrj convention.
                qe_mat = -qe_mat
            vasp_stress = np.asarray(job["vasp_stress"], dtype=float)
            if job.get("stress_unit", "kbar") == "eV_A3":
                vasp_mat = voigt_to_matrix(vasp_stress) if vasp_stress.shape == (6,) else vasp_stress
            else:
                vasp_mat = vasp_stress * KBAR_TO_EV_A3
                if vasp_mat.shape == (6,):
                    vasp_mat = voigt_to_matrix(vasp_mat)
            sdiff = qe_mat - vasp_mat
            res["stress_MAE_eV_A3"] = float(np.mean(np.abs(sdiff)))
            res["stress_max_err_eV_A3"] = float(np.max(np.abs(sdiff)))
    except Exception:  # noqa: BLE001
        res["qe_stress_eV_A3"] = None
    return res


def summarize(records: list[dict]) -> dict:
    ok = [r for r in records if r.get("status") == "ok"]
    if not ok:
        return {"n_frames": 0}
    per_material = {}
    for mp_id in sorted({r["mp_id"] for r in ok}):
        group = sorted([r for r in ok if r["mp_id"] == mp_id],
                       key=lambda r: r["vasp_energy_per_atom_eV"])
        qe_e = np.array([r["qe_energy_per_atom_eV"] for r in group])
        vasp_e = np.array([r["vasp_energy_per_atom_eV"] for r in group])
        rel_err = (qe_e - qe_e.min()) - (vasp_e - vasp_e.min())
        stress_vals = [r["stress_MAE_eV_A3"] for r in group if "stress_MAE_eV_A3" in r]
        per_material[mp_id] = {
            "formula": group[0]["formula"], "n_frames": len(group),
            "energy_offset_eV_atom": float(np.mean(qe_e - vasp_e)),
            "rel_energy_MAE_meV_atom": float(np.mean(np.abs(rel_err)) * 1000.0),
            "rel_energy_max_meV_atom": float(np.max(np.abs(rel_err)) * 1000.0),
            "force_MAE_eV_A": float(np.mean([r["force_mae_eV_A"] for r in group])),
            "force_RMSE_eV_A": float(np.sqrt(np.mean([r["force_rmse_eV_A"] ** 2 for r in group]))),
            "force_max_err_eV_A": float(np.max([r["force_max_err_eV_A"] for r in group])),
            "stress_MAE_eV_A3": float(np.mean(stress_vals)) if stress_vals else None,
            "mp_bandgap_eV": group[0].get("mp_bandgap_eV"),
        }
    stress_vals = [r["stress_MAE_eV_A3"] for r in ok if "stress_MAE_eV_A3" in r]
    stress_max_vals = [r["stress_max_err_eV_A3"] for r in ok if "stress_max_err_eV_A3" in r]
    gap_recs = [r for r in ok if r.get("mp_bandgap_eV") and r.get("qe_band_gap_eV") is not None]
    return {
        "n_frames": len(ok), "n_materials": len(per_material),
        "force_MAE_eV_A": float(np.mean([r["force_mae_eV_A"] for r in ok])),
        "force_RMSE_eV_A": float(np.sqrt(np.mean([r["force_rmse_eV_A"] ** 2 for r in ok]))),
        "force_max_err_eV_A": float(np.max([r["force_max_err_eV_A"] for r in ok])),
        "stress_MAE_eV_A3": float(np.mean(stress_vals)) if stress_vals else None,
        "stress_max_err_eV_A3": float(np.max(stress_max_vals)) if stress_max_vals else None,
        "gap_MAE_eV": float(np.mean([abs(r["gap_diff_eV"]) for r in gap_recs])) if gap_recs else None,
        "n_gap_materials": len(gap_recs),
        "per_material": per_material,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--mptraj-dir", default=None)
    ap.add_argument("--tag", default="mptraj20")
    ap.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    ap.add_argument("--gpu-mode", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    p = get_paths(args.config)
    mptraj = Path(args.mptraj_dir) if args.mptraj_dir else p.work / "mptraj_20"
    kind = "gpu" if args.gpu_mode == "gpu" else "cpu"
    ctx = {"pw": p.pw(kind), "mpirun": p.mpirun(), "np": args.np,
           "env": p.env(), "timeout": 6 * 3600,
           "pseudo_dir": p.pseudo_dir("sssp_efficiency"), "base_dir": mptraj}
    jobs = json.loads((mptraj / "jobs.json").read_text())
    devices = [int(x) for x in args.devices.split(",")]
    gpu_values = devices if args.gpu_mode == "gpu" else [None]

    records = []
    print(f"running {len(jobs)} MPtrj SCF jobs ({args.gpu_mode})", flush=True)
    with ThreadPoolExecutor(max_workers=len(gpu_values)) as ex:
        futs = {ex.submit(run_scf, job, gpu_values[i % len(gpu_values)], args.force, ctx): job
                for i, job in enumerate(jobs)}
        for fut in as_completed(futs):
            job = futs[fut]
            try:
                res = fut.result()
            except Exception as exc:  # noqa: BLE001
                res = dict(job)
                res["status"] = f"exception: {exc}"
            records.append(res)
            print(job["sid"], res.get("status"), res.get("force_rmse_eV_A"), flush=True)

    records = sorted(records, key=lambda r: (r["mp_id"], r["traj_id"]))
    p.results.mkdir(parents=True, exist_ok=True)
    (p.results / f"{args.tag}_frames.json").write_text(json.dumps(records, indent=2))
    df = pd.DataFrame([
        {k: v for k, v in r.items()
         if k not in {"vasp_forces", "qe_forces", "vasp_stress", "qe_stress_eV_A3"}}
        for r in records
    ])
    df.to_csv(p.results / f"{args.tag}_frames.csv", index=False)
    summary = summarize(records)
    (p.results / f"{args.tag}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
