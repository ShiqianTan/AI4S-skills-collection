#!/usr/bin/env python3
"""Run MatPES-R2SCAN SCF jobs with the CPU libxc-enabled QE build."""
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


def run_job(job: dict, slot, force: bool, ctx: dict) -> dict:
    d = Path(job["dir"])
    if not d.is_absolute():
        d = (Path(ctx.get("base_dir", ".")) / d).resolve()
    out = d / "scf.out"
    res = dict(job)
    res["dir"] = str(d)
    res["slot"] = slot
    if not force and runner.is_done(out):
        res["exit"] = 0
    else:
        _ensure_pseudo_dir(d / "scf.in", ctx.get("pseudo_dir"))
        runner.prepare_tmp(d)
        run = runner.run_pw(d / "scf.in", out, pw=ctx["pw"], mpirun=ctx["mpirun"],
                            np=ctx["np"], env=ctx["env"], timeout=ctx["timeout"])
        res["exit"], res["wall_s"] = run["exit"], run["wall_s"]
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
    res.update({
        "nat": nat, "qe_total_energy_eV": qe_energy,
        "qe_energy_per_atom_eV": qe_energy / nat,
        "qe_forces": qe_forces.tolist(),
        "force_mae_eV_A": float(np.mean(np.abs(diff))),
        "force_rmse_eV_A": float(np.sqrt(np.mean(diff ** 2))),
        "force_max_err_eV_A": float(np.max(np.abs(diff))),
        "qe_max_force_eV_A": float(np.max(np.linalg.norm(qe_forces, axis=1))),
        "vasp_max_force_eV_A": float(np.max(np.linalg.norm(vasp_forces, axis=1))),
        "status": "ok",
    })
    try:
        qe_stress = np.asarray(atoms.get_stress(), dtype=float)
        res["qe_stress_eV_A3"] = qe_stress.tolist()
        if job.get("vasp_stress") is not None:
            qe_mat = -voigt_to_matrix(qe_stress)
            vasp_mat = voigt_to_matrix(np.asarray(job["vasp_stress"], dtype=float) * KBAR_TO_EV_A3)
            sdiff = qe_mat - vasp_mat
            res["stress_MAE_eV_A3"] = float(np.mean(np.abs(sdiff)))
            res["stress_max_err_eV_A3"] = float(np.max(np.abs(sdiff)))
    except Exception:  # noqa: BLE001
        res["qe_stress_eV_A3"] = None
    return res


def summarize(records: list[dict]) -> dict:
    ok = [r for r in records if r.get("status") == "ok"]
    stress = [r["stress_MAE_eV_A3"] for r in ok if "stress_MAE_eV_A3" in r]
    per_material = {}
    for mp_id in sorted({r["mp_id"] for r in ok}):
        group = sorted([r for r in ok if r["mp_id"] == mp_id],
                       key=lambda r: r["vasp_energy_per_atom_eV"])
        qe = np.array([r["qe_energy_per_atom_eV"] for r in group])
        va = np.array([r["vasp_energy_per_atom_eV"] for r in group])
        rel = (qe - qe.min()) - (va - va.min())
        per_material[mp_id] = {
            "formula": group[0]["formula"], "n_frames": len(group),
            "rel_energy_MAE_meV_atom": float(np.mean(np.abs(rel)) * 1000.0),
            "rel_energy_max_meV_atom": float(np.max(np.abs(rel)) * 1000.0),
        }
    return {
        "n_frames": len(ok), "n_materials": len(per_material),
        "force_MAE_eV_A": float(np.mean([r["force_mae_eV_A"] for r in ok])) if ok else None,
        "force_RMSE_eV_A": float(np.sqrt(np.mean([r["force_rmse_eV_A"] ** 2 for r in ok]))) if ok else None,
        "force_max_err_eV_A": float(np.max([r["force_max_err_eV_A"] for r in ok])) if ok else None,
        "stress_MAE_eV_A3": float(np.mean(stress)) if stress else None,
        "stress_max_err_eV_A3": float(np.max([r["stress_max_err_eV_A3"] for r in ok if "stress_max_err_eV_A3" in r])) if stress else None,
        "per_material": per_material,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--r2scan-dir", default=None)
    ap.add_argument("--tag", default="matpes_r2scan")
    ap.add_argument("--np", type=int, default=4, help="MPI ranks per job")
    ap.add_argument("--concurrent", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    p = get_paths(args.config)
    r2 = Path(args.r2scan_dir) if args.r2scan_dir else p.work / "matpes_r2scan"
    jobs = json.loads((r2 / "jobs.json").read_text())
    env = p.env()
    env["OMP_NUM_THREADS"] = str(p.cfg.get("qe", {}).get("omp_num_threads", 4))
    ctx = {"pw": p.pw("libxc"), "mpirun": p.mpirun(), "np": args.np,
           "env": env, "timeout": 6 * 3600,
           "pseudo_dir": p.pseudo_dir("pseudodojo"), "base_dir": r2}

    records = []
    print(f"running {len(jobs)} R2SCAN jobs ({args.np} MPI x {args.concurrent} concurrent)",
          flush=True)
    with ThreadPoolExecutor(max_workers=args.concurrent) as ex:
        futs = {ex.submit(run_job, job, i, args.force, ctx): job
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
    records.sort(key=lambda r: (r["mp_id"], r["traj_id"]))
    p.results.mkdir(parents=True, exist_ok=True)
    (p.results / f"{args.tag}_frames.json").write_text(json.dumps(records, indent=2))
    pd.DataFrame([{k: v for k, v in r.items()
                   if k not in {"vasp_forces", "qe_forces", "vasp_stress", "qe_stress_eV_A3"}}
                  for r in records]).to_csv(p.results / f"{args.tag}_frames.csv", index=False)
    summary = summarize(records)
    (p.results / f"{args.tag}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
