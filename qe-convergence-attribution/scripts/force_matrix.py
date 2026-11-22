#!/usr/bin/env python3
"""High-force MPtrj convergence matrix (one factor at a time).

Selects the worst-force frames from a prior MPtrj run and re-evaluates them with
variations of k-point density, plane-wave cutoff, smearing and pseudopotential
set. Forces are compared with the MPtrj/VASP reference.
"""
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
from qe_bench_core import parse, pseudos, qeio, runner  # noqa: E402
from qe_bench_core.structure import read_structure  # noqa: E402

RY_TO_EV = 13.605693122994
KBAR_TO_EV_A3 = 0.0006241509
DEGAUSS = 0.05 / RY_TO_EV


def _pseudo(p, name: str, elements: list[str]):
    if name == "sssp_eff":
        upf = p.pseudo_dir("sssp_efficiency")
        return pseudos.sssp_map(json.loads(p.pseudo_json("efficiency").read_text()), elements, upf), upf
    if name == "sssp_prec":
        upf = p.pseudo_dir("sssp_precision")
        return pseudos.sssp_map(json.loads(p.pseudo_json("precision").read_text()), elements, upf), upf
    if name == "pdojo":
        upf = p.pseudo_dir("pseudodojo")
        return pseudos.pseudodojo_map(elements, upf), upf
    raise ValueError(name)


def select_frames(args, p, mptraj: Path) -> list[dict]:
    jobs = {j["sid"]: j for j in json.loads((mptraj / "jobs.json").read_text())}
    if args.sids:
        sids = [x for x in args.sids.split(",") if x]
    else:
        csv = Path(args.frames_csv) if args.frames_csv else p.results / "mptraj20_frames.csv"
        if not csv.is_file():
            raise SystemExit(f"need --sids or an existing --frames-csv ({csv})")
        df = pd.read_csv(csv)
        df = df[df.get("status", "ok") == "ok"] if "status" in df else df
        sids = list(df.sort_values("force_max_err_eV_A", ascending=False)
                    .drop_duplicates("mp_id").head(args.top_n)["sid"])
    out = []
    for sid in sids:
        job = jobs[sid]
        d = Path(job["dir"])
        frame = json.loads((d / "vasp_frame.json").read_text())
        out.append({
            "mp_id": job["mp_id"], "formula": job["formula"], "traj_id": job["traj_id"],
            "sid": sid, "dir": str(d), "structure": read_structure(d / "structure.cif"),
            "vasp_forces": np.asarray(frame["vasp_forces_eV_A"], dtype=float),
            "vasp_stress": frame.get("vasp_stress"),
            "vasp_energy_per_atom_eV": frame.get("vasp_energy_per_atom_eV"),
            "vasp_max_force_eV_A": float(np.max(np.linalg.norm(
                np.asarray(frame["vasp_forces_eV_A"], dtype=float), axis=1))),
            "base_ecutwfc": float(job["ecutwfc"]), "base_ecutrho": float(job["ecutrho"]),
            "base_kpts": tuple(int(x) for x in job["kpts"]),
            "k_shift": tuple(int(x) for x in job.get("k_shift", (0, 0, 0))),
            "mp_bandgap_eV": job.get("mp_bandgap_eV"),
        })
    return out


def variants_for(p, frame: dict) -> list[dict]:
    st = frame["structure"]
    elements = sorted(str(e) for e in st.composition.elements)
    variants = []

    def add(name, pseudo="sssp_eff", kpts=None, ecut_scale=1.0, smearing="gaussian",
            occupations="smearing", degauss=DEGAUSS):
        pmap, pseudo_dir = _pseudo(p, pseudo, elements)
        base_ecut = max(pmap[e]["cutoff_wfc"] for e in elements)
        base_rho = max(pmap[e]["cutoff_rho"] for e in elements)
        if pseudo == "sssp_eff":
            base_ecut, base_rho = frame["base_ecutwfc"], frame["base_ecutrho"]
        variants.append({
            "variant": name, "pseudo": pseudo, "pseudo_dir": str(pseudo_dir), "pmap": pmap,
            "ecutwfc": round(base_ecut * ecut_scale, 1), "ecutrho": round(base_rho * ecut_scale, 1),
            "kpts": list(kpts if kpts is not None else frame["base_kpts"]),
            "k_shift": list(frame["k_shift"]), "occupations": occupations,
            "smearing": smearing, "degauss": degauss,
        })

    add("baseline")
    for spacing in (0.30, 0.25, 0.20, 0.15):
        add(f"k_{spacing:.2f}", kpts=qeio.kgrid(st, spacing))
    for scale in (0.9, 1.2, 1.5):
        add(f"ecut_{scale:.1f}", ecut_scale=scale)
    add("smear_mv", smearing="mv")
    add("smear_tetra", smearing=None, occupations="tetrahedra", degauss=None)
    add("pseudo_prec", pseudo="sssp_prec")
    add("pseudo_pdojo", pseudo="pdojo")
    add("combined_prec_k020_ecut1.2", pseudo="sssp_prec", kpts=qeio.kgrid(st, 0.20),
        ecut_scale=1.2)
    return variants


def prepare(args, p, mptraj: Path) -> list[dict]:
    out = Path(args.out) if args.out else p.work / "force_matrix"
    out.mkdir(parents=True, exist_ok=True)
    frames = select_frames(args, p, mptraj)
    jobs = []
    for frame in frames:
        for var in variants_for(p, frame):
            d = out / frame["mp_id"] / frame["sid"] / var["variant"]
            d.mkdir(parents=True, exist_ok=True)
            qeio.write_qe_input(
                d / "scf.in", frame["structure"],
                f"mat_{frame['sid']}_{var['variant']}".replace("-", "m"), "scf",
                var["ecutwfc"], var["ecutrho"], tuple(var["kpts"]), var["pmap"],
                pseudo_dir=var["pseudo_dir"], occupations=var["occupations"],
                smearing=var["smearing"] or "gaussian",
                degauss=var["degauss"] if var["degauss"] is not None else 0.02,
                conv_thr=1.0e-8, mixing_beta=0.7, electron_maxstep=200,
                k_shift=tuple(var["k_shift"]),
            )
            jobs.append({
                "mp_id": frame["mp_id"], "formula": frame["formula"],
                "traj_id": frame["traj_id"], "sid": frame["sid"],
                "variant": var["variant"], "dir": str(d), "pseudo": var["pseudo"],
                "ecutwfc": var["ecutwfc"], "ecutrho": var["ecutrho"],
                "kpts": var["kpts"], "k_shift": var["k_shift"],
                "occupations": var["occupations"], "smearing": var["smearing"],
                "vasp_forces": frame["vasp_forces"].tolist(),
                "vasp_stress": frame["vasp_stress"],
                "vasp_energy_per_atom_eV": frame["vasp_energy_per_atom_eV"],
                "vasp_max_force_eV_A": frame["vasp_max_force_eV_A"],
                "mp_bandgap_eV": frame["mp_bandgap_eV"],
            })
    (out / "jobs.json").write_text(json.dumps(jobs, indent=2))
    print(f"prepared {len(jobs)} matrix jobs for {len(frames)} frames in {out}")
    for f in frames:
        print(f["mp_id"], f["formula"], f["sid"], f["vasp_max_force_eV_A"])
    return jobs


def run_job(job: dict, gpu, ctx: dict) -> dict:
    d = Path(job["dir"])
    out = d / "scf.out"
    res = dict(job)
    res["gpu"] = gpu
    if not runner.is_done(out):
        runner.prepare_tmp(d)
        run = runner.run_pw(d / "scf.in", out, pw=ctx["pw"], mpirun=ctx["mpirun"],
                            np=ctx["np"], env=ctx["env"], gpu=gpu, timeout=ctx["timeout"])
        res["exit"] = run["exit"]
        if run["exit"] != 0 or not run["done"]:
            res["status"] = "scf_failed"
            return res
    parsed = parse.parse_scf_text(out, len(job["vasp_forces"]))
    qe_forces = parsed["forces"]
    if len(qe_forces) != len(job["vasp_forces"]):
        res["status"] = f"no_forces ({len(qe_forces)} vs {len(job['vasp_forces'])})"
        return res
    vasp_forces = np.asarray(job["vasp_forces"], dtype=float)
    diff = qe_forces - vasp_forces
    nat = len(qe_forces)
    res.update({
        "nat": nat,
        "qe_total_energy_eV": parsed["energy_eV"],
        "qe_forces": qe_forces.tolist(),
        "force_mae_eV_A": float(np.mean(np.abs(diff))),
        "force_rmse_eV_A": float(np.sqrt(np.mean(diff ** 2))),
        "force_max_err_eV_A": float(np.max(np.abs(diff))),
        "qe_max_force_eV_A": float(np.max(np.linalg.norm(qe_forces, axis=1))),
        "status": "ok",
    })
    if job.get("vasp_stress") is not None and parsed.get("stress_eV_A3") is not None:
        try:
            qe_mat = -parsed["stress_eV_A3"]
            vasp_mat = np.asarray(job["vasp_stress"], dtype=float) * KBAR_TO_EV_A3
            res["stress_MAE_eV_A3"] = float(np.mean(np.abs(qe_mat - vasp_mat)))
        except Exception:  # noqa: BLE001
            res["stress_MAE_eV_A3"] = None
    return res


def summarize(records: list[dict]) -> dict:
    ok = [r for r in records if r.get("status") == "ok"]
    out = {"n_jobs": len(records), "n_ok": len(ok), "variants": {}, "frames": {}, "factor_effects": {}}
    for variant in sorted({r["variant"] for r in ok}):
        rows = [r for r in ok if r["variant"] == variant]
        out["variants"][variant] = {
            "n": len(rows),
            "force_MAE_eV_A": float(np.mean([r["force_mae_eV_A"] for r in rows])),
            "force_RMSE_eV_A": float(np.sqrt(np.mean([r["force_rmse_eV_A"] ** 2 for r in rows]))),
            "force_max_err_eV_A": float(np.max([r["force_max_err_eV_A"] for r in rows])),
        }
    for sid in sorted({r["sid"] for r in ok}):
        rows = [r for r in ok if r["sid"] == sid]
        baseline = next((r for r in rows if r["variant"] == "baseline"), None)
        best = min(rows, key=lambda r: r["force_mae_eV_A"])
        out["frames"][sid] = {
            "formula": rows[0]["formula"], "mp_id": rows[0]["mp_id"],
            "baseline_MAE_eV_A": baseline["force_mae_eV_A"] if baseline else None,
            "best_variant": best["variant"], "best_MAE_eV_A": best["force_mae_eV_A"],
        }
    groups = {
        "k": ["k_0.30", "k_0.25", "k_0.20", "k_0.15"],
        "ecut": ["ecut_0.9", "ecut_1.2", "ecut_1.5"],
        "smearing": ["smear_mv", "smear_tetra"],
        "pseudo": ["pseudo_prec", "pseudo_pdojo"],
        "combined": ["combined_prec_k020_ecut1.2"],
    }
    for group, names in groups.items():
        out["factor_effects"][group] = {}
        for name in names:
            rows = [r for r in ok if r["variant"] == name]
            if not rows:
                continue
            out["factor_effects"][group][name] = {
                "n": len(rows),
                "force_MAE_eV_A": float(np.mean([r["force_mae_eV_A"] for r in rows])),
                "delta_vs_baseline_eV_A": float(np.mean([
                    r["force_mae_eV_A"] - next(
                        (b["force_mae_eV_A"] for b in ok
                         if b["sid"] == r["sid"] and b["variant"] == "baseline"),
                        r["force_mae_eV_A"])
                    for r in rows
                ])),
            }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--mptraj-dir", default=None)
    ap.add_argument("--frames-csv", default=None)
    ap.add_argument("--sids", default="")
    ap.add_argument("--top-n", type=int, default=5)
    ap.add_argument("--out", default=None)
    ap.add_argument("--tag", default="force_matrix")
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    ap.add_argument("--gpu-mode", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--np", type=int, default=1)
    args = ap.parse_args()

    p = get_paths(args.config)
    mptraj = Path(args.mptraj_dir) if args.mptraj_dir else p.work / "mptraj_20"
    jobs = prepare(args, p, mptraj)
    if args.prepare_only:
        return
    kind = "gpu" if args.gpu_mode == "gpu" else "cpu"
    ctx = {"pw": p.pw(kind), "mpirun": p.mpirun(), "np": args.np,
           "env": p.env(), "timeout": 6 * 3600}
    devices = [int(x) for x in args.devices.split(",")]
    gpu_values = devices if args.gpu_mode == "gpu" else [None]
    records = []
    print(f"running {len(jobs)} matrix jobs ({args.gpu_mode})", flush=True)
    with ThreadPoolExecutor(max_workers=len(gpu_values)) as ex:
        futs = {ex.submit(run_job, job, gpu_values[i % len(gpu_values)], ctx): job
                for i, job in enumerate(jobs)}
        for fut in as_completed(futs):
            job = futs[fut]
            try:
                res = fut.result()
            except Exception as exc:  # noqa: BLE001
                res = dict(job)
                res["status"] = f"exception: {exc}"
            records.append(res)
            print(job["sid"], job["variant"], res.get("status"),
                  res.get("force_mae_eV_A"), res.get("force_max_err_eV_A"), flush=True)
    records.sort(key=lambda r: (r["mp_id"], r["variant"]))
    p.results.mkdir(parents=True, exist_ok=True)
    (p.results / f"{args.tag}_frames.json").write_text(json.dumps(records, indent=2))
    pd.DataFrame([{k: v for k, v in r.items() if k not in {"vasp_forces", "qe_forces"}}
                  for r in records]).to_csv(p.results / f"{args.tag}_frames.csv", index=False)
    summary = summarize(records)
    (p.results / f"{args.tag}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
