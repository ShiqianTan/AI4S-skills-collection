#!/usr/bin/env python3
"""Controlled band-gap comparison on fixed WBM structures.

Modes:
  path   - SCF + high-symmetry band path (seekpath); shows path/off-path effects
  dense  - SCF + dense uniform nscf (0.08-0.10 A^-1); recommended baseline

Pseudopotential sets: sssp_eff, sssp_prec, pdojo. Functionals: pbe, pbesol,
and r2scan (libxc build + PseudoDojo required).
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core import parse, pseudos, qeio, runner  # noqa: E402
from qe_bench_core.structure import read_structure  # noqa: E402

RY_TO_EV = 13.605693122994


def _pseudo_set(p, name: str, elements: list[str]) -> tuple[dict, Path]:
    if name == "sssp_eff":
        upf = p.pseudo_dir("sssp_efficiency")
        sssp = json.loads(p.pseudo_json("efficiency").read_text())
        return pseudos.sssp_map(sssp, elements, upf), upf
    if name == "sssp_prec":
        upf = p.pseudo_dir("sssp_precision")
        sssp = json.loads(p.pseudo_json("precision").read_text())
        return pseudos.sssp_map(sssp, elements, upf), upf
    if name == "pdojo":
        upf = p.pseudo_dir("pseudodojo")
        return pseudos.pseudodojo_map(elements, upf), upf
    raise ValueError(name)


def _kpath(structure, path_file: Path | None) -> list[tuple[float, float, float]]:
    if path_file is not None and path_file.is_file():
        return [tuple(k) for k in json.loads(path_file.read_text())["kpoints"]]
    import seekpath

    path = seekpath.get_explicit_k_path(
        (structure.lattice.matrix, structure.frac_coords,
         [int(x.Z) for x in structure.species]),
        with_time_reversal=True,
    )
    return [tuple(float(x) for x in k) for k in path["explicit_kpoints_rel"]]


def prepare(args, p, wbm: Path) -> list[dict]:
    manifest = json.loads((wbm / "manifest.json").read_text())
    systems = manifest.get("systems", manifest)
    ids = [x for x in args.ids.split(",") if x] if args.ids else sorted(systems)
    out = Path(args.out) if args.out else p.work / "gap_controlled"
    out.mkdir(parents=True, exist_ok=True)
    jobs = []
    for mid in ids:
        d0 = wbm / mid
        meta_path = d0 / "meta.json"
        meta = json.loads(meta_path.read_text()) if meta_path.is_file() else systems.get(mid, {})
        cif = d0 / ("initial.cif" if args.start == "initial" else "vasp_final.cif")
        if not cif.is_file():
            cif = d0 / "vasp_final.cif"
        structure = read_structure(cif)
        elements = sorted(str(e) for e in structure.composition.elements)
        kpts = qeio.kgrid(structure, args.scf_density)
        ref_gap = meta.get("vasp_bandgap_pbe_eV", meta.get("vasp_bandgap_pbe"))
        for pseudo in [x for x in args.pseudo_sets.split(",") if x]:
            pmap, pseudo_dir = _pseudo_set(p, pseudo, elements)
            ecutwfc = max(pmap[e]["cutoff_wfc"] for e in elements)
            ecutrho = max(pmap[e]["cutoff_rho"] for e in elements)
            nelec = pseudos.n_electrons(structure, pmap)
            nocc = int(nelec // 2)
            nbnd = nocc + 12
            d = out / mid / pseudo
            d.mkdir(parents=True, exist_ok=True)
            if args.start == "initial":
                structure.to(filename=str(d / "start.cif"))
            functional = {} if args.functional == "pbe" else {
                "input_dft": "pbesol" if args.functional == "pbesol" else "r2scan"
            }
            scf_in = d / "scf.in"
            qeio.write_qe_input(
                scf_in, structure, f"gap_{mid}_{pseudo}", "scf",
                ecutwfc, ecutrho, kpts, pmap, pseudo_dir=pseudo_dir,
                occupations="smearing", smearing="gaussian",
                degauss=0.05 / RY_TO_EV, conv_thr=1.0e-8, mixing_beta=0.7,
                electron_maxstep=200, k_shift=(0, 0, 0), **functional,
            )
            if args.mode == "dense":
                calc_in = d / "nscf.in"
                qeio.write_qe_input(
                    calc_in, structure, f"gap_{mid}_{pseudo}", "nscf",
                    ecutwfc, ecutrho, qeio.kgrid(structure, args.nscf_density), pmap,
                    pseudo_dir=pseudo_dir, occupations="fixed", nbnd=nbnd,
                    conv_thr=1.0e-8, mixing_beta=0.7, electron_maxstep=200,
                    nosym=True, diagonalization="david", diago_full_acc=True,
                    diago_thr_init=1.0e-6, **functional,
                )
            else:
                calc_in = d / "bands.in"
                kpath = _kpath(structure, d0 / "kpath.json")
                qeio.write_qe_input(
                    calc_in, structure, f"gap_{mid}_{pseudo}", "bands",
                    ecutwfc, ecutrho, (1, 1, 1), pmap, pseudo_dir=pseudo_dir,
                    occupations="fixed", nbnd=nbnd, conv_thr=1.0e-8,
                    mixing_beta=0.7, electron_maxstep=200, nosym=True,
                    kpoints_path=kpath, diagonalization="david",
                    diago_full_acc=True, diago_thr_init=1.0e-6, **functional,
                )
            jobs.append({
                "material_id": mid, "formula": meta.get("formula"), "pseudo": pseudo,
                "mode": args.mode, "functional": args.functional, "dir": str(d),
                "nocc": nocc, "nbnd": nbnd, "nelec": nelec,
                "vasp_bandgap_pbe_eV": ref_gap,
                "ecutwfc": ecutwfc, "ecutrho": ecutrho, "kpts": list(kpts),
            })
    (out / "jobs.json").write_text(json.dumps(jobs, indent=2))
    print(f"prepared {len(jobs)} gap jobs ({args.mode}, {args.functional}) in {out}")
    return jobs


def run_job(job: dict, gpu, ctx: dict) -> dict:
    d = Path(job["dir"])
    res = dict(job)
    res["gpu"] = gpu
    kind = "libxc" if job["functional"] == "r2scan" else ctx["kind"]
    pw = ctx["p"][kind]
    scf_out = d / "scf.out"
    if not runner.is_done(scf_out):
        runner.prepare_tmp(d, f"tmp_{job['pseudo']}")
        run = runner.run_pw(d / "scf.in", scf_out, pw=pw, mpirun=ctx["mpirun"],
                            np=ctx["np"], env=ctx["env"], gpu=gpu, timeout=ctx["timeout"])
        res["scf_exit"] = run["exit"]
        if run["exit"] != 0 or not run["done"]:
            res["status"] = "scf_failed"
            return res
    calc_name = "nscf" if job["mode"] == "dense" else "bands"
    calc_out = d / f"{calc_name}.out"
    if not runner.is_done(calc_out):
        run = runner.run_pw(d / f"{calc_name}.in", calc_out, pw=pw, mpirun=ctx["mpirun"],
                            np=ctx["np"], env=ctx["env"], gpu=gpu, timeout=ctx["timeout"])
        res[f"{calc_name}_exit"] = run["exit"]
        if run["exit"] != 0 or not run["done"]:
            res["status"] = f"{calc_name}_failed"
            return res
    if job["mode"] == "dense":
        gap = parse.gap_from_scf(calc_out)
    else:
        gap = parse.gap_from_bands(calc_out, job["nocc"]).get("gap_eV")
    res["gap_eV"] = gap
    if gap is not None and job.get("vasp_bandgap_pbe_eV") is not None:
        res["gap_error_eV"] = gap - job["vasp_bandgap_pbe_eV"]
        res["status"] = "ok"
    else:
        res["status"] = "parse_failed"
    return res


def summarize(records: list[dict], pseudo_sets: list[str]) -> dict:
    ok = [r for r in records if r.get("status") == "ok"]
    out = {"n_jobs": len(records), "n_ok": len(ok), "per_pseudo": {}}
    for pseudo in pseudo_sets:
        rows = [r for r in ok if r["pseudo"] == pseudo and r.get("gap_error_eV") is not None]
        if not rows:
            continue
        err = np.array([r["gap_error_eV"] for r in rows], dtype=float)
        out["per_pseudo"][pseudo] = {
            "n": len(rows), "gap_MAE_eV": float(np.mean(np.abs(err))),
            "gap_RMSE_eV": float(np.sqrt(np.mean(err ** 2))),
            "gap_bias_eV": float(np.mean(err)),
            "gap_max_abs_eV": float(np.max(np.abs(err))),
        }
    if ok:
        err = np.array([r["gap_error_eV"] for r in ok if r.get("gap_error_eV") is not None], dtype=float)
        if err.size:
            out["all"] = {"n": int(err.size), "gap_MAE_eV": float(np.mean(np.abs(err))),
                          "gap_RMSE_eV": float(np.sqrt(np.mean(err ** 2))),
                          "gap_bias_eV": float(np.mean(err)),
                          "gap_max_abs_eV": float(np.max(np.abs(err)))}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--wbm-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--ids", default="")
    ap.add_argument("--start", choices=["final", "initial"], default="final")
    ap.add_argument("--mode", choices=["path", "dense"], default="dense")
    ap.add_argument("--pseudo-sets", default="sssp_eff,sssp_prec,pdojo")
    ap.add_argument("--functional", choices=["pbe", "pbesol", "r2scan"], default="pbe")
    ap.add_argument("--scf-density", type=float, default=0.20)
    ap.add_argument("--nscf-density", type=float, default=0.10)
    ap.add_argument("--tag", default="gap_controlled")
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    ap.add_argument("--gpu-mode", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--np", type=int, default=1)
    args = ap.parse_args()

    p = get_paths(args.config)
    wbm = Path(args.wbm_dir) if args.wbm_dir else p.work / "wbm_gap"
    jobs = prepare(args, p, wbm)
    if args.prepare_only:
        return
    kind = "gpu" if args.gpu_mode == "gpu" else "cpu"
    ctx = {
        "kind": kind, "p": {"gpu": p.pw("gpu"), "cpu": p.pw("cpu"), "libxc": p.pw("libxc")},
        "mpirun": p.mpirun(), "np": args.np, "env": p.env(), "timeout": 6 * 3600,
    }
    devices = [int(x) for x in args.devices.split(",")]
    gpu_values = devices if args.gpu_mode == "gpu" else [None]
    records = []
    print(f"running {len(jobs)} gap jobs ({args.gpu_mode})", flush=True)
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
            print(job["material_id"], job["pseudo"], res.get("status"),
                  res.get("vasp_bandgap_pbe_eV"), res.get("gap_eV"), res.get("gap_error_eV"),
                  flush=True)
    records.sort(key=lambda r: (r["material_id"], r["pseudo"]))
    p.results.mkdir(parents=True, exist_ok=True)
    (p.results / f"{args.tag}_frames.json").write_text(json.dumps(records, indent=2))
    import pandas as pd

    pd.DataFrame([{k: v for k, v in r.items() if k != "kpts"} for r in records]) \
        .to_csv(p.results / f"{args.tag}_frames.csv", index=False)
    summary = summarize(records, [x for x in args.pseudo_sets.split(",") if x])
    (p.results / f"{args.tag}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
