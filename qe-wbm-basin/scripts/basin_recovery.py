#!/usr/bin/env python3
"""WBM basin-recovery benchmark.

For each WBM system, vc-relax is started from
  * the VASP final structure,
  * the WBM initial structure,
  * a perturbed VASP final structure.
For anomalous systems, additional perturbations and a volume scan are added.
Recovery is measured against the VASP final structure with StructureMatcher.
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
from qe_bench_core import mprelax, qeio, runner  # noqa: E402

RY_TO_EV = 13.605693122994
DEG = 0.05 / RY_TO_EV
IMP = dict(
    smearing="gaussian", degauss=DEG, conv_thr=1.0e-6, mixing_beta=0.7,
    electron_maxstep=100, diagonalization="david", diago_thr_init=1.0e-6,
    diago_full_acc=True, nstep=60, forc_conv_thr=5.0e-4, etot_conv_thr=1.0e-5,
    nosym=True, ion_dynamics="bfgs", cell_dynamics="bfgs", press_conv_thr=0.5,
)


def case_structures(mid: str, wbm: Path, anomalies: set[str]) -> dict:
    from pymatgen.core import Structure

    d0 = wbm / mid
    vasp = Structure.from_file(d0 / "vasp_final.cif")
    init = Structure.from_file(d0 / "initial.cif")
    rng = np.random.default_rng(abs(hash(mid)) % (2 ** 32))
    pert = vasp.copy().scale_lattice(vasp.volume * 1.03)
    disp = rng.normal(0.0, 0.08, size=(len(pert), 3))
    for i, site in enumerate(pert):
        pert.translate_sites(i, disp[i], frac_coords=False)
    cases = {"vaspfinal": vasp, "initial": init, "perturb": pert}
    if mid in anomalies:
        rng2 = np.random.default_rng((abs(hash(mid)) + 17) % (2 ** 32))
        pert2 = Structure.from_file(d0 / "vasp_final.cif").scale_lattice(vasp.volume * 0.97)
        disp2 = rng2.normal(0.0, 0.12, size=(len(pert2), 3))
        for i, site in enumerate(pert2):
            pert2.translate_sites(i, disp2[i], frac_coords=False)
        cases["perturb2"] = pert2
        for factor in (0.95, 0.97, 1.03, 1.05):
            s = Structure.from_file(d0 / "vasp_final.cif")
            cases[f"vol_{factor:.2f}"] = s.scale_lattice(s.volume * factor)
    return cases


def prepare(args, p, wbm: Path, out: Path) -> list[dict]:
    from pymatgen.core import Structure

    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((wbm / "manifest.json").read_text())["systems"]
    sssp = json.loads(p.pseudo_json("efficiency").read_text())
    pseudo_dir = p.pseudo_dir("sssp_efficiency")
    anomalies = {x for x in args.anomalies.split(",") if x}
    materials = {x for x in args.materials.split(",") if x}
    cases_filter = {x for x in args.cases.split(",") if x}
    jobs = []
    for mid, meta in manifest.items():
        if materials and mid not in materials:
            continue
        d0 = wbm / mid
        vasp = Structure.from_file(d0 / "vasp_final.cif")
        ecut_i, ecutrho_i = mprelax.mp_cutoffs(vasp, sssp)
        ecut = max(45.0, round(ecut_i * args.ecut_scale, 1))
        ecutrho = max(8.0 * ecut, round(ecutrho_i * args.ecut_scale, 1))
        kpts = qeio.kgrid(vasp, args.k_spacing)
        for case, structure in case_structures(mid, wbm, anomalies).items():
            if cases_filter and case not in cases_filter:
                continue
            d = out / mid / case
            d.mkdir(parents=True, exist_ok=True)
            structure.to(filename=str(d / "start.cif"))
            prefix = f"basin_{mid}_{case}".replace("-", "m")
            imp = dict(IMP)
            imp.update(nstep=args.nstep, forc_conv_thr=args.forc_conv,
                       press_conv_thr=args.press_conv, conv_thr=args.conv_thr,
                       mixing_beta=args.mixing_beta,
                       electron_maxstep=args.electron_maxstep,
                       diagonalization=args.diagonalization)
            if args.mixing_mode:
                imp["mixing_mode"] = args.mixing_mode
            qeio.write_qe_input(d / "vc_relax.in", structure, prefix, "vc-relax",
                                ecut, ecutrho, kpts, sssp, pseudo_dir=pseudo_dir, **imp)
            jobs.append({
                "material_id": mid, "formula": meta["formula"], "case": case,
                "dir": str(d), "nat": len(structure),
                "start_volume_A3": float(structure.volume),
                "vasp_volume_A3": float(vasp.volume),
                "ecutwfc": ecut, "ecutrho": ecutrho, "kpts": list(kpts),
            })
    (out / "jobs.json").write_text(json.dumps(jobs, indent=2))
    print(f"prepared {len(jobs)} basin jobs for {len(manifest)} systems", flush=True)
    return jobs


def run_job(job: dict, gpu, ctx: dict) -> dict:
    d = Path(job["dir"])
    out = d / "vc_relax.out"
    res = dict(job)
    res["gpu"] = gpu
    if not runner.is_done(out):
        runner.prepare_tmp(d)
        run = runner.run_pw(d / "vc_relax.in", out, pw=ctx["pw"], mpirun=ctx["mpirun"],
                            np=ctx["np"], env=ctx["env"], gpu=gpu, timeout=ctx["timeout"])
        res["exit"], res["wall_s"] = run["exit"], run["wall_s"]
        if run["exit"] != 0 or not run["done"]:
            res["status"] = "vc_failed"
            return res
    else:
        res["exit"], res["wall_s"] = 0, None
    try:
        from ase.io import read
        from pymatgen.io.ase import AseAtomsAdaptor

        atoms = read(str(out), format="espresso-out")
        qe = AseAtomsAdaptor.get_structure(atoms)
        energy = float(atoms.get_potential_energy())
    except Exception as exc:  # noqa: BLE001
        res["status"] = f"parse_failed: {exc}"
        return res
    res.update(_compare(qe, ctx["wbm"] / job["material_id"] / "vasp_final.cif", energy))
    res["status"] = "ok"
    return res


def _compare(qe, vasp_path: Path, energy: float) -> dict:
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    vasp = Structure.from_file(vasp_path)
    out = {
        "qe_energy_eV": energy,
        "qe_energy_per_atom_eV": energy / len(qe),
        "qe_volume_A3": float(qe.volume),
        "volume_diff_pct": float(100.0 * (qe.volume - vasp.volume) / vasp.volume),
        "lattice_abc_mae_A": float(np.mean(np.abs(np.array(qe.lattice.abc) - np.array(vasp.lattice.abc)))),
        "lattice_angle_mae_deg": float(np.mean(np.abs(np.array(qe.lattice.angles) - np.array(vasp.lattice.angles)))),
        "qe_spacegroup": int(SpacegroupAnalyzer(qe, symprec=0.1).get_space_group_number()),
        "vasp_spacegroup": int(SpacegroupAnalyzer(vasp, symprec=0.1).get_space_group_number()),
    }
    try:
        sm = StructureMatcher(primitive_cell=False, scale=False, attempt_supercell=False,
                              ltol=0.2, stol=0.3, angle_tol=5)
        out["structure_match"] = bool(sm.fit(qe, vasp))
        rms = sm.get_rms_dist(qe, vasp)
        out["structure_rms"] = float(rms[0]) if rms else None
        out["structure_max_dist"] = float(rms[1]) if rms else None
    except Exception as exc:  # noqa: BLE001
        out["structure_match"] = None
        out["structure_error"] = str(exc)
    return out


def collect(p, out: Path, wbm: Path, tag: str) -> None:
    from ase.io import read
    from pymatgen.io.ase import AseAtomsAdaptor

    jobs = json.loads((out / "jobs.json").read_text())
    rows = []
    for job in jobs:
        d = Path(job["dir"])
        outfile = d / "vc_relax.out"
        if not runner.is_done(outfile):
            rows.append({**job, "status": "missing"})
            continue
        try:
            atoms = read(str(outfile), format="espresso-out")
            qe = AseAtomsAdaptor.get_structure(atoms)
            energy = float(atoms.get_potential_energy())
        except Exception as exc:  # noqa: BLE001
            rows.append({**job, "status": f"parse_failed: {exc}"})
            continue
        rows.append({**job, **_compare(qe, wbm / job["material_id"] / "vasp_final.cif", energy),
                     "status": "ok"})

    ok = [r for r in rows if r.get("status") == "ok"]
    by_mid: dict[str, list] = {}
    for r in ok:
        by_mid.setdefault(r["material_id"], []).append(r)
    for rs in by_mid.values():
        emin = min(r["qe_energy_per_atom_eV"] for r in rs)
        for r in rs:
            r["rel_energy_meV_atom"] = (r["qe_energy_per_atom_eV"] - emin) * 1000.0
            r["is_lowest_energy"] = r["rel_energy_meV_atom"] <= 1.0

    p.results.mkdir(parents=True, exist_ok=True)
    (p.results / f"{tag}_frames.json").write_text(json.dumps(rows, indent=2))
    import pandas as pd

    pd.DataFrame([{k: v for k, v in r.items() if k != "kpts"} for r in rows]) \
        .to_csv(p.results / f"{tag}_frames.csv", index=False)

    def rate(subset):
        vals = [bool(r["structure_match"]) for r in subset if r.get("structure_match") is not None]
        return float(np.mean(vals)) if vals else None

    summary = {"n_systems": len(by_mid), "n_jobs": len(rows), "n_ok": len(ok), "by_case": {}}
    for case in sorted({r["case"] for r in ok}):
        sub = [r for r in ok if r["case"] == case]
        summary["by_case"][case] = {
            "n": len(sub), "match_rate": rate(sub),
            "mean_volume_diff_pct": float(np.mean([r["volume_diff_pct"] for r in sub])),
            "mean_rel_energy_meV_atom": float(np.mean([r["rel_energy_meV_atom"] for r in sub])),
        }
    summary["vaspfinal_recovery_rate"] = rate([r for r in ok if r["case"] == "vaspfinal"])
    summary["initial_recovery_rate"] = rate([r for r in ok if r["case"] == "initial"])
    summary["perturb_recovery_rate"] = rate([r for r in ok if r["case"] == "perturb"])
    summary["any_start_matches_rate"] = float(np.mean([
        any(r.get("structure_match") for r in rs) for rs in by_mid.values()]))
    summary["best_energy_matches_rate"] = float(np.mean([
        any(r.get("structure_match") and r.get("is_lowest_energy") for r in rs)
        for rs in by_mid.values()]))
    (p.results / f"{tag}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--wbm-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--tag", default="wbm_basin")
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--collect-only", action="store_true")
    ap.add_argument("--cases", default="vaspfinal,initial,perturb")
    ap.add_argument("--anomalies", default="")
    ap.add_argument("--materials", default="")
    ap.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    ap.add_argument("--backend", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--k-spacing", type=float, default=0.25)
    ap.add_argument("--ecut-scale", type=float, default=1.0)
    ap.add_argument("--nstep", type=int, default=60)
    ap.add_argument("--forc-conv", type=float, default=5.0e-4)
    ap.add_argument("--press-conv", type=float, default=0.5)
    ap.add_argument("--conv-thr", type=float, default=1.0e-6)
    ap.add_argument("--mixing-mode", default="")
    ap.add_argument("--mixing-beta", type=float, default=0.7)
    ap.add_argument("--diagonalization", default="david")
    ap.add_argument("--electron-maxstep", type=int, default=100)
    args = ap.parse_args()

    p = get_paths(args.config)
    wbm = Path(args.wbm_dir) if args.wbm_dir else p.work / "wbm"
    out = Path(args.out) if args.out else p.work / "wbm_basin"
    if args.collect_only:
        collect(p, out, wbm, args.tag)
        return
    jobs = prepare(args, p, wbm, out)
    if args.prepare_only:
        return
    kind = "gpu" if args.backend == "gpu" else "cpu"
    ctx = {"pw": p.pw(kind), "mpirun": p.mpirun(), "np": args.np,
           "env": p.env(), "timeout": 6 * 3600, "wbm": wbm}
    devices = [int(x) for x in args.devices.split(",")]
    gpu_values = devices if args.backend == "gpu" else [None]
    records = []
    print(f"running {len(jobs)} basin jobs ({args.backend})", flush=True)
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
            print(job["material_id"], job["case"], res.get("status"),
                  res.get("structure_match"), res.get("volume_diff_pct"), flush=True)
    (out / "run_records.json").write_text(json.dumps(records, indent=2))
    collect(p, out, wbm, args.tag)


if __name__ == "__main__":
    main()
