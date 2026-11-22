#!/usr/bin/env python3
"""Run WBM initial-structure relaxations with QE and compare to VASP finals."""
from __future__ import annotations

import argparse
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core import runner  # noqa: E402

RY_TO_EV = 13.605693122994


def first_energy_ev(outfile: Path) -> float | None:
    text = outfile.read_text(errors="ignore")
    m = re.search(r"!\s+total energy\s+=\s+([-\d.]+)\s+Ry", text)
    return float(m.group(1)) * RY_TO_EV if m else None


def last_energy_ev(outfile: Path) -> float | None:
    vals = [float(x) for x in re.findall(r"!\s+total energy\s+=\s+([-\d.]+)\s+Ry",
                                          outfile.read_text(errors="ignore"))]
    return vals[-1] * RY_TO_EV if vals else None


def process(mid: str, meta: dict, gpu, ctx: dict) -> dict:
    from ase.io import read
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    d = Path(meta["dir"])
    vc_out, sf_out = d / "vc_relax.out", d / "scf_vasp_final.out"
    res = {"material_id": mid, **meta, "gpu": gpu}
    if not runner.is_done(vc_out):
        runner.prepare_tmp(d)
        run = runner.run_pw(d / "vc_relax.in", vc_out, pw=ctx["pw"],
                            mpirun=ctx["mpirun"], np=ctx["np"], env=ctx["env"],
                            gpu=gpu, timeout=ctx["timeout"])
        res["vc_exit"], res["vc_wall_s"] = run["exit"], run["wall_s"]
        if run["exit"] != 0 or not run["done"]:
            res["status"] = "vc_failed"
            return res
    if not runner.is_done(sf_out):
        run = runner.run_pw(d / "scf_vasp_final.in", sf_out, pw=ctx["pw"],
                            mpirun=ctx["mpirun"], np=ctx["np"], env=ctx["env"],
                            gpu=gpu, timeout=ctx["timeout"])
        res["sf_exit"], res["sf_wall_s"] = run["exit"], run["wall_s"]
        if run["exit"] != 0 or not run["done"]:
            res["status"] = "sf_failed"
            return res

    qe_final = AseAtomsAdaptor.get_structure(read(str(vc_out), format="espresso-out"))
    qe_final.to(filename=str(d / "qe_relaxed.cif"))
    vasp_final = Structure.from_file(d / "vasp_final.cif")
    scf_atoms = read(str(sf_out), format="espresso-out")
    forces = np.asarray(scf_atoms.calc.results.get("forces", []), dtype=float)
    fnorms = np.linalg.norm(forces, axis=1) if forces.size else np.array([])

    res.update({
        "qe_final_energy_eV": last_energy_ev(vc_out),
        "qe_initial_energy_eV": first_energy_ev(vc_out),
        "qe_energy_at_vasp_final_eV": last_energy_ev(sf_out),
        "qe_final_volume_A3": float(qe_final.volume),
        "qe_final_abc": list(qe_final.lattice.abc),
        "qe_final_angles": list(qe_final.lattice.angles),
        "vasp_final_volume_A3": float(vasp_final.volume),
        "vasp_final_abc": list(vasp_final.lattice.abc),
        "vasp_final_angles": list(vasp_final.lattice.angles),
        "volume_diff_pct": float(100.0 * (qe_final.volume - vasp_final.volume) / vasp_final.volume),
        "lattice_abc_mae": float(np.mean(np.abs(np.array(qe_final.lattice.abc) - np.array(vasp_final.lattice.abc)))),
        "lattice_angle_mae": float(np.mean(np.abs(np.array(qe_final.lattice.angles) - np.array(vasp_final.lattice.angles)))),
        "qe_final_spacegroup": int(SpacegroupAnalyzer(qe_final, symprec=0.1).get_space_group_number()),
        "vasp_final_spacegroup": int(SpacegroupAnalyzer(vasp_final, symprec=0.1).get_space_group_number()),
        "qe_force_max_at_vasp_final_eV_A": float(fnorms.max()) if fnorms.size else None,
        "qe_force_mean_at_vasp_final_eV_A": float(fnorms.mean()) if fnorms.size else None,
        "status": "ok",
    })
    try:
        sm = StructureMatcher(primitive_cell=False, scale=False, attempt_supercell=False,
                              ltol=0.2, stol=0.3, angle_tol=5)
        fit = sm.fit(qe_final, vasp_final)
        rms = sm.get_rms_dist(qe_final, vasp_final)
        res["structure_match"] = bool(fit)
        res["structure_rms"] = float(rms[0]) if rms else None
        res["structure_max_dist"] = float(rms[1]) if rms else None
    except Exception as exc:  # noqa: BLE001
        res["structure_match"] = None
        res["structure_error"] = str(exc)
    try:
        res["qe_relaxation_energy_eV_atom"] = (
            (res["qe_initial_energy_eV"] - res["qe_final_energy_eV"]) / len(qe_final))
        res["qe_vaspfinal_energy_diff_meV_atom"] = (
            (res["qe_energy_at_vasp_final_eV"] - res["qe_final_energy_eV"]) / len(qe_final) * 1000.0)
    except Exception:  # noqa: BLE001
        pass
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--wbm-dir", default=None)
    ap.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    ap.add_argument("--gpu-mode", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--tag", default="wbm")
    args = ap.parse_args()

    p = get_paths(args.config)
    wbm = Path(args.wbm_dir) if args.wbm_dir else p.work / "wbm"
    kind = "gpu" if args.gpu_mode == "gpu" else "cpu"
    ctx = {"pw": p.pw(kind), "mpirun": p.mpirun(), "np": args.np,
           "env": p.env(), "timeout": 6 * 3600}
    manifest = json.loads((wbm / "manifest.json").read_text())["systems"]
    devices = [int(x) for x in args.devices.split(",")]
    gpu_values = devices if args.gpu_mode == "gpu" else [None]

    records = []
    print(f"running {len(manifest)} WBM systems ({args.gpu_mode})", flush=True)
    with ThreadPoolExecutor(max_workers=len(gpu_values)) as ex:
        futs = {ex.submit(process, mid, meta, gpu_values[i % len(gpu_values)], ctx): mid
                for i, (mid, meta) in enumerate(manifest.items())}
        for fut in as_completed(futs):
            mid = futs[fut]
            try:
                res = fut.result()
            except Exception as exc:  # noqa: BLE001
                res = {"material_id": mid, "status": f"exception: {exc}"}
            records.append(res)
            print(mid, res.get("status"), res.get("volume_diff_pct"),
                  res.get("structure_rms"), res.get("qe_force_max_at_vasp_final_eV_A"), flush=True)
    records.sort(key=lambda r: r["material_id"])
    p.results.mkdir(parents=True, exist_ok=True)
    (p.results / f"{args.tag}_frames.json").write_text(json.dumps(records, indent=2))
    import pandas as pd

    pd.DataFrame([{k: v for k, v in r.items() if k != "potcar_symbols"} for r in records]) \
        .to_csv(p.results / f"{args.tag}_frames.csv", index=False)

    ok = [r for r in records if r.get("status") == "ok"]
    summary = {
        "n_systems": len(records), "n_ok": len(ok),
        "volume_MAE_pct": float(np.mean([abs(r["volume_diff_pct"]) for r in ok])) if ok else None,
        "volume_max_abs_pct": float(np.max([abs(r["volume_diff_pct"]) for r in ok])) if ok else None,
        "lattice_abc_MAE_A": float(np.mean([r["lattice_abc_mae"] for r in ok])) if ok else None,
        "lattice_angle_MAE_deg": float(np.mean([r["lattice_angle_mae"] for r in ok])) if ok else None,
        "structure_match_rate": float(np.mean([bool(r.get("structure_match")) for r in ok])) if ok else None,
        "structure_rms_median": float(np.median([r["structure_rms"] for r in ok if r.get("structure_rms") is not None])) if ok else None,
        "qe_force_max_at_vasp_final_MAE_eV_A": float(np.mean([r["qe_force_max_at_vasp_final_eV_A"] for r in ok])) if ok else None,
        "qe_vaspfinal_energy_diff_MAE_meV_atom": float(np.mean([abs(r["qe_vaspfinal_energy_diff_meV_atom"]) for r in ok if "qe_vaspfinal_energy_diff_meV_atom" in r])) if ok else None,
    }
    (p.results / f"{args.tag}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
