#!/usr/bin/env python3
"""Run the MP20 QE workflow (vc-relax + optional fixed SCF + dense nscf)."""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _bootstrap import add_config_arg, get_paths  # noqa: E402
from mp20_prepare import _filename_map, load_sssp, write_qe_input  # noqa: E402
from qe_bench_core import parse, qeio, runner  # noqa: E402
from qe_bench_core.structure import structure_from_ase  # noqa: E402


def _read_relaxed(path: Path):
    from ase.io import read

    return read(str(path), format="espresso-out")


def relax_and_nscf(mp_id: str, meta: dict, gpu, force: bool, out_dir: Path, ctx: dict) -> dict:
    system_dir = out_dir / f"{mp_id}_{meta['formula']}"
    result = {"mp_id": mp_id, "formula": meta["formula"], "gpu": gpu}

    if force or not runner.is_done(system_dir / "vc_relax.out"):
        runner.prepare_tmp(system_dir)
        res = runner.run_pw(
            system_dir / "vc_relax.in", system_dir / "vc_relax.out",
            pw=ctx["pw"], mpirun=ctx["mpirun"], np=ctx["np"],
            env=ctx["base_env"], gpu=gpu, timeout=ctx["timeout"],
        )
        result["vc_relax_exit"] = res["exit"]
        result["vc_relax_wall_s"] = res["wall_s"]
        if res["exit"] != 0 or not res["done"]:
            result["status"] = "vc_relax_failed"
            return result
    else:
        result["vc_relax_exit"] = 0
        result["vc_relax_wall_s"] = None

    atoms = _read_relaxed(system_dir / "vc_relax.out")
    structure = structure_from_ase(atoms)
    result["relaxed_volume_A3"] = round(float(structure.volume), 6)
    result["relaxed_abc"] = [round(x, 6) for x in structure.lattice.abc]
    result["relaxed_angles"] = [round(x, 6) for x in structure.lattice.angles]
    result["vc_relax_energy_Ry"] = runner.last_energy_ry(system_dir / "vc_relax.out")

    pmap = _filename_map(ctx["sssp"], [str(e) for e in structure.composition.elements])

    if meta.get("pre_nscf_scf", False) and not meta.get("metal", False):
        write_qe_input(
            system_dir / "scf_fixed.in", structure, f"{mp_id}_relax", "scf",
            meta["ecutwfc"], meta["ecutrho"], tuple(meta.get("k_scf", meta["k_relax"])),
            pmap, pseudo_dir=ctx["pseudo_dir"], occupations="fixed",
            conv_thr=1.0e-8, diagonalization=meta.get("nscf_diagonalization"),
            k_shift=tuple(meta.get("k_shift", (0, 0, 0))),
        )
        if force or not runner.is_done(system_dir / "scf_fixed.out"):
            res = runner.run_pw(
                system_dir / "scf_fixed.in", system_dir / "scf_fixed.out",
                pw=ctx["pw"], mpirun=ctx["mpirun"], np=ctx["np"],
                env=ctx["base_env"], gpu=gpu, timeout=ctx["timeout"],
            )
            result["scf_fixed_exit"] = res["exit"]
            result["scf_fixed_wall_s"] = res["wall_s"]
            if res["exit"] != 0 or not res["done"]:
                result["status"] = "scf_fixed_failed"
                return result

    write_qe_input(
        system_dir / "nscf_run.in", structure, f"{mp_id}_relax", "nscf",
        meta["ecutwfc"], meta["ecutrho"], tuple(meta["k_nscf"]), pmap,
        pseudo_dir=ctx["pseudo_dir"],
        occupations=meta.get("nscf_occupations", meta.get("occupations", "fixed")),
        nbnd=meta["nbnd"], diagonalization=meta.get("nscf_diagonalization"),
    )

    nscf_out = system_dir / "nscf.out"
    if force or not runner.is_done(nscf_out):
        res = runner.run_pw(
            system_dir / "nscf_run.in", nscf_out,
            pw=ctx["pw"], mpirun=ctx["mpirun"], np=ctx["np"],
            env=ctx["base_env"], gpu=gpu, timeout=ctx["timeout"],
        )
        result["nscf_exit"] = res["exit"]
        result["nscf_wall_s"] = res["wall_s"]
        if res["exit"] != 0 or not res["done"]:
            result["status"] = "nscf_failed"
            return result
    else:
        result["nscf_exit"] = 0
        result["nscf_wall_s"] = None

    homo = lumo = None
    for line in nscf_out.read_text(errors="ignore").splitlines():
        if "highest occupied, lowest unoccupied level" in line:
            vals = line.split(":")[1].split()
            homo, lumo = float(vals[0]), float(vals[1])
        elif "highest occupied level" in line and "lowest unoccupied" not in line:
            homo = float(line.split(":")[1].split()[0])
    result["homo_eV"] = homo
    result["lumo_eV"] = lumo
    result["qe_band_gap_eV"] = round(lumo - homo, 4) if homo is not None and lumo is not None else 0.0
    result["mp_band_gap_eV"] = meta["mp_band_gap"]
    result["gap_error_eV"] = round(result["qe_band_gap_eV"] - meta["mp_band_gap"], 4)
    result["mp_formation_energy_per_atom"] = meta["mp_formation_energy_per_atom"]
    result["cif_volume_A3"] = meta["cif_volume"]
    result["volume_error_pct"] = round(
        100.0 * (result["relaxed_volume_A3"] - meta["cif_volume"]) / meta["cif_volume"], 4
    )
    result["status"] = "ok"
    return result


def run_reference(el: str, meta: dict, gpu, force: bool, out_dir: Path, ctx: dict) -> dict:
    d = out_dir / "references" / el
    outfile = d / "vc_relax.out"
    result = {"element": el, "mp_id": meta["mp_id"], "gpu": gpu}
    if force or not runner.is_done(outfile):
        runner.prepare_tmp(d)
        res = runner.run_pw(
            d / "vc_relax.in", outfile, pw=ctx["pw"], mpirun=ctx["mpirun"],
            np=ctx["np"], env=ctx["base_env"], gpu=gpu, timeout=ctx["timeout"],
        )
        result["exit"] = res["exit"]
        result["wall_s"] = res["wall_s"]
        if res["exit"] != 0 or not res["done"]:
            result["status"] = "failed"
            return result
    else:
        result["exit"] = 0
        result["wall_s"] = None
    atoms = _read_relaxed(outfile)
    result["nat"] = len(atoms)
    result["energy_Ry"] = runner.last_energy_ry(outfile)
    result["energy_per_atom_Ry"] = result["energy_Ry"] / len(atoms)
    result["volume_A3"] = round(float(atoms.get_volume()), 6)
    result["status"] = "ok"
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--skip-refs", action="store_true")
    ap.add_argument("--mp20-dir", default=None)
    ap.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    ap.add_argument("--gpu-mode", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--timeout", type=int, default=6 * 3600)
    args = ap.parse_args()

    p = get_paths(args.config)
    out_dir = Path(args.mp20_dir) if args.mp20_dir else p.work / "mp20_20"
    manifest = json.loads((out_dir / "manifest.json").read_text())
    kind = "gpu" if args.gpu_mode == "gpu" else "cpu"
    ctx = {
        "pw": p.pw(kind),
        "mpirun": p.mpirun(),
        "np": args.np,
        "base_env": p.env(),
        "pseudo_dir": p.pseudo_dir("sssp_efficiency"),
        "sssp": load_sssp(p.pseudo_json("efficiency")),
        "timeout": args.timeout,
    }
    devices = [int(x) for x in args.devices.split(",")]
    gpu_values = devices if args.gpu_mode == "gpu" else [None]

    systems = list(manifest["systems"])
    results = {}
    print(f"running {len(systems)} MP20 systems ({args.gpu_mode})", flush=True)
    with ThreadPoolExecutor(max_workers=len(gpu_values)) as ex:
        futs = {
            ex.submit(
                relax_and_nscf, mp_id, manifest["systems"][mp_id],
                gpu_values[i % len(gpu_values)], args.force, out_dir, ctx,
            ): mp_id
            for i, mp_id in enumerate(systems)
        }
        for fut in as_completed(futs):
            mp_id = futs[fut]
            try:
                res = fut.result()
            except Exception as exc:  # noqa: BLE001
                res = {"mp_id": mp_id, "status": f"exception: {exc}"}
            results[mp_id] = res
            print(mp_id, res.get("status"), res.get("qe_band_gap_eV"), flush=True)

    refs = {}
    if not args.skip_refs:
        print("running elemental references", flush=True)
        for i, (el, meta) in enumerate(manifest.get("references", {}).items()):
            try:
                res = run_reference(
                    el, meta, gpu_values[i % len(gpu_values)], args.force, out_dir, ctx
                )
            except Exception as exc:  # noqa: BLE001
                res = {"element": el, "status": f"exception: {exc}"}
            refs[el] = res
            print(el, res.get("status"), flush=True)

    payload = {"systems": results, "references": refs}
    (out_dir / "run_status.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
