#!/usr/bin/env python3
"""Re-run selected outliers with the SSSP precision pseudopotentials."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _bootstrap import add_config_arg, get_paths  # noqa: E402
from mp20_prepare import kgrid, load_sssp, write_qe_input, z_valence  # noqa: E402
from qe_bench_core import runner  # noqa: E402
from qe_bench_core.parse import gap_from_scf  # noqa: E402
from qe_bench_core.structure import read_structure, structure_from_ase  # noqa: E402

SYSTEMS = {"mp-2172": "AlAs", "mp-22925": "AgI"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--devices", default="0,1")
    ap.add_argument("--gpu-mode", choices=["gpu", "cpu"], default="gpu")
    args = ap.parse_args()

    from ase.io import read

    p = get_paths(args.config)
    sssp = load_sssp(p.pseudo_json("precision"))
    pseudo_dir = p.pseudo_dir("sssp_precision")
    data_dir = p.skill_root / "data" / "mp20_20"
    manifest = json.loads((data_dir / "manifest.json").read_text())["systems"]
    out = p.work / "precision"
    kind = "gpu" if args.gpu_mode == "gpu" else "cpu"
    devices = [int(x) for x in args.devices.split(",")]
    gpu_values = devices if args.gpu_mode == "gpu" else [None]

    tasks = []
    for i, (mp_id, formula) in enumerate(SYSTEMS.items()):
        meta = manifest[mp_id]
        structure = read_structure(data_dir / f"{mp_id}_{formula}.cif")
        ecutwfc, ecutrho = float(meta["ecutwfc"]), float(meta["ecutrho"])
        k_relax = kgrid(structure, spacing=0.25)
        k_nscf = kgrid(structure, spacing=0.15)
        nelec = sum(z_valence(sssp[str(s)]["filename"], pseudo_dir) for s in structure.species)
        nbnd = int(round(nelec)) // 2 + 10
        d = out / f"{mp_id}_{formula}"
        d.mkdir(parents=True, exist_ok=True)
        write_qe_input(d / "vc_relax.in", structure, f"{mp_id}_relax", "vc-relax",
                       ecutwfc, ecutrho, k_relax, sssp, pseudo_dir=pseudo_dir,
                       occupations="fixed")
        tasks.append({"mp_id": mp_id, "formula": formula, "dir": str(d),
                      "gpu": gpu_values[i % len(gpu_values)], "ecutwfc": ecutwfc,
                      "ecutrho": ecutrho, "k_nscf": k_nscf, "nbnd": nbnd})

    for t in tasks:
        d = Path(t["dir"])
        gpu = t["gpu"] if args.gpu_mode == "gpu" else None
        if not runner.is_done(d / "vc_relax.out"):
            runner.prepare_tmp(d)
            res = runner.run_pw(d / "vc_relax.in", d / "vc_relax.out",
                                pw=p.pw(kind), mpirun=p.mpirun(), np=1,
                                env=p.env(), gpu=gpu, timeout=6 * 3600)
            print(t["mp_id"], "vc-relax exit", res["exit"], flush=True)
        atoms = read(str(d / "vc_relax.out"), format="espresso-out")
        structure = structure_from_ase(atoms)
        write_qe_input(d / "nscf.in", structure, f"{t['mp_id']}_relax", "nscf",
                       t["ecutwfc"], t["ecutrho"], tuple(t["k_nscf"]), sssp,
                       pseudo_dir=pseudo_dir, occupations="fixed", nbnd=t["nbnd"])
        if not runner.is_done(d / "nscf.out"):
            res = runner.run_pw(d / "nscf.in", d / "nscf.out",
                                pw=p.pw(kind), mpirun=p.mpirun(), np=1,
                                env=p.env(), gpu=gpu, timeout=6 * 3600)
            print(t["mp_id"], "nscf exit", res["exit"], flush=True)
        t.update({
            "gap_eV": gap_from_scf(d / "nscf.out"),
            "volume_A3": float(structure.volume),
            "energy_Ry": runner.last_energy_ry(d / "vc_relax.out"),
        })
        t.pop("gpu", None)
    p.results.mkdir(parents=True, exist_ok=True)
    (p.results / "precision_results.json").write_text(json.dumps(tasks, indent=2))
    for t in tasks:
        print(t["mp_id"], t["formula"], "ecut", t["ecutwfc"], "gap", t["gap_eV"],
              "V", t["volume_A3"])


if __name__ == "__main__":
    main()
