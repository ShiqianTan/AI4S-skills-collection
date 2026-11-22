#!/usr/bin/env python3
"""k-point and cutoff convergence test for selected MP20 systems."""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _bootstrap import add_config_arg, get_paths  # noqa: E402
from mp20_prepare import _filename_map, load_sssp, write_qe_input  # noqa: E402
from qe_bench_core import parse, runner  # noqa: E402
from qe_bench_core.structure import structure_from_ase  # noqa: E402

DEFAULT_SYSTEMS = {
    "mp-2172": {"formula": "AlAs", "k_list": [6, 8, 10, 12, 14], "e_list": [35, 42, 49]},
    "mp-22925": {"formula": "AgI", "k_list": [6, 8, 10, 12], "e_list": None},
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--mp20-dir", default=None)
    ap.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    ap.add_argument("--gpu-mode", choices=["gpu", "cpu"], default="gpu")
    args = ap.parse_args()

    from ase.io import read

    p = get_paths(args.config)
    mp20 = Path(args.mp20_dir) if args.mp20_dir else p.work / "mp20_20"
    manifest = json.loads((mp20 / "manifest.json").read_text())
    sssp = load_sssp(p.pseudo_json("efficiency"))
    pseudo_dir = p.pseudo_dir("sssp_efficiency")
    kind = "gpu" if args.gpu_mode == "gpu" else "cpu"
    devices = [int(x) for x in args.devices.split(",")]
    gpu_values = devices if args.gpu_mode == "gpu" else [None]

    cases = []
    for mp_id, cfg in DEFAULT_SYSTEMS.items():
        if mp_id not in manifest["systems"]:
            continue
        meta = manifest["systems"][mp_id]
        d = mp20 / f"{mp_id}_{cfg['formula']}"
        atoms = read(str(d / "vc_relax.out"), format="espresso-out")
        structure = structure_from_ase(atoms)
        pmap = _filename_map(sssp, [str(e) for e in structure.composition.elements])
        ref_e = meta["ecutwfc"]
        e_list = cfg["e_list"] or [round(ref_e * f, 1) for f in (1.0, 1.15, 1.30)]
        for k in cfg["k_list"]:
            tag = f"conv_k{k}_e{ref_e:.0f}"
            write_qe_input(d / f"{tag}.in", structure, f"{mp_id}_relax", "nscf",
                           ref_e, meta["ecutrho"], (k, k, k), pmap,
                           pseudo_dir=pseudo_dir, occupations=meta.get("occupations", "fixed"),
                           nbnd=meta["nbnd"])
            cases.append({"mp_id": mp_id, "formula": cfg["formula"], "kind": "k",
                          "k": k, "ecut": ref_e, "tag": tag, "dir": str(d)})
        for e in e_list:
            ratio = e / ref_e
            tag = f"conv_k12_e{e:.0f}"
            write_qe_input(d / f"{tag}.in", structure, f"{mp_id}_relax", "nscf",
                           e, meta["ecutrho"] * ratio, (12, 12, 12), pmap,
                           pseudo_dir=pseudo_dir, occupations=meta.get("occupations", "fixed"),
                           nbnd=meta["nbnd"])
            cases.append({"mp_id": mp_id, "formula": cfg["formula"], "kind": "ecut",
                          "k": 12, "ecut": e, "tag": tag, "dir": str(d)})

    results = []

    def run_case(case, gpu):
        d = Path(case["dir"])
        out = d / case["tag"]
        res = runner.run_pw(
            d / f"{case['tag']}.in", out, pw=p.pw(kind), mpirun=p.mpirun(),
            np=1, env=p.env(), gpu=gpu, timeout=3600,
        )
        return {**case, "exit": res["exit"], "gap_eV": parse.gap_from_scf(out),
                "wall_s": res["wall_s"], "output": str(out)}

    print(f"running {len(cases)} convergence cases", flush=True)
    with ThreadPoolExecutor(max_workers=len(gpu_values)) as ex:
        futs = {ex.submit(run_case, c, gpu_values[i % len(gpu_values)]): c
                for i, c in enumerate(cases)}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            print(r["mp_id"], r["kind"], "k", r["k"], "ecut", r["ecut"],
                  "gap", r["gap_eV"], flush=True)
    results.sort(key=lambda r: (r["mp_id"], r["kind"], r["k"], r["ecut"]))
    p.results.mkdir(parents=True, exist_ok=True)
    (p.results / "convergence_results.json").write_text(json.dumps(results, indent=2))
    import pandas as pd

    pd.DataFrame(results).to_csv(p.results / "convergence_results.csv", index=False)


if __name__ == "__main__":
    main()
