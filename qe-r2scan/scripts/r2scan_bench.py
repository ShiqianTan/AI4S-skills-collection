#!/usr/bin/env python3
"""Benchmark R2SCAN CPU configurations (MPI x OpenMP, cutoff, k-grid)."""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core import runner  # noqa: E402

CONFIGS = [
    {"tag": "e60_k9_4x4", "ecut": 60, "k": 9, "np": 4, "omp": 4},
    {"tag": "e60_k9_8x4", "ecut": 60, "k": 9, "np": 8, "omp": 4},
    {"tag": "e60_k9_4x8", "ecut": 60, "k": 9, "np": 4, "omp": 8},
    {"tag": "e60_k9_2x8", "ecut": 60, "k": 9, "np": 2, "omp": 8},
    {"tag": "e60_k9_1x16", "ecut": 60, "k": 9, "np": 1, "omp": 16},
    {"tag": "e45_k7_4x4", "ecut": 45, "k": 7, "np": 4, "omp": 4},
    {"tag": "e45_k7_8x4", "ecut": 45, "k": 7, "np": 8, "omp": 4},
]


def make_input(src: Path, dst: Path, ecut: int, k: int) -> None:
    text = src.read_text()
    text = re.sub(r"ecutwfc\s*=\s*[\d.]+", f"ecutwfc = {ecut}.0", text)
    text = re.sub(r"ecutrho\s*=\s*[\d.]+", f"ecutrho = {4*ecut}.0", text)
    text = re.sub(r"K_POINTS automatic\s*\n\s*\d+\s+\d+\s+\d+",
                  f"K_POINTS automatic\n  {k} {k} {k}", text)
    text = re.sub(r"prefix\s*=\s*'[^']*'", "prefix = 'bench'", text)
    text = re.sub(r"outdir\s*=\s*'[^']*'", "outdir = './tmp'", text)
    dst.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--base-input", required=True, help="R2SCAN scf.in to perturb")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    from ase.io import read

    p = get_paths(args.config)
    out = Path(args.out) if args.out else p.work / "r2scan_bench"
    out.mkdir(parents=True, exist_ok=True)
    results = []
    for cfg in CONFIGS:
        d = out / cfg["tag"]
        runner.prepare_tmp(d)
        make_input(Path(args.base_input), d / "scf.in", cfg["ecut"], cfg["k"])
        env = p.env()
        env["OMP_NUM_THREADS"] = str(cfg["omp"])
        run = runner.run_pw(d / "scf.in", d / "scf.out", pw=p.pw("libxc"),
                            mpirun=p.mpirun(), np=cfg["np"], env=env, timeout=3600)
        res = {"tag": cfg["tag"], **cfg, "wall_s": run["wall_s"], "ok": run["done"]}
        if run["done"]:
            atoms = read(str(d / "scf.out"), format="espresso-out")
            forces = np.asarray(atoms.get_forces(), dtype=float)
            res["energy"] = float(atoms.get_potential_energy())
            res["fmax"] = float(np.max(np.linalg.norm(forces, axis=1)))
            res["forces"] = forces.tolist()
        results.append(res)
        print(res["tag"], "wall", res.get("wall_s"), "ok", res.get("ok"),
              "E", round(res.get("energy", float("nan")), 6),
              "Fmax", round(res.get("fmax", float("nan")), 4), flush=True)

    base = next((r for r in results if r["tag"] == "e60_k9_4x4" and r.get("ok")), None)
    if base is not None:
        for r in results:
            if not r.get("ok"):
                continue
            r["dE_meV"] = round((r["energy"] - base["energy"]) * 1000.0, 3)
            r["dFmax"] = round(r["fmax"] - base["fmax"], 4)
            r["force_MAE"] = round(float(np.mean(np.abs(
                np.asarray(r["forces"]) - np.asarray(base["forces"])))), 4)
    (out / "results.json").write_text(json.dumps(
        results, indent=2, default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x)))
    print(f"wrote {out/'results.json'}")


if __name__ == "__main__":
    main()
