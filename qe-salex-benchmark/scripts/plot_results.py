#!/usr/bin/env python3
"""Plots for sAlex force/stress/relative-energy validation."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--tag", default="salex")
    args = ap.parse_args()
    p = get_paths(args.config)
    recs = [r for r in json.loads((p.results / f"{args.tag}_frames.json").read_text())
            if r.get("status") == "ok"]
    if not recs:
        raise SystemExit("no successful frames")
    vf = np.concatenate([np.asarray(r["vasp_forces"], float).ravel() for r in recs])
    qf = np.concatenate([np.asarray(r["qe_forces"], float).ravel() for r in recs])

    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    lo, hi = min(vf.min(), qf.min()), max(vf.max(), qf.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6)
    ax.scatter(vf, qf, s=12, c="#1f77b4", alpha=0.75, edgecolor="none")
    ax.set_xlabel("sAlex force component (eV/Å)")
    ax.set_ylabel("QE force component (eV/Å)")
    ax.set_title(f"{args.tag} forces ({len(recs)} frames)\n"
                 f"MAE = {np.mean(np.abs(qf-vf)):.3f} eV/Å", fontsize=8)
    fig.tight_layout(); fig.savefig(p.results / f"{args.tag}_force_parity.png", bbox_inches="tight")
    plt.close(fig)

    vmax = [np.max(np.linalg.norm(np.asarray(r["vasp_forces"], float), axis=1)) for r in recs]
    qmax = [np.max(np.linalg.norm(np.asarray(r["qe_forces"], float), axis=1)) for r in recs]
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    hi = max(max(vmax), max(qmax)) * 1.05
    ax.plot([0, hi], [0, hi], "k--", lw=0.8, alpha=0.6)
    ax.scatter(vmax, qmax, s=30, c="#d62728", edgecolor="k", linewidth=0.3)
    ax.set_xlabel("sAlex max |F| (eV/Å)")
    ax.set_ylabel("QE max |F| (eV/Å)")
    ax.set_title("sAlex force magnitudes", fontsize=8)
    ax.set_xlim(0, hi); ax.set_ylim(0, hi); ax.set_aspect("equal", adjustable="box")
    fig.tight_layout(); fig.savefig(p.results / f"{args.tag}_force_magnitude_parity.png", bbox_inches="tight")
    plt.close(fig)

    rq, rv, colors = [], [], []
    palette = plt.cm.tab20(np.linspace(0, 1, 20))
    groups = sorted({r["mp_id"] for r in recs})
    for i, g in enumerate(groups):
        group = [r for r in recs if r["mp_id"] == g]
        qe = np.array([r["qe_energy_per_atom_eV"] for r in group])
        va = np.array([r["vasp_energy_per_atom_eV"] for r in group])
        rq.extend(((qe - qe.min()) * 1000).tolist())
        rv.extend(((va - va.min()) * 1000).tolist())
        colors.extend([palette[i]] * len(group))
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    lo, hi = min(min(rq), min(rv)), max(max(rq), max(rv))
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6)
    ax.scatter(rv, rq, s=30, c=colors, edgecolor="k", linewidth=0.3)
    ax.set_xlabel("sAlex relative energy (meV/atom)")
    ax.set_ylabel("QE relative energy (meV/atom)")
    ax.set_title(f"{args.tag} relative energies\n"
                 f"MAE = {np.mean(np.abs(np.array(rq)-np.array(rv))):.2f} meV/atom", fontsize=8)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout(); fig.savefig(p.results / f"{args.tag}_relative_energy_parity.png", bbox_inches="tight")
    plt.close(fig)

    vstr = np.concatenate([np.asarray(r["vasp_stress"], float).ravel() for r in recs
                           if r.get("vasp_stress") is not None])
    qstr = np.concatenate([np.asarray(r["qe_stress_eV_A3"], float).ravel() for r in recs
                           if r.get("qe_stress_eV_A3") is not None])
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    lo, hi = min(vstr.min(), qstr.min()), max(vstr.max(), qstr.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6)
    ax.scatter(vstr, qstr, s=14, c="#2ca02c", alpha=0.75, edgecolor="none")
    ax.set_xlabel("sAlex stress component (eV/Å³)")
    ax.set_ylabel("QE stress component (eV/Å³)")
    ax.set_title(f"{args.tag} stress ({len(vstr)} components)\n"
                 f"MAE = {np.mean(np.abs(qstr-vstr)):.2e}", fontsize=8)
    fig.tight_layout(); fig.savefig(p.results / f"{args.tag}_stress_parity.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote sAlex plots to", p.results)


if __name__ == "__main__":
    main()
