#!/usr/bin/env python3
"""Plots for the MatPES-PBE vs R2SCAN comparison."""
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
from qe_bench_core.parse import KBAR_TO_EV_A3, voigt_to_matrix  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--tag", default="matpes_r2scan")
    ap.add_argument("--pbe-summary", default=None)
    args = ap.parse_args()
    p = get_paths(args.config)
    recs = [r for r in json.loads((p.results / f"{args.tag}_frames.json").read_text())
            if r.get("status") == "ok"]
    if not recs:
        raise SystemExit("no successful frames")

    vf = np.concatenate([np.asarray(r["vasp_forces"], float).ravel() for r in recs])
    qf = np.concatenate([np.asarray(r["qe_forces"], float).ravel() for r in recs])
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    lo, hi = min(vf.min(), qf.min()), max(vf.max(), qf.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6)
    ax.scatter(vf, qf, s=10, c="#9467bd", alpha=0.75, edgecolor="none")
    ax.set_xlabel("MatPES-R2SCAN force (eV/Å)")
    ax.set_ylabel("QE R2SCAN force (eV/Å)")
    ax.set_title(f"{args.tag} forces ({len(recs)} frames)\n"
                 f"MAE = {np.mean(np.abs(qf-vf)):.3f} eV/Å", fontsize=8)
    fig.tight_layout(); fig.savefig(p.results / f"{args.tag}_force_parity.png", bbox_inches="tight")
    plt.close(fig)

    vstr, qstr = [], []
    for r in recs:
        if r.get("vasp_stress") is None or r.get("qe_stress_eV_A3") is None:
            continue
        v = voigt_to_matrix(np.asarray(r["vasp_stress"], float) * KBAR_TO_EV_A3)
        q = -voigt_to_matrix(np.asarray(r["qe_stress_eV_A3"], float))
        vstr.extend(v.ravel().tolist())
        qstr.extend(q.ravel().tolist())
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    lo, hi = min(min(vstr), min(qstr)), max(max(vstr), max(qstr))
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6)
    ax.scatter(vstr, qstr, s=10, c="#2ca02c", alpha=0.75, edgecolor="none")
    ax.set_xlabel("MatPES-R2SCAN stress (eV/Å³)")
    ax.set_ylabel("QE R2SCAN stress (eV/Å³)")
    ax.set_title(f"{args.tag} stress ({len(vstr)} components)\n"
                 f"MAE = {np.mean(np.abs(np.array(qstr)-np.array(vstr))):.2e}", fontsize=8)
    fig.tight_layout(); fig.savefig(p.results / f"{args.tag}_stress_parity.png", bbox_inches="tight")
    plt.close(fig)

    summary = json.loads((p.results / f"{args.tag}_summary.json").read_text())
    pbe_path = Path(args.pbe_summary) if args.pbe_summary else p.results / "matpes_summary.json"
    labels = ["force MAE (eV/Å)", "stress MAE (eV/Å³)"]
    r2 = [summary.get("force_MAE_eV_A"), summary.get("stress_MAE_eV_A3")]
    pbe = [None, None]
    if pbe_path.is_file():
        s = json.loads(pbe_path.read_text())
        pbe = [s.get("force_MAE_eV_A"), s.get("stress_MAE_eV_A3")]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.2))
    for ax, label, a, b in zip(axes, labels, pbe, r2):
        vals = [v if v is not None else 0 for v in (a, b)]
        ax.bar(["PBE", "R2SCAN"], vals, color=["#1f77b4", "#9467bd"])
        ax.set_title(label, fontsize=8)
        for i, v in enumerate(vals):
            ax.text(i, v, f"{v:.4g}", ha="center", va="bottom", fontsize=7)
    fig.tight_layout(); fig.savefig(p.results / "r2scan_vs_pbe.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote R2SCAN plots to", p.results)


if __name__ == "__main__":
    main()
