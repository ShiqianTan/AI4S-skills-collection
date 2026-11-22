#!/usr/bin/env python3
"""Plots for the WBM structure/formation-energy and basin-recovery results."""
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
    ap.add_argument("--tag", default="wbm")
    args = ap.parse_args()
    p = get_paths(args.config)
    frames_path = p.results / f"{args.tag}_frames.json"
    if not frames_path.is_file():
        raise SystemExit(f"missing {frames_path}")
    recs = [r for r in json.loads(frames_path.read_text()) if r.get("status") == "ok"]

    vv = np.array([r["vasp_final_volume_A3"] for r in recs], dtype=float)
    qv = np.array([r["qe_final_volume_A3"] for r in recs], dtype=float)
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    lo, hi = min(vv.min(), qv.min()), max(vv.max(), qv.max())
    pad = 0.05 * (hi - lo or 1.0)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.8)
    ax.scatter(vv, qv, s=32, edgecolor="k", linewidth=0.3)
    ax.set_xlabel("VASP final volume (Å³)")
    ax.set_ylabel("QE relaxed volume (Å³)")
    ax.set_title(f"{args.tag} volume\nMAE = {np.mean(np.abs(qv-vv)):.3f} Å³", fontsize=8)
    fig.tight_layout(); fig.savefig(p.results / f"{args.tag}_volume_parity.png", bbox_inches="tight")
    plt.close(fig)

    form_csv = p.results / "wbm_formation_energies.csv"
    if form_csv.is_file():
        import pandas as pd

        df = pd.read_csv(form_csv)
        fig, ax = plt.subplots(figsize=(4.4, 4.0))
        x = df.wbm_Ef_MP2020_eV_atom.to_numpy()
        y = df.qe_Ef_MP2020_eV_atom.to_numpy()
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        pad = 0.05 * (hi - lo or 1.0)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.8)
        ax.scatter(x, y, s=32, edgecolor="k", linewidth=0.3)
        ax.set_xlabel("WBM Ef MP2020 (eV/atom)")
        ax.set_ylabel("QE Ef MP2020 (eV/atom)")
        ax.set_title(f"{args.tag} formation energy\nMAE = {np.mean(np.abs(y-x)):.4f} eV/atom",
                     fontsize=8)
        fig.tight_layout(); fig.savefig(p.results / f"{args.tag}_formation_energy.png",
                                        bbox_inches="tight")
        plt.close(fig)

    basin_path = p.results / "wbm_basin_summary.json"
    if basin_path.is_file():
        s = json.loads(basin_path.read_text())
        labels = ["vaspfinal", "initial", "perturb", "any_start"]
        keys = ["vaspfinal_recovery_rate", "initial_recovery_rate",
                "perturb_recovery_rate", "any_start_matches_rate"]
        vals = [s.get(k) or 0.0 for k in keys]
        fig, ax = plt.subplots(figsize=(4.6, 3.6))
        ax.bar(labels, [100 * v for v in vals], color="#1f77b4")
        ax.set_ylabel("recovery rate (%)")
        ax.set_ylim(0, 105)
        ax.set_title("WBM basin recovery", fontsize=8)
        for i, v in enumerate(vals):
            ax.text(i, 100 * v + 2, f"{100*v:.1f}%", ha="center", fontsize=7)
        fig.tight_layout(); fig.savefig(p.results / "wbm_basin_recovery.png", bbox_inches="tight")
        plt.close(fig)
    print("wrote WBM plots to", p.results)


if __name__ == "__main__":
    main()
