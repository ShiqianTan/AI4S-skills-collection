#!/usr/bin/env python3
"""Plots for the controlled gap and high-force convergence matrix."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--gap-tag", default="gap_controlled")
    ap.add_argument("--force-tag", default="force_matrix")
    args = ap.parse_args()
    p = get_paths(args.config)

    gap_summary = p.results / f"{args.gap_tag}_summary.json"
    if gap_summary.is_file():
        s = json.loads(gap_summary.read_text())
        per = s.get("per_pseudo", {})
        names = list(per)
        maes = [per[k]["gap_MAE_eV"] for k in names]
        fig, ax = plt.subplots(figsize=(4.8, 3.6))
        ax.bar(names, maes, color="#1f77b4")
        ax.set_ylabel("gap MAE (eV)")
        ax.set_title("Controlled band gap by pseudopotential", fontsize=8)
        for i, v in enumerate(maes):
            ax.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=7)
        fig.tight_layout(); fig.savefig(p.results / f"{args.gap_tag}_mae.png", bbox_inches="tight")
        plt.close(fig)

    force_summary = p.results / f"{args.force_tag}_summary.json"
    if force_summary.is_file():
        s = json.loads(force_summary.read_text())
        variants = s.get("variants", {})
        names = list(variants)
        maes = [variants[k]["force_MAE_eV_A"] for k in names]
        fig, ax = plt.subplots(figsize=(max(6.0, 0.45 * len(names)), 3.8))
        ax.bar(range(len(names)), maes, color="#d62728")
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=60, ha="right", fontsize=6)
        ax.set_ylabel("force MAE (eV/Å)")
        ax.set_title("High-force convergence matrix", fontsize=8)
        fig.tight_layout(); fig.savefig(p.results / f"{args.force_tag}_mae.png", bbox_inches="tight")
        plt.close(fig)
    print("wrote plots to", p.results)


if __name__ == "__main__":
    main()
