#!/usr/bin/env python3
"""Force, stress and relative-energy parity plots for MPtrj frames."""
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
    ap.add_argument("--tag", default="mptraj20")
    args = ap.parse_args()
    p = get_paths(args.config)
    records = [r for r in json.loads((p.results / f"{args.tag}_frames.json").read_text())
               if r.get("status") == "ok"]
    if not records:
        raise SystemExit("no successful frames")

    xs, ys, cs = [], [], []
    for r in records:
        q = np.asarray(r["qe_forces"], dtype=float).ravel()
        v = np.asarray(r["vasp_forces"], dtype=float).ravel()
        xs.extend(v.tolist())
        ys.extend(q.tolist())
        cs.extend([r["vasp_max_force_eV_A"]] * len(v))
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    lo, hi = min(min(xs), min(ys)), max(max(xs), max(ys))
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6)
    sc = ax.scatter(xs, ys, c=cs, s=10, cmap="viridis", vmin=0, vmax=3, alpha=0.75,
                    edgecolor="none")
    ax.set_xlabel("VASP force component (eV/Å)")
    ax.set_ylabel("QE force component (eV/Å)")
    ax.set_title(f"{args.tag} forces ({len(records)} frames)\n"
                 f"MAE = {np.mean(np.abs(np.array(ys) - np.array(xs))):.3f} eV/Å", fontsize=8)
    fig.colorbar(sc, ax=ax, label="VASP max |F| (eV/Å)")
    fig.tight_layout()
    fig.savefig(p.results / f"{args.tag}_force_parity.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    for r in records:
        ax.scatter(r["vasp_max_force_eV_A"], r["force_rmse_eV_A"], s=28,
                   c="#d62728", edgecolor="k", linewidth=0.3)
    ax.axhline(0.1, color="gray", ls=":", lw=0.8, label="0.1 eV/Å")
    ax.set_xlabel("VASP max |F| (eV/Å)")
    ax.set_ylabel("QE vs VASP force RMSE (eV/Å)")
    ax.set_title("Force error grows with distortion", fontsize=8)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(p.results / f"{args.tag}_force_error_vs_Fmax.png", bbox_inches="tight")
    plt.close(fig)

    rel_qe, rel_va, colors = [], [], []
    palette = plt.cm.tab20(np.linspace(0, 1, 20))
    for i, mp_id in enumerate(sorted({r["mp_id"] for r in records})):
        group = [r for r in records if r["mp_id"] == mp_id]
        qe = np.array([r["qe_energy_per_atom_eV"] for r in group])
        va = np.array([r["vasp_energy_per_atom_eV"] for r in group])
        rel_qe.extend(((qe - qe.min()) * 1000).tolist())
        rel_va.extend(((va - va.min()) * 1000).tolist())
        colors.extend([palette[i]] * len(group))
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    lo, hi = min(min(rel_qe), min(rel_va)), max(max(rel_qe), max(rel_va))
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6)
    ax.scatter(rel_va, rel_qe, s=28, c=colors, edgecolor="k", linewidth=0.3)
    ax.set_xlabel("VASP relative energy (meV/atom)")
    ax.set_ylabel("QE relative energy (meV/atom)")
    ax.set_title(f"{args.tag} relative energies\n"
                 f"MAE = {np.mean(np.abs(np.array(rel_qe) - np.array(rel_va))):.2f} meV/atom",
                 fontsize=8)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(p.results / f"{args.tag}_relative_energy_parity.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote plots to", p.results)


if __name__ == "__main__":
    main()
