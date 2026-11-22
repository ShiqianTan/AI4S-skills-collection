"""Shared parity plots for the benchmark skills."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def parity_plot(
    reference,
    predicted,
    outfile: str | Path,
    xlabel: str = "Reference",
    ylabel: str = "QE",
    title: str | None = None,
    unit: str = "",
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ref = np.asarray(reference, dtype=float)
    pred = np.asarray(predicted, dtype=float)
    diff = np.abs(pred - ref)
    mae = float(np.mean(diff)) if diff.size else 0.0
    lo = float(min(ref.min(), pred.min()))
    hi = float(max(ref.max(), pred.max()))
    pad = 0.05 * (hi - lo or 1.0)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(ref, pred, s=28, alpha=0.8)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.8)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel(f"{xlabel} ({unit})" if unit else xlabel)
    ax.set_ylabel(f"{ylabel} ({unit})" if unit else ylabel)
    if title:
        ax.set_title(title)
    ax.text(
        0.04,
        0.96,
        f"MAE = {mae:.4g} {unit}".strip(),
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    fig.tight_layout()
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=160)
    plt.close(fig)
    return outfile
