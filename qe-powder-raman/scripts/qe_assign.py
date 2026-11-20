#!/usr/bin/env python3
"""qe_assign — mode assignment and experimental comparison (E-04).

The assignment is activity based: every optical mode is labelled
``Raman+IR`` / ``Raman`` / ``IR`` / ``silent`` (or ``soft`` for imaginary
frequencies), and near-degenerate modes are grouped within
``degeneracy_tol`` cm^-1.  Irreducible representations are NOT computed
(that would require eigenvectors plus space-group character tables).

An optional experimental file can be supplied:

  * a peak table ``freq_cm1,intensity[,error]`` (few rows), or
  * a full spectrum (many rows), in which case peaks are detected and each
    peak is fitted with a Lorentzian (scipy if available, otherwise a
    moment estimate).

Outputs: ``mode_assignment.csv``, ``assignment.md`` and
``assignment_overlay.{png,svg}`` (predicted sticks + experimental data +
fitted peaks), plus ``experimental_fit.csv`` when a spectrum was fitted.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from qe_parse import classify_modes

plt.rcParams["font.sans-serif"] = [
    "Songti SC", "PingFang HK", "Hiragino Sans GB",
    "Arial Unicode MS", "STHeiti", "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False

SPECTRUM_ROWS = 50  # more rows than this -> treat the file as a spectrum
ACTIVITY_LABELS = {
    "both": "Raman+IR", "raman": "仅 Raman", "ir": "仅 IR",
    "silent": "静默", "soft": "软模",
}


# ---------------------------------------------------------------------------
# assignment
# ---------------------------------------------------------------------------

def _activity(mode: dict, threshold: float) -> str:
    if mode["freq"] < 0.0:
        return "soft"
    raman = mode.get("raman") or 0.0
    ir = mode.get("ir") or 0.0
    if raman > threshold and ir > threshold:
        return "both"
    if raman > threshold:
        return "raman"
    if ir > threshold:
        return "ir"
    return "silent"


def assign_modes(modes: list[dict], threshold: float = 1e-4,
                 degeneracy_tol: float = 1.0, table: str = "to") -> dict:
    """Label optical modes by activity and group near-degenerate partners."""
    cls = classify_modes(modes)
    optical = cls["optical"]
    entries: list[dict] = []
    used: set[int] = set()
    group = 0
    for i, mode in enumerate(optical):
        if i in used:
            continue
        partners = [i]
        for j in range(i + 1, len(optical)):
            if j in used:
                continue
            if abs(optical[j]["freq"] - mode["freq"]) <= degeneracy_tol:
                partners.append(j)
        used.update(partners)
        group += 1
        for j in partners:
            item = optical[j]
            entries.append({
                "mode": item["mode"],
                "freq": item["freq"],
                "thz": item["thz"],
                "ir": item.get("ir") or 0.0,
                "raman": item.get("raman") or 0.0,
                "depol": item.get("depol"),
                "activity": _activity(item, threshold),
                "table": table,
                "group": group,
                "degeneracy": len(partners),
                "partners": [optical[k]["mode"] for k in partners],
            })
    entries.sort(key=lambda e: (e["freq"], e["mode"]))
    counts: dict[str, int] = {}
    for e in entries:
        counts[e["activity"]] = counts.get(e["activity"], 0) + 1
    return {
        "entries": entries,
        "acoustic": [m["mode"] for m in cls["acoustic"]],
        "counts": counts,
        "threshold": threshold,
        "degeneracy_tol": degeneracy_tol,
    }


# ---------------------------------------------------------------------------
# experimental data
# ---------------------------------------------------------------------------

def read_experimental(path: str | Path) -> list[dict]:
    """Read ``freq,intensity[,error]`` rows (CSV or whitespace separated)."""
    rows: list[dict] = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",") if "," in line else line.split()
        parts = [p.strip() for p in parts if p.strip()]
        try:
            vals = [float(p) for p in parts[:3]]
        except ValueError:
            continue  # header
        if len(vals) < 2:
            continue
        rows.append({
            "freq": vals[0],
            "intensity": vals[1],
            "error": vals[2] if len(vals) > 2 else None,
        })
    return rows


def _lorentz(x, x0, amp, fwhm):
    g = 0.5 * fwhm
    return amp * g * g / ((x - x0) ** 2 + g * g)


def detect_peaks(x, y, *, min_prominence: float = 0.05,
                 max_peaks: int = 30) -> list[float]:
    """Return local maxima above ``min_prominence * max(y)``."""
    ymax = float(np.max(y)) if len(y) else 0.0
    if ymax <= 0.0:
        return []
    peaks: list[float] = []
    for i in range(1, len(y) - 1):
        if y[i] >= y[i - 1] and y[i] > y[i + 1] and y[i] >= min_prominence * ymax:
            peaks.append(float(x[i]))
    peaks.sort(key=lambda f: -float(np.interp(f, x, y)))
    return sorted(peaks[:max_peaks])


def _moment_fit(x, y, amp):
    """Fallback fit: intensity centroid + half-maximum width."""
    w = np.clip(y, 0.0, None)
    total = float(w.sum())
    dx = float(np.median(np.diff(x))) if len(x) > 1 else 1.0
    if total <= 0.0:
        return None
    center = float((x * w).sum() / total)
    half = 0.5 * amp
    above = np.where(w >= half)[0]
    if len(above) >= 2:
        fwhm = float(x[above[-1]] - x[above[0]])
    else:
        fwhm = 2.0 * dx
    fwhm = max(fwhm, dx)
    n_eff = max(total / max(amp, 1e-12), 1.0)
    return {
        "center": center,
        "amplitude": amp,
        "fwhm": fwhm,
        "center_err": fwhm / (2.355 * math.sqrt(n_eff)),
        "fwhm_err": fwhm / math.sqrt(n_eff),
        "method": "moment",
    }


def fit_single_peak(x, y, x0: float, *, window: float | None = None) -> dict | None:
    """Fit one Lorentzian around ``x0`` (scipy if available)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        return None
    dx = float(np.median(np.diff(x))) if len(x) > 1 else 1.0
    half = window if window is not None else max(20.0 * dx, 5.0)
    mask = np.abs(x - x0) <= half
    if mask.sum() < 3:
        mask = np.ones_like(x, dtype=bool)
    amp = float(np.max(y[mask]))
    if amp <= 0.0:
        return None
    try:
        from scipy.optimize import curve_fit  # optional dependency

        p0 = [x0, amp, max(2.0 * dx, 1.0)]
        bounds = (
            [float(x[mask].min()), 0.0, 1e-3],
            [float(x[mask].max()), np.inf, float(x[mask].max() - x[mask].min())],
        )
        popt, pcov = curve_fit(_lorentz, x[mask], y[mask], p0=p0, bounds=bounds,
                               maxfev=20000)
        err = np.sqrt(np.clip(np.diag(pcov), 0.0, None))
        return {
            "center": float(popt[0]),
            "amplitude": float(popt[1]),
            "fwhm": float(popt[2]),
            "center_err": float(err[0]),
            "fwhm_err": float(err[2]),
            "method": "scipy",
        }
    except Exception:  # noqa: BLE001
        return _moment_fit(x[mask], y[mask], amp)


def fit_spectrum(rows: list[dict], *, min_prominence: float = 0.05) -> list[dict]:
    """Detect and fit all peaks of an experimental spectrum."""
    x = np.array([r["freq"] for r in rows], dtype=float)
    y = np.array([r["intensity"] for r in rows], dtype=float)
    fits = []
    for x0 in detect_peaks(x, y, min_prominence=min_prominence):
        fit = fit_single_peak(x, y, x0)
        if fit:
            fits.append(fit)
    return fits


def match_experimental(entries: list[dict], exp_peaks: list[dict],
                       tol_cm1: float = 20.0) -> dict:
    """Nearest-predicted-mode matching with a frequency tolerance."""
    predicted = [e for e in entries if e["freq"] > 0.0]
    matches: list[dict] = []
    deltas: list[float] = []
    for row in exp_peaks:
        if not predicted:
            matches.append({**row, "matched": False})
            continue
        best = min(predicted, key=lambda e: abs(e["freq"] - row["freq"]))
        delta = float(row["freq"] - best["freq"])
        within = abs(delta) <= tol_cm1
        if within:
            deltas.append(delta)
        matches.append({
            **row,
            "matched": within,
            "predicted_mode": best["mode"],
            "predicted_table": best.get("table", "to"),
            "predicted_freq": best["freq"],
            "predicted_activity": best["activity"],
            "delta": delta,
        })
    return {
        "matches": matches,
        "tolerance": tol_cm1,
        "n_matched": len(deltas),
        "rms_delta": (math.sqrt(sum(d * d for d in deltas) / len(deltas))
                      if deltas else None),
    }


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------

def _fmt(value, nd: int = 4) -> str:
    if value is None:
        return "—"
    return ("%." + str(nd) + "f") % value


def write_assignment(outdir, assignment: dict, matching: dict | None = None,
                     fits: list[dict] | None = None,
                     exp_rows: list[dict] | None = None,
                     formats=("png", "svg"), dpi: int = 160) -> list[Path]:
    outdir = Path(outdir)
    entries = assignment["entries"]
    written: list[Path] = []

    csv_path = outdir / "mode_assignment.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["table", "mode", "freq_cm1", "thz", "ir", "raman", "depol",
                    "activity", "group", "degeneracy"])
        for e in entries:
            w.writerow([
                e.get("table", "to"), e["mode"], "%.6f" % e["freq"], "%.6f" % e["thz"],
                "%.6f" % e["ir"], "%.6f" % e["raman"],
                "" if e["depol"] is None else "%.6f" % e["depol"],
                e["activity"], e["group"], e["degeneracy"],
            ])
    written.append(csv_path)

    if fits:
        fit_path = outdir / "experimental_fit.csv"
        with open(fit_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["center_cm1", "amplitude", "fwhm_cm1",
                        "center_err", "fwhm_err", "method"])
            for f in fits:
                w.writerow(["%.6f" % f["center"], "%.6f" % f["amplitude"],
                            "%.6f" % f["fwhm"], "%.6f" % f["center_err"],
                            "%.6f" % f["fwhm_err"], f["method"]])
        written.append(fit_path)

    md = [
        "# 模式指认与实验对照",
        "",
        "- 光学模数：%d" % len(entries),
        "- 声学模：%s" % (", ".join(str(m) for m in assignment["acoustic"]) or "—"),
        "- 活性阈值：%g；简并容差：%g cm⁻¹"
        % (assignment["threshold"], assignment["degeneracy_tol"]),
    ]
    counts = assignment["counts"]
    md.append("- 活性分类：" + "，".join(
        "%s %d" % (ACTIVITY_LABELS.get(k, k), counts[k])
        for k in ("both", "raman", "ir", "silent", "soft") if counts.get(k)))
    md.extend([
        "",
        "| 模 | 频率 (cm⁻¹) | 活性 | Raman | IR | 退偏 | 简并 | 组 |",
        "|---|---|---|---|---|---|---|---|",
    ])
    for e in entries:
        md.append("| %s | %s | %s | %s | %s | %s | %s | %s |" % (
            e["mode"], _fmt(e["freq"], 2),
            ACTIVITY_LABELS.get(e["activity"], e["activity"]),
            _fmt(e["raman"]), _fmt(e["ir"]), _fmt(e["depol"]),
            e["degeneracy"], e["group"]))
    if matching and matching.get("matches"):
        md.extend([
            "",
            "## 实验峰匹配（容差 %s cm⁻¹）" % _fmt(matching["tolerance"], 1),
            "",
            "| 实验 (cm⁻¹) | 强度 | 表 | 匹配模 | 预测 (cm⁻¹) | 偏差 | 状态 |",
            "|---|---|---|---|---|---|---|",
        ])
        for m in matching["matches"]:
            md.append("| %s | %s | %s | %s | %s | %s | %s |" % (
                _fmt(m["freq"], 2), _fmt(m.get("intensity"), 3),
                m.get("predicted_table", "—"), m.get("predicted_mode", "—"),
                _fmt(m.get("predicted_freq"), 2),
                _fmt(m.get("delta"), 2), "匹配" if m["matched"] else "未匹配"))
        md.append("")
        md.append("- 匹配数：%d；RMS 偏差：%s cm⁻¹"
                  % (matching["n_matched"], _fmt(matching["rms_delta"], 2)))
    if fits:
        md.extend([
            "",
            "## 实验谱峰拟合",
            "",
            "| 中心 (cm⁻¹) | 振幅 | FWHM (cm⁻¹) | 中心误差 | 方法 |",
            "|---|---|---|---|---|",
        ])
        for f in fits:
            md.append("| %s | %s | %s | %s | %s |" % (
                _fmt(f["center"], 2), _fmt(f["amplitude"], 3),
                _fmt(f["fwhm"], 2), _fmt(f["center_err"], 2), f["method"]))
        md.append("")
        md.append("- 误差为 1σ 拟合不确定度；无 scipy 时使用矩估计（method=moment）。")
    md.extend([
        "",
        "> 活性指认由 Raman/IR 张量给出；未计算不可约表示。",
    ])
    (outdir / "assignment.md").write_text("\n".join(md) + "\n")
    written.append(outdir / "assignment.md")

    written.extend(_plot_overlay(outdir, entries, fits, exp_rows,
                                 matching, formats=formats, dpi=dpi))
    return written


def _plot_overlay(outdir, entries, fits, exp_rows, matching, *,
                  formats, dpi):
    pred = [e for e in entries if e["freq"] > 0.0 and e["raman"] > 0.0]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    if exp_rows:
        x = np.array([r["freq"] for r in exp_rows], dtype=float)
        y = np.array([r["intensity"] for r in exp_rows], dtype=float)
        if len(exp_rows) > SPECTRUM_ROWS:
            ax.plot(x, y, color="0.4", lw=1.0, label="实验谱")
        else:
            err = [r["error"] for r in exp_rows]
            if any(e is not None for e in err):
                ax.errorbar(x, y, yerr=[e if e is not None else 0.0 for e in err],
                            fmt="o", ms=4, color="0.3", capsize=3,
                            label="实验峰")
            else:
                ax.plot(x, y, "o", ms=4, color="0.3", label="实验峰")
    if pred:
        amax = max(e["raman"] for e in pred) or 1.0
        expmax = max((r["intensity"] for r in exp_rows), default=1.0) if exp_rows else 1.0
        for e in pred:
            ax.vlines(e["freq"], 0.0, 0.9 * expmax * e["raman"] / amax,
                      color="tab:blue", lw=1.2, alpha=0.8)
        ax.vlines([], [], [], color="tab:blue", lw=1.2, label="预测 Raman 模")
    if fits:
        if exp_rows and len(exp_rows) > SPECTRUM_ROWS:
            xg = np.linspace(min(r["freq"] for r in exp_rows),
                             max(r["freq"] for r in exp_rows), 1200)
        else:
            xg = np.linspace(min(f["center"] for f in fits) - 50,
                             max(f["center"] for f in fits) + 50, 1200)
        for f in fits:
            ax.plot(xg, _lorentz(xg, f["center"], f["amplitude"], f["fwhm"]),
                    color="tab:red", lw=1.0, ls="--", alpha=0.8)
        ax.plot([], [], color="tab:red", lw=1.0, ls="--", label="拟合峰")
    ax.set_title("模式指认与实验对照", fontsize=13)
    ax.set_xlabel("波数 (cm$^{-1}$)")
    ax.set_ylabel("强度 / Raman 活性 (a.u.)")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    written = []
    for ext in ("png", "svg"):
        if ext in formats:
            path = Path(outdir) / ("assignment_overlay." + ext)
            fig.savefig(path, dpi=dpi)
            written.append(path)
    plt.close(fig)
    return written


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="qe_assign.py",
                                description="mode assignment + experimental matching")
    p.add_argument("--outdir", required=True)
    p.add_argument("--experimental", default=None,
                   help="CSV/whitespace freq,intensity[,error]; >50 rows = spectrum")
    p.add_argument("--tol", type=float, default=20.0,
                   help="experimental matching tolerance (cm^-1)")
    p.add_argument("--degeneracy-tol", type=float, default=1.0)
    p.add_argument("--threshold", type=float, default=1e-4)
    p.add_argument("--min-prominence", type=float, default=0.05)
    p.add_argument("--formats", default="png,svg")
    p.add_argument("--dpi", type=int, default=160)
    args = p.parse_args(argv)

    from qe_spectrum import load_modes

    outdir = Path(args.outdir)
    tables = load_modes(outdir)
    modes = tables.get("to", [])
    assignment = assign_modes(modes, threshold=args.threshold,
                              degeneracy_tol=args.degeneracy_tol)
    match_entries = list(assignment["entries"])
    for key in sorted(k for k in tables if k != "to"):
        match_entries.extend(assign_modes(
            tables[key], threshold=args.threshold,
            degeneracy_tol=args.degeneracy_tol, table=key)["entries"])
    matching = None
    fits = None
    exp_rows = None
    if args.experimental:
        exp_rows = read_experimental(args.experimental)
        if len(exp_rows) > SPECTRUM_ROWS:
            fits = fit_spectrum(exp_rows, min_prominence=args.min_prominence)
            exp_peaks = [{"freq": f["center"], "intensity": f["amplitude"],
                          "error": f["center_err"]} for f in fits]
        else:
            exp_peaks = exp_rows
        if exp_peaks:
            matching = match_experimental(match_entries, exp_peaks,
                                          tol_cm1=args.tol)
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    written = write_assignment(outdir, assignment, matching, fits, exp_rows,
                               formats=formats, dpi=args.dpi)
    for e in assignment["entries"]:
        print("  mode %-3d %8.2f cm-1  %-9s deg=%d"
              % (e["mode"], e["freq"],
                 ACTIVITY_LABELS.get(e["activity"], e["activity"]),
                 e["degeneracy"]))
    if matching:
        print("assign: %d/%d experimental peaks matched (RMS %s cm-1)"
              % (matching["n_matched"], len(matching["matches"]),
                 _fmt(matching["rms_delta"], 2)))
    print("assign: wrote " + ", ".join(p.name for p in written))
    return 0


if __name__ == "__main__":
    sys.exit(main())
