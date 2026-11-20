#!/usr/bin/env python3
"""qe_spectrum — broaden parsed Raman modes into a powder-Raman spectrum.

Reads ``raman_modes.csv`` in a run directory and produces a broadened
(Lorentzian/Gaussian) spectrum:

  * ``activity_spectrum`` (default): raw Raman activity (A^4/amu) vs cm^-1.
  * ``stokes_spectrum`` (optional): Stokes intensity with the
    ``(omega_L - omega)^4 / omega`` frequency factor and the Bose-Einstein
    thermal population factor ``n(omega)+1``.

Output files: ``spectrum_activity.{csv,png,svg}`` and/or
``spectrum_stokes.{csv,png,svg}``.  E-03 adds:

  * ``spectrum_ir.{csv,png,svg}`` — IR activity of the Γ-point modes.
  * ``spectrum_raman_ir.png`` — Raman and IR activity on twin axes.
  * ``spectrum_lo[_N].{csv,png,svg}`` — one spectrum per requested LO
    direction (title/label carries ``q=(x, y, z)`` from inputs_manifest.json).
  * ``spectrum_to_vs_lo[_N].png`` — the same spectrum overlaid on the Γ-point
    TO spectrum, making the LO-TO shift visible.
"""
from __future__ import annotations

import argparse
import csv
import json
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

HC_OVER_KB = 1.438776877  # h c / k_B in cm*K


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------

def load_modes(outdir: Path) -> dict[str, list[dict]]:
    """Load the mode tables from raman_modes.csv -> {table: [mode dict, ...]}."""
    tables: dict[str, list[dict]] = {}
    path = Path(outdir) / "raman_modes.csv"
    with open(path) as fh:
        for r in csv.DictReader(fh):
            t = r["table"]
            m = {
                "mode": int(r["mode"]),
                "freq": float(r["freq_cm1"]),
                "thz": float(r["freq_thz"]),
                "ir": float(r["ir"]),
                "raman": float(r["raman"]),
                "depol": float(r["depol"]) if r["depol"] != "" else None,
            }
            tables.setdefault(t, []).append(m)
    return tables


# ---------------------------------------------------------------------------
# physics helpers
# ---------------------------------------------------------------------------

def bose_plus_one(freq: float, temperature: float) -> float:
    """Bose-Einstein occupation n(omega)+1 (Stokes emission factor)."""
    if temperature <= 0.0 or freq <= 0.0:
        return 1.0
    x = HC_OVER_KB * freq / temperature
    if x > 50.0:
        return 1.0
    return 1.0 / (1.0 - math.exp(-x))


def stokes_prefactor(freq: float, laser_wavelength: float) -> float:
    """(omega_L - omega)^4 / omega, with laser wavelength in nm."""
    laser_cm1 = 1.0e7 / laser_wavelength
    w = max(freq, 1e-3)
    return (laser_cm1 - freq) ** 4 / w


def lorentz(x, x0: float, amp: float, gamma: float):
    return amp * (gamma / 2.0) ** 2 / ((x - x0) ** 2 + (gamma / 2.0) ** 2)


def gauss(x, x0: float, amp: float, fwhm: float):
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    return amp * np.exp(-((x - x0) ** 2) / (2.0 * sigma ** 2))


# ---------------------------------------------------------------------------
# spectrum construction
# ---------------------------------------------------------------------------

def validate_spectrum_params(*, fwhm: float, temperature: float,
                             laser_wavelength: float, threshold: float,
                             npoints: int, formats) -> None:
    """Raise ValueError on invalid broadening/Stokes parameters (P1-04)."""
    problems: list[str] = []
    if not (isinstance(fwhm, (int, float)) and math.isfinite(fwhm) and fwhm > 0):
        problems.append("fwhm must be > 0 (got %r)" % (fwhm,))
    if not (math.isfinite(temperature) and temperature >= 0):
        problems.append("temperature must be >= 0 K (got %r)" % (temperature,))
    if not (math.isfinite(laser_wavelength) and laser_wavelength > 0):
        problems.append("laser_wavelength must be > 0 nm (got %r)" % (laser_wavelength,))
    if not (math.isfinite(threshold) and threshold >= 0):
        problems.append("threshold must be >= 0 (got %r)" % (threshold,))
    if int(npoints) < 2:
        problems.append("npoints must be >= 2")
    bad = [f for f in formats if f not in ("png", "svg", "csv")]
    if bad:
        problems.append("unsupported formats: %s" % ", ".join(map(str, bad)))
    if problems:
        raise ValueError("; ".join(problems))


def select_peaks(modes: list[dict], threshold: float = 1e-4,
                 min_freq: float = 1.0) -> list[tuple[float, float]]:
    """Return (freq, raman_activity) for Raman-active optical modes.

    Uses the shared classifier so the same acoustic modes are excluded here as
    in the CSV/report; genuinely low-frequency optical modes are retained.
    """
    cls = classify_modes(modes)
    return [(m["freq"], m["raman"]) for m in cls["optical"]
            if m["freq"] > min_freq and m["raman"] > threshold]


def select_ir_peaks(modes: list[dict], threshold: float = 1e-4,
                    min_freq: float = 1.0) -> list[tuple[float, float]]:
    """Return (freq, IR_activity) for IR-active optical modes (E-03)."""
    cls = classify_modes(modes)
    return [(m["freq"], m["ir"]) for m in cls["optical"]
            if m["freq"] > min_freq and m["ir"] > threshold]


def lo_direction_labels(outdir, keys) -> dict[str, str]:
    """Map LO table keys ('lo', 'lo_1', ...) to 'q=(x, y, z)' labels.

    The requested directions live in inputs_manifest.json; historical imports
    without that manifest simply fall back to the bare table key.
    """
    labels: dict[str, str] = {}
    try:
        manifest = json.loads((Path(outdir) / "inputs_manifest.json").read_text())
        dirs = manifest.get("lo_dirs") or []
    except (OSError, ValueError, AttributeError):
        dirs = []
    for key in keys:
        idx = None
        if key == "lo":
            idx = 0
        elif key.startswith("lo_"):
            try:
                idx = int(key.split("_")[-1]) - 1
            except ValueError:
                idx = None
        if idx is not None and 0 <= idx < len(dirs):
            labels[key] = "q=(%s)" % ", ".join("%.3f" % float(x) for x in dirs[idx])
    return labels


def build_spectrum(peaks, *, broadening: str = "lorentzian", fwhm: float = 8.0,
                   spectrum_type: str = "activity", temperature: float = 300.0,
                   laser_wavelength: float = 532.0, npoints: int = 4000,
                   pad: float = 4.0):
    """Return (x, y, scaled_peaks); x/y are None when there are no peaks."""
    if not (isinstance(fwhm, (int, float)) and math.isfinite(fwhm) and fwhm > 0):
        raise ValueError("fwhm must be > 0 for broadening (got %r)" % (fwhm,))
    scaled: list[tuple[float, float]] = []
    for f, a in peaks:
        if spectrum_type == "stokes":
            a = a * stokes_prefactor(f, laser_wavelength) * bose_plus_one(f, temperature)
        scaled.append((f, a))
    if not scaled:
        return None, None, []
    xmin = max(0.0, min(f for f, _ in scaled) - pad * fwhm)
    xmax = max(f for f, _ in scaled) + pad * fwhm
    x = np.linspace(xmin, xmax, npoints)
    y = np.zeros_like(x)
    for f, a in scaled:
        y += gauss(x, f, a, fwhm) if broadening == "gaussian" else lorentz(x, f, a, fwhm)
    return x, y, scaled


def plot_spectrum(x, y, peaks, *, title: str, xlabel: str, ylabel: str,
                  outpath: Path, dpi: int = 160,
                  blank_text: str = "无一级拉曼活性") -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    if x is None:
        ax.text(0.5, 0.5, blank_text, ha="center", va="center",
                transform=ax.transAxes, fontsize=14)
    else:
        for f, a in peaks:
            ax.vlines(f, 0, a, color="gray", lw=1, alpha=0.5)
        ax.plot(x, y, lw=1.4)
        fmax, amax = max(peaks, key=lambda t: t[1])
        ax.annotate("%.1f" % fmax, (fmax, amax), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=9)
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def plot_joint_spectrum(x_raman, y_raman, x_ir, y_ir, *,
                        outpath: Path, dpi: int = 160) -> None:
    """Raman (left axis) and IR (right axis) powder spectra in one figure."""
    fig, ax = plt.subplots(figsize=(8, 5))
    lines = []
    if x_raman is not None:
        line, = ax.plot(x_raman, y_raman, color="tab:blue", lw=1.4,
                        label="Raman 活性")
        lines.append(line)
        ax.set_ylabel("Raman 活性 (Å$^4$/amu)", color="tab:blue")
        ax.tick_params(axis="y", labelcolor="tab:blue")
    else:
        ax.text(0.5, 0.5, "无一级拉曼活性", ha="center", va="center",
                transform=ax.transAxes, fontsize=14)

    ax2 = ax.twinx()
    if x_ir is not None:
        line, = ax2.plot(x_ir, y_ir, color="tab:red", lw=1.2, label="IR 活性")
        lines.append(line)
        ax2.set_ylabel("IR 活性 (D/Å)$^2$/amu", color="tab:red")
        ax2.tick_params(axis="y", labelcolor="tab:red")

    ax.set_title("粉末 Raman 与 IR 活性对比", fontsize=13)
    ax.set_xlabel("波数 (cm$^{-1}$)")
    if lines:
        ax.legend(lines, [ln.get_label() for ln in lines], loc="upper right",
                  fontsize=9)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def plot_to_lo_comparison(x_to, y_to, x_lo, y_lo, *, label: str = "",
                          outpath: Path, dpi: int = 160) -> None:
    """Overlay the Γ-point TO spectrum with one LO-direction spectrum."""
    fig, ax = plt.subplots(figsize=(8, 5))
    if x_to is not None:
        ax.plot(x_to, y_to, color="gray", lw=1.0, ls="--", label="Γ 点 TO")
    if x_lo is not None:
        ax.plot(x_lo, y_lo, color="tab:blue", lw=1.4,
                label=("LO " + label).strip())
    if x_to is None and x_lo is None:
        ax.text(0.5, 0.5, "无一级拉曼活性", ha="center", va="center",
                transform=ax.transAxes, fontsize=14)
    ax.set_title("TO 与 LO 方向谱对比" + (("  " + label) if label else ""),
                 fontsize=13)
    ax.set_xlabel("波数 (cm$^{-1}$)")
    ax.set_ylabel("Raman 活性 (Å$^4$/amu)")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def _write_series(base: Path, x, y, peaks, *, title: str, ylabel: str,
                  formats, dpi: int, blank_text: str =
                  "无一级拉曼活性") -> list[Path]:
    """Write one broadened series as CSV/PNG/SVG under ``base``."""
    written: list[Path] = []
    if "csv" in formats:
        csv_path = base.with_suffix(".csv")
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["freq_cm1", "intensity"])
            if x is not None:
                for xi, yi in zip(x, y):
                    w.writerow([xi, yi])
        written.append(csv_path)
    for ext in ("png", "svg"):
        if ext in formats:
            path = base.with_suffix("." + ext)
            plot_spectrum(x, y, peaks, title=title, xlabel="波数 (cm$^{-1}$)",
                          ylabel=ylabel, outpath=path, dpi=dpi,
                          blank_text=blank_text)
            written.append(path)
    return written


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def generate_spectra(outdir, *, spectrum_type: str = "activity",
                     broadening: str = "lorentzian", fwhm: float = 8.0,
                     temperature: float = 300.0, laser_wavelength: float = 532.0,
                     threshold: float = 1e-4, formats=("png", "csv"),
                     dpi: int = 160, include_ir: bool = True,
                     include_lo: bool = True) -> list[Path]:
    outdir = Path(outdir)
    validate_spectrum_params(fwhm=fwhm, temperature=temperature,
                             laser_wavelength=laser_wavelength, threshold=threshold,
                             npoints=4000, formats=formats)
    tables = load_modes(outdir)
    to_modes = tables.get("to", [])
    peaks = select_peaks(to_modes, threshold=threshold)

    types = ["activity", "stokes"] if spectrum_type == "both" else [spectrum_type]
    written: list[Path] = []
    for t in types:
        use_peaks = peaks
        if t == "stokes":
            laser_cm1 = 1.0e7 / laser_wavelength
            kept = [p for p in peaks if 0.0 < p[0] < laser_cm1]
            dropped = len(peaks) - len(kept)
            if dropped:
                print("spectrum: WARN excluded %d peak(s) outside Stokes window "
                      "0 < shift < %.1f cm-1" % (dropped, laser_cm1), file=sys.stderr)
            use_peaks = kept
        x, y, scaled = build_spectrum(
            use_peaks, broadening=broadening, fwhm=fwhm, spectrum_type=t,
            temperature=temperature, laser_wavelength=laser_wavelength,
        )
        is_stokes = t == "stokes"
        written.extend(_write_series(
            outdir / ("spectrum_" + t), x, y, scaled,
            title=("粉末拉曼光谱（Stokes）" if is_stokes else "粉末拉曼光谱（活性）"),
            ylabel=("强度 (a.u.)" if is_stokes else "Raman 活性 (Å$^4$/amu)"),
            formats=formats, dpi=dpi,
        ))

    if include_ir:
        ir_peaks = select_ir_peaks(to_modes, threshold=threshold)
        x_ir, y_ir, s_ir = build_spectrum(
            ir_peaks, broadening=broadening, fwhm=fwhm,
            spectrum_type="activity",
        )
        written.extend(_write_series(
            outdir / "spectrum_ir", x_ir, y_ir, s_ir,
            title="粉末红外光谱（IR 活性）",
            ylabel="IR 活性 (D/Å)$^2$/amu",
            formats=formats, dpi=dpi, blank_text="无 IR 活性",
        ))
        if "png" in formats:
            x_act, y_act, _ = build_spectrum(
                peaks, broadening=broadening, fwhm=fwhm,
                spectrum_type="activity",
            )
            joint = outdir / "spectrum_raman_ir.png"
            plot_joint_spectrum(x_act, y_act, x_ir, y_ir,
                                outpath=joint, dpi=dpi)
            written.append(joint)

    if include_lo:
        lo_keys = sorted(k for k in tables if k != "to")
        labels = lo_direction_labels(outdir, lo_keys)
        x_to, y_to, _ = build_spectrum(
            peaks, broadening=broadening, fwhm=fwhm, spectrum_type="activity",
        )
        for key in lo_keys:
            lo_peaks = select_peaks(tables[key], threshold=threshold)
            x, y, scaled = build_spectrum(
                lo_peaks, broadening=broadening, fwhm=fwhm,
                spectrum_type="activity",
            )
            suffix = ("  " + labels[key]) if key in labels else ""
            written.extend(_write_series(
                outdir / ("spectrum_" + key), x, y, scaled,
                title="LO 方向谱%s" % suffix,
                ylabel="Raman 活性 (Å$^4$/amu)",
                formats=formats, dpi=dpi,
            ))
            if "png" in formats:
                cmp_path = outdir / ("spectrum_to_vs_%s.png" % key)
                plot_to_lo_comparison(
                    x_to, y_to, x, y,
                    label=labels.get(key, ""), outpath=cmp_path, dpi=dpi,
                )
                written.append(cmp_path)
    return written


def _print_peaks(peaks: list[tuple[float, float]]) -> None:
    if not peaks:
        print("spectrum: no Raman-active optical modes (blank spectrum)")
        return
    for f, a in sorted(peaks, key=lambda t: t[0]):
        print("  peak %8.2f cm-1  activity %.6f" % (f, a))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="qe_spectrum.py",
                                description="broaden Raman modes into a spectrum")
    p.add_argument("--outdir", required=True, help="run directory with raman_modes.csv")
    p.add_argument("--type", default="activity",
                   choices=["activity", "stokes", "both"],
                   help="spectrum type (default activity)")
    p.add_argument("--broadening", default="lorentzian",
                   choices=["lorentzian", "gaussian"])
    p.add_argument("--fwhm", type=float, default=8.0, help="broadening FWHM (cm-1)")
    p.add_argument("--temperature", type=float, default=300.0, help="Kelvin (Stokes)")
    p.add_argument("--laser-wavelength", type=float, default=532.0,
                   help="laser wavelength nm (Stokes)")
    p.add_argument("--threshold", type=float, default=1e-4,
                   help="Raman activity cutoff")
    p.add_argument("--formats", default="png,csv",
                   help="comma-separated: png,svg,csv (default png,csv)")
    p.add_argument("--no-ir", dest="ir", action="store_false", default=True,
                   help="skip the IR activity spectrum and Raman/IR joint plot")
    p.add_argument("--no-lo-spectra", dest="lo_spectra", action="store_false",
                   default=True, help="skip the per-LO-direction spectra")
    p.add_argument("--dpi", type=int, default=160)
    args = p.parse_args(argv)

    outdir = Path(args.outdir)
    if not (outdir / "raman_modes.csv").exists():
        print("spectrum: FAIL (no raman_modes.csv in %s; run `parse` first)"
              % outdir, file=__import__("sys").stderr)
        return 1

    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    peaks = select_peaks(load_modes(outdir).get("to", []), threshold=args.threshold)
    _print_peaks(peaks)

    written = generate_spectra(
        outdir, spectrum_type=args.type, broadening=args.broadening,
        fwhm=args.fwhm, temperature=args.temperature,
        laser_wavelength=args.laser_wavelength, threshold=args.threshold,
        formats=formats, dpi=args.dpi, include_ir=args.ir,
        include_lo=args.lo_spectra,
    )
    for w in written:
        print("  wrote %s" % w)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
