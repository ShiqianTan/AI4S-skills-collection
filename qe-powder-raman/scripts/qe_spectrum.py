#!/usr/bin/env python3
"""qe_spectrum — broaden parsed Raman modes into a powder-Raman spectrum.

Reads ``raman_modes.csv`` in a run directory and produces a broadened
(Lorentzian/Gaussian) spectrum:

  * ``activity_spectrum`` (default): raw Raman activity (A^4/amu) vs cm^-1.
  * ``stokes_spectrum`` (optional): Stokes intensity with the
    ``(omega_L - omega)^4 / omega`` frequency factor and the Bose-Einstein
    thermal population factor ``n(omega)+1``.

Output files: ``spectrum_activity.{csv,png,svg}`` and/or
``spectrum_stokes.{csv,png,svg}``.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

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

def select_peaks(modes: list[dict], threshold: float = 1e-4,
                 min_freq: float = 1.0) -> list[tuple[float, float]]:
    """Return (freq, raman_activity) for Raman-active optical modes."""
    return [(m["freq"], m["raman"]) for m in modes
            if m["freq"] > min_freq and m["raman"] > threshold]


def build_spectrum(peaks, *, broadening: str = "lorentzian", fwhm: float = 8.0,
                   spectrum_type: str = "activity", temperature: float = 300.0,
                   laser_wavelength: float = 532.0, npoints: int = 4000,
                   pad: float = 4.0):
    """Return (x, y, scaled_peaks); x/y are None when there are no peaks."""
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
                  outpath: Path, dpi: int = 160) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    if x is None:
        ax.text(0.5, 0.5, "无一级拉曼活性", ha="center", va="center",
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


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def generate_spectra(outdir, *, spectrum_type: str = "activity",
                     broadening: str = "lorentzian", fwhm: float = 8.0,
                     temperature: float = 300.0, laser_wavelength: float = 532.0,
                     threshold: float = 1e-4, formats=("png", "csv"),
                     dpi: int = 160) -> list[Path]:
    outdir = Path(outdir)
    tables = load_modes(outdir)
    to_modes = tables.get("to", [])
    peaks = select_peaks(to_modes, threshold=threshold)

    types = ["activity", "stokes"] if spectrum_type == "both" else [spectrum_type]
    written: list[Path] = []
    for t in types:
        x, y, scaled = build_spectrum(
            peaks, broadening=broadening, fwhm=fwhm, spectrum_type=t,
            temperature=temperature, laser_wavelength=laser_wavelength,
        )
        base = outdir / ("spectrum_" + t)
        is_stokes = t == "stokes"
        title = ("粉末拉曼光谱（Stokes）" if is_stokes else "粉末拉曼光谱（活性）")
        ylabel = "强度 (a.u.)" if is_stokes else "Raman 活性 (Å$^4$/amu)"

        if "csv" in formats:
            csv_path = base.with_suffix(".csv")
            with open(csv_path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["freq_cm1", "intensity"])
                if x is not None:
                    for xi, yi in zip(x, y):
                        w.writerow([xi, yi])
            written.append(csv_path)
        if "png" in formats:
            png_path = base.with_suffix(".png")
            plot_spectrum(x, y, scaled, title=title, xlabel="波数 (cm$^{-1}$)",
                          ylabel=ylabel, outpath=png_path, dpi=dpi)
            written.append(png_path)
        if "svg" in formats:
            svg_path = base.with_suffix(".svg")
            plot_spectrum(x, y, scaled, title=title, xlabel="波数 (cm$^{-1}$)",
                          ylabel=ylabel, outpath=svg_path, dpi=dpi)
            written.append(svg_path)
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
        formats=formats, dpi=args.dpi,
    )
    for w in written:
        print("  wrote %s" % w)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
