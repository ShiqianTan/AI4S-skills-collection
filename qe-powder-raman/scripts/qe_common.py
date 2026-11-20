#!/usr/bin/env python3
"""
Shared constants and helpers for the qe-powder-raman skill.

Everything here is import-safe and side-effect free. Modules import this to
find the skill root, load defaults, and resolve QE binaries/pseudopotentials.
"""

from __future__ import annotations

import os
from pathlib import Path

# Root of the skill directory (two levels up from this script).
SKILL_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = SKILL_ROOT / "scripts"
TEMPLATE_DIR = SKILL_ROOT / "templates"
DATA_DIR = SKILL_ROOT / "data"
PSEUDO_DIR = SKILL_ROOT / "pseudo"
REFERENCES_DIR = SKILL_ROOT / "references"
EXAMPLES_DIR = SKILL_ROOT / "examples"
TESTS_DIR = SKILL_ROOT / "tests"

# ---------------------------------------------------------------------------
# Defaults (mirrors config.example.yaml; config file overrides these).
# ---------------------------------------------------------------------------
DEFAULT_CONFIG: dict = {
    # SCF
    "ecutwfc": 60.0,          # Ry
    "ecutrho": 240.0,         # Ry (4 x ecutwfc for NC pseudopotentials)
    "conv_thr": 1.0e-10,
    "mixing_beta": 0.5,
    "occupations": "fixed",
    "electron_maxstep": 200,
    "diagonalization": "david",
    # ph.x
    "tr2_ph": 1.0e-14,
    "alpha_mix": 0.5,
    # Raman / spectrum
    "fwhm": 8.0,              # cm^-1
    "temperature": 300.0,     # K (for optional Stokes spectrum)
    "laser_wavelength": 532.0,  # nm (for optional Stokes spectrum)
    "spectrum_broadening": "lorentzian",
    # k-point density (per reciprocal Angstrom): k_i = max(1, round(density/len_i))
    "kpt_density": 30.0,
    "kpt_shift": [0, 0, 0],
    # LO-TO (optional)
    "lo_to_directions": [],   # list of (x, y, z); empty = TO only
    # Run control
    "np": 4,
    "mpi_command": "mpirun",
    "omp_num_threads": 1,
    # Prefix / outdir defaults
    "prefix": "raman",
    "outdir": "./tmp",
    "pseudo_dir": "auto",     # "auto" -> bundled pseudo/ or --pseudo-dir
    "verbosity": "high",
}

# Default QE binary directory guess (machine-specific; override via config/env/CLI).
DEFAULT_BIN_DIR = Path(
    os.environ.get(
        "QE_BIN",
        "/Users/shiqian/Desktop/HefeiAILab/quantum_espresso/q-e/install_cpu/bin",
    )
)

QE_BINARIES = ("pw.x", "ph.x", "dynmat.x")

# XC functional strings accepted for LDA + norm-conserving Raman (lraman).
# HGH pseudopotentials use PZ (SLA-PZ-NOGX-NOGC). Anything else is rejected.
ACCEPTED_XC_LDA = ("PZ", "SLA-PZ-NOGX-NOGC", "SLA-PZ-NOGC-NOGX")

# Elements whose presence warns about possible magnetism (MVP is nonmagnetic).
# This is only a hint; a real magnetic check requires an SCF run.
MAGNETIC_HINT_ELEMENTS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Th", "U", "Np", "Pu",
}

# Minimum interatomic distance (Angstrom) below which atoms are "overlapping".
OVERLAP_MIN_DIST = 0.5

# Hexagonal-cell angle-noise threshold (degrees). CIFs from MP-20 sometimes
# carry gamma ~ 119.99998586 deg; spglib.refine_cell should be applied first.
# Any deviation from exactly 120 deg beyond ~1e-5 is treated as numerical noise.
HEX_ANGLE_NOISE_DEG = 1.0e-5

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def load_defaults() -> dict:
    """Return a deep copy of DEFAULT_CONFIG (plus any YAML defaults if present)."""
    cfg = dict(DEFAULT_CONFIG)
    cfg["kpt_shift"] = list(DEFAULT_CONFIG["kpt_shift"])
    cfg["lo_to_directions"] = list(DEFAULT_CONFIG["lo_to_directions"])
    # Optionally merge data/defaults.yaml (kept minimal; config.example.yaml is
    # the canonical override). Imported lazily to avoid a hard PyYAML dep here.
    defaults_yaml = DATA_DIR / "defaults.yaml"
    if defaults_yaml.exists():
        try:
            import yaml  # type: ignore

            extra = yaml.safe_load(defaults_yaml.read_text()) or {}
            if isinstance(extra, dict):
                for k, v in extra.items():
                    if k in cfg and isinstance(cfg[k], (list, tuple)):
                        cfg[k] = list(v)
                    else:
                        cfg[k] = v
        except Exception:  # pragma: no cover - defaults must never break import
            pass
    return cfg


def resolve_bin_dir(bin_dir: str | Path | None) -> Path:
    """Resolve the QE binary directory from CLI > env > default."""
    if bin_dir:
        return Path(bin_dir).expanduser()
    return Path(DEFAULT_BIN_DIR).expanduser()


def resolve_pseudo_dir(pseudo_dir: str | Path | None) -> Path:
    """Resolve the pseudopotential directory from CLI > bundled default."""
    if pseudo_dir and str(pseudo_dir).strip() not in ("", "auto"):
        return Path(pseudo_dir).expanduser()
    return PSEUDO_DIR


def atomic_number(symbol: str) -> int:
    """Atomic number for an element symbol (used for lightweight checks)."""
    table = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
        "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
        "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
        "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
        "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
        "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
        "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
        "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
        "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
        "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
        "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
        "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
        "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92,
        "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99,
    }
    return table[symbol.capitalize()]
