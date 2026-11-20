#!/usr/bin/env python3
"""
Structure reading and preflight checks for the qe-powder-raman skill.

Step 0 scope:
  - read a structure (CIF / POSCAR / PW input) via ASE
  - detect partial occupancies and overlapping atoms
  - run spglib symmetry detection
  - flag near-hexagonal cells with numerical angle noise (recommend refine_cell)

Step 1 will extend this with spglib.refine_cell, validated.cif output, and a
pseudopotential manifest.
"""

from __future__ import annotations

import sys
from pathlib import Path

from qe_common import (
    HEX_ANGLE_NOISE_DEG,
    MAGNETIC_HINT_ELEMENTS,
    OVERLAP_MIN_DIST,
)


class StructureCheckResult:
    """Container for the outcome of structure checks."""

    def __init__(self) -> None:
        self.messages: list[tuple[str, str, str]] = []  # (level, key, text)
        self.n_atoms = 0
        self.elements: list[str] = []
        self.spacegroup_number: int | None = None
        self.spacegroup_symbol: str | None = None
        self.n_sym_ops = 0
        self.cell = None
        self.angles: tuple[float, float, float] | None = None

    def add(self, level: str, key: str, text: str) -> None:
        self.messages.append((level, key, text))

    @property
    def has_fail(self) -> bool:
        return any(lvl == "FAIL" for lvl, _, _ in self.messages)

    @property
    def has_warn(self) -> bool:
        return any(lvl == "WARN" for lvl, _, _ in self.messages)


def read_structure(path: str | Path):
    """Read a structure with ASE. Returns an ase.Atoms or raises."""
    from ase import io

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"structure file not found: {p}")
    fmt = None
    if p.suffix.lower() == ".cif":
        fmt = "cif"
    elif p.suffix.lower() in (".vasp", ".poscar", ".contcar"):
        fmt = "vasp"
    elif p.name.lower().startswith(("poscar", "contcar")):
        fmt = "vasp"
    atoms = io.read(str(p), format=fmt)
    return atoms


def check_structure(path: str | Path) -> StructureCheckResult:
    """Run structure-level preflight checks. Never raises on checkable problems."""
    res = StructureCheckResult()
    try:
        atoms = read_structure(path)
    except Exception as exc:  # noqa: BLE001
        res.add("FAIL", "read", f"cannot read structure: {exc}")
        return res

    res.n_atoms = len(atoms)
    res.elements = list(atoms.get_chemical_symbols())
    res.cell = atoms.get_cell()
    try:
        cellpar = atoms.cell.cellpar()
        res.angles = tuple(round(float(a), 6) for a in cellpar[3:])
    except Exception:  # noqa: BLE001
        res.angles = None

    # ---- partial occupancy -------------------------------------------------
    _check_partial_occupancy(atoms, res)

    # ---- overlap -----------------------------------------------------------
    _check_overlap(atoms, res)

    # ---- symmetry ----------------------------------------------------------
    _check_symmetry(atoms, res)

    # ---- magnetic hint -----------------------------------------------------
    mag = [e for e in set(res.elements) if e in MAGNETIC_HINT_ELEMENTS]
    if mag:
        res.add(
            "WARN",
            "magnetic_hint",
            "elements %s may be magnetic; MVP assumes nonmagnetic "
            "(verify with a spin-polarized SCF if unsure)" % ", ".join(sorted(mag)),
        )

    if res.n_atoms == 0:
        res.add("FAIL", "empty", "structure has 0 atoms")
    return res


def _check_partial_occupancy(atoms, res: StructureCheckResult) -> None:
    # ASE CIF reader stores occupancies in atoms.info["occupancy"] (per-site).
    occ = atoms.info.get("occupancy")
    if occ is None:
        res.add("PASS", "occupancy", "no partial occupancies detected")
        return
    try:
        frac = [float(o) for o in occ if float(o) > 0 and float(o) < 1.0 - 1e-8]
    except Exception:  # noqa: BLE001
        frac = []
    if frac:
        res.add(
            "FAIL",
            "occupancy",
            "%d site(s) have partial occupancy (site occupancy must be 1.0)" % len(frac),
        )
    else:
        res.add("PASS", "occupancy", "all sites fully occupied")


def _check_overlap(atoms, res: StructureCheckResult) -> None:
    import numpy as np

    n = len(atoms)
    if n < 2:
        res.add("PASS", "overlap", "fewer than 2 atoms; no overlap check needed")
        return
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    pbc = atoms.get_pbc().any()
    from ase.geometry import get_distances  # type: ignore

    _, dmat = get_distances(pos, cell=cell, pbc=atoms.get_pbc())
    np.fill_diagonal(dmat, np.inf)
    if dmat.size and float(dmat.min()) < OVERLAP_MIN_DIST:
        i, j = np.unravel_index(np.argmin(dmat), dmat.shape)
        res.add(
            "FAIL",
            "overlap",
            "atoms %d and %d are %.3f A apart (< %.2f A): overlapping"
            % (i, j, float(dmat.min()), OVERLAP_MIN_DIST),
        )
    else:
        res.add("PASS", "overlap", "no overlapping atoms")


def _check_symmetry(atoms, res: StructureCheckResult) -> None:
    import spglib

    try:
        cell = (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
        ds = spglib.get_symmetry_dataset(cell, symprec=1e-3)
        if ds is None:
            res.add("WARN", "symmetry", "spglib could not determine symmetry")
            return
        res.spacegroup_number = ds.number
        res.spacegroup_symbol = ds.international
        res.n_sym_ops = len(ds.rotations)
        res.add(
            "PASS",
            "symmetry",
            "space group #%d %s (%d sym ops)"
            % (ds.number, ds.international, len(ds.rotations)),
        )
    except Exception as exc:  # noqa: BLE001
        res.add("WARN", "symmetry", "spglib symmetry detection failed: %s" % exc)
        return

    # ---- hexagonal-angle noise --------------------------------------------
    if res.angles is not None:
        a, b, g = res.angles
        # hexagonal family: a==b, alpha==beta==90, gamma~120
        lengths = atoms.cell.cellpar()[:3]
        if (
            abs(lengths[0] - lengths[1]) < 1e-3
            and abs(a - 90.0) < 1e-2
            and abs(b - 90.0) < 1e-2
            and abs(g - 120.0) > HEX_ANGLE_NOISE_DEG
            and abs(g - 120.0) < 1.0
        ):
            res.add(
                "WARN",
                "hex_angle_noise",
                "gamma=%.8f deg (hexagonal expects 120.0); numerical angle noise "
                "may break ph.x symmetry analysis -> run spglib.refine_cell first" % g,
            )
        # rhombohedral family: a==b==c, alpha==beta==gamma ~ 60
        if (
            abs(lengths[0] - lengths[1]) < 1e-3
            and abs(lengths[1] - lengths[2]) < 1e-3
            and abs(a - 60.0) < 1.0
            and abs(b - 60.0) < 1.0
            and abs(g - 60.0) < 1.0
        ):
            res.add(
                "INFO",
                "rhombohedral",
                "primitive rhombohedral cell detected (angles ~60 deg); QE may "
                "only resolve a subgroup of the full space group",
            )


def format_check_result(res: StructureCheckResult) -> str:
    lines = ["structure: %d atoms, elements=%s" % (res.n_atoms, ",".join(res.elements))]
    if res.angles is not None:
        lines.append("cell angles (deg): %s" % ", ".join("%.6f" % a for a in res.angles))
    for level, key, text in res.messages:
        lines.append("[%s] %-16s %s" % (level, key, text))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Step 1: symmetry refinement (spglib.refine_cell) + validated.cif output
# ---------------------------------------------------------------------------


def refine_structure(atoms, symprec: float = 1e-3):
    """Idealize a structure via spglib.refine_cell.

    Fixes small numerical angle/length noise (e.g. hexagonal gamma=119.99998586
    deg) WITHOUT changing the atom count or forcing a different space group.
    It is a "cleanup", not a symmetry search: a genuinely P1-degenerate CIF
    (e.g. MP-20 Al2O3) stays low-symmetry.

    Returns (refined_atoms, info). `info` records original vs refined space
    group, cell parameters, and whether the cell actually changed.
    """
    import numpy as np  # noqa: F401  (kept for clarity of fallback checks)
    import spglib
    from ase import Atoms

    cell = (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())

    def _ds(c):
        d = spglib.get_symmetry_dataset(c, symprec=symprec)
        return d

    ds0 = _ds(cell)
    original = {
        "spacegroup_number": ds0.number if ds0 else None,
        "spacegroup_symbol": ds0.international if ds0 else None,
        "n_sym_ops": len(ds0.rotations) if ds0 else 0,
    }

    refined = atoms
    refined_info = original
    changed = False
    cp0 = list(atoms.cell.cellpar())

    try:
        lat, pos, nums = spglib.refine_cell(cell, symprec=symprec)
        if lat is not None and len(nums) == len(atoms):
            refined = Atoms(numbers=nums, scaled_positions=pos, cell=lat, pbc=True)
            ds1 = _ds((lat, pos, nums))
            refined_info = {
                "spacegroup_number": ds1.number if ds1 else None,
                "spacegroup_symbol": ds1.international if ds1 else None,
                "n_sym_ops": len(ds1.rotations) if ds1 else 0,
            }
            cp1 = list(refined.cell.cellpar())
            changed = (
                max(abs(a - b) for a, b in zip(cp0[3:], cp1[3:])) > 1e-6
                or max(abs(a - b) for a, b in zip(cp0[:3], cp1[:3])) > 1e-6
            )
            cp0 = cp1
    except Exception as exc:  # noqa: BLE001
        refined_info = dict(original)
        refined_info["refine_error"] = str(exc)

    info = {
        "original": original,
        "refined": refined_info,
        "standardized": bool(changed),
        "cellpar": [round(float(x), 8) for x in cp0],
        "n_atoms": len(refined),
        "elements": list(refined.get_chemical_symbols()),
    }
    return refined, info


def write_validated_cif(atoms, outdir: str | Path) -> Path:
    """Write the (possibly refined) structure as a P1 CIF for QE consumption."""
    from ase import io

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "validated.cif"
    io.write(str(path), atoms, format="cif")
    return path


if __name__ == "__main__":  # pragma: no cover
    if len(sys.argv) < 2:
        print("usage: qe_struct.py STRUCTURE")
        sys.exit(2)
    r = check_structure(sys.argv[1])
    print(format_check_result(r))
    sys.exit(1 if r.has_fail else 0)
