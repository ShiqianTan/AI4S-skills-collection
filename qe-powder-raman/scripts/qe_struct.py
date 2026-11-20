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

import re
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
        self.layers: dict | None = None

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

    # ---- basic structural validity ----------------------------------------
    validate_structure(atoms, res)

    # ---- partial occupancy -------------------------------------------------
    _check_partial_occupancy(atoms, res)

    # ---- overlap -----------------------------------------------------------
    _check_overlap(atoms, res)

    # ---- symmetry ----------------------------------------------------------
    _check_symmetry(atoms, res)

    # ---- layered / 2D hints (E-05) -----------------------------------------
    _check_layers(atoms, res)

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


def validate_structure(atoms, res: StructureCheckResult) -> None:
    """Check that the structure is a valid 3D periodic crystal with finite,
    non-degenerate cell and recognizable elements. Appends FAIL/WARN to `res`."""
    import numpy as np

    if not atoms.get_pbc().all():
        res.add(
            "FAIL",
            "non_3d",
            "structure is not 3D periodic (pbc must be True on all axes)",
        )

    try:
        cell = atoms.get_cell()
        pos = atoms.get_positions()
    except Exception as exc:  # noqa: BLE001
        res.add("FAIL", "non_finite", "cannot read cell/positions: %s" % exc)
        return
    if not (np.isfinite(cell).all() and np.isfinite(pos).all()):
        res.add("FAIL", "non_finite", "cell or coordinates contain NaN/Inf")

    try:
        vol = float(atoms.get_volume())
    except Exception:  # noqa: BLE001
        vol = 0.0
    if not np.isfinite(vol) or vol <= 0.0:
        res.add("FAIL", "cell_volume", "cell volume is zero, negative, or undefined")

    if not atoms.get_chemical_symbols():
        res.add("FAIL", "no_elements", "structure has no recognizable elements")


def analyze_layers(atoms, *, slab_vacuum_ang: float = 8.0,
                   gap_min_ang: float = 2.5) -> dict:
    """Detect slab vacuum and layered/vdW gaps (E-05).

    For each lattice direction the atoms are projected onto the fractional
    coordinate and the largest empty gap is measured in Angstrom using the
    perpendicular width ``V / |a_j x a_k|``.  A gap >= 8 A is treated as slab
    vacuum; otherwise a gap >= 2.5 A bordered by the same element on both sides
    is treated as a possible van der Waals interlayer gap.
    """
    import numpy as np

    n = len(atoms)
    if n == 0:
        return {"is_slab": False, "is_layered": False, "layer_count": 0}
    cell = np.asarray(atoms.get_cell(), dtype=float)
    vol = abs(float(np.linalg.det(cell)))
    frac = atoms.get_scaled_positions(wrap=True)
    symbols = atoms.get_chemical_symbols()

    best = None
    heights: list[float] = []
    for axis in range(3):
        j, k = (axis + 1) % 3, (axis + 2) % 3
        area = float(np.linalg.norm(np.cross(cell[j], cell[k])))
        height = vol / area if area > 0.0 else 0.0
        heights.append(round(height, 6))
        s = frac[:, axis] % 1.0
        order = np.argsort(s)
        ss = s[order]
        if n == 1:
            gaps = [(height, 1.0, symbols[order[0]], symbols[order[0]],
                     int(order[0]), int(order[0]))]
        else:
            gaps = []
            for t in range(n):
                if t < n - 1:
                    g = float(ss[t + 1] - ss[t])
                else:
                    g = float(1.0 + ss[0] - ss[-1])
                lo = int(order[t])
                hi = int(order[(t + 1) % n])
                gaps.append((g * height, g, symbols[lo], symbols[hi], lo, hi))
        gap = max(gaps, key=lambda item: item[0])
        if best is None or gap[0] > best[0]:
            best = (gap[0], axis, height, gap, gaps)

    max_gap, axis, _height, gap, gaps = best
    is_slab = max_gap >= slab_vacuum_ang
    same_element = gap[2] == gap[3]
    is_layered = (not is_slab) and max_gap >= gap_min_ang and same_element

    internal = [
        g for g in gaps
        if g[0] >= gap_min_ang and g[2] == g[3] and g is not gap
    ]
    if is_slab:
        layer_count = max(1, len(internal))
    elif is_layered:
        layer_count = max(1, len(internal) + 1)
    else:
        layer_count = 1

    return {
        "axis": "abc"[axis],
        "axis_index": axis,
        "vacuum_thickness": round(float(max_gap), 4),
        "peaks": heights,
        "gap_elements": [gap[2], gap[3]],
        "layer_count": int(layer_count),
        "is_slab": bool(is_slab),
        "is_layered": bool(is_layered),
        "slab_vacuum_ang": slab_vacuum_ang,
        "gap_min_ang": gap_min_ang,
    }


def _check_layers(atoms, res: StructureCheckResult) -> None:
    """Emit E-05 warnings for slab vacuum / layered structures."""
    try:
        info = analyze_layers(atoms)
    except Exception as exc:  # noqa: BLE001
        res.add("WARN", "layers", "layer analysis failed: %s" % exc)
        return
    res.layers = info
    if info.get("is_slab"):
        res.add(
            "WARN",
            "slab_vacuum",
            "quasi-2D/slab detected: %.2f A vacuum along %s (%d layer(s)); "
            "Gamma-only bulk workflow, no interlayer coupling in SCF; relax "
            "in-plane parameters and check for soft interlayer modes"
            % (info["vacuum_thickness"], info["axis"], info["layer_count"]),
        )
    elif info.get("is_layered"):
        res.add(
            "WARN",
            "layered_gap",
            "possible layered/vdW structure: %.2f A %s-%s gap along %s "
            "(%d layer(s)); interlayer modes may be soft -> consider vc-relax "
            "and dispersion correction; Raman assumes a 3D crystal"
            % (info["vacuum_thickness"], info["gap_elements"][0],
               info["gap_elements"][1], info["axis"], info["layer_count"]),
        )


def _check_partial_occupancy(atoms, res: StructureCheckResult) -> None:
    # ASE CIF reader stores occupancies in atoms.info["occupancy"] as either a
    # flat {site: occ} or nested {site: {element: occ}}. Iterate the VALUES
    # (never the dict keys) so partial occupancies are actually detected.
    occ = atoms.info.get("occupancy")
    if occ is None:
        res.add("PASS", "occupancy", "no partial occupancies detected")
        return

    values = list(occ.values()) if isinstance(occ, dict) else list(occ)
    partial = 0
    mixed = 0
    for v in values:
        if isinstance(v, dict):
            if len(v) > 1:
                mixed += 1
            for o in v.values():
                try:
                    f = float(o)
                except (TypeError, ValueError):
                    f = 1.0
                if 0.0 < f < 1.0 - 1e-8:
                    partial += 1
        else:
            try:
                f = float(v)
            except (TypeError, ValueError):
                f = 1.0
            if 0.0 < f < 1.0 - 1e-8:
                partial += 1

    if mixed:
        res.add(
            "FAIL",
            "mixed_occupancy",
            "%d site(s) have mixed-species occupancy (not supported)" % mixed,
        )
    if partial:
        res.add(
            "FAIL",
            "occupancy",
            "%d site(s) have partial occupancy (site occupancy must be 1.0)" % partial,
        )
    if not mixed and not partial:
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
    fallback_reason = None
    cp0 = list(atoms.cell.cellpar())

    try:
        lat, pos, nums = spglib.refine_cell(cell, symprec=symprec)
        if lat is None:
            fallback_reason = "spglib.refine_cell returned no lattice"
        elif len(nums) != len(atoms):
            fallback_reason = (
                "refine_cell changed atom count (%d -> %d); kept original"
                % (len(atoms), len(nums))
            )
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
        fallback_reason = "refine_cell raised: %s" % exc

    info = {
        "original": original,
        "refined": refined_info,
        "standardized": bool(changed),
        "fallback_reason": fallback_reason,
        "cellpar": [round(float(x), 8) for x in cp0],
        "n_atoms": len(refined),
        "elements": list(refined.get_chemical_symbols()),
    }
    return refined, info


_BOHR_TO_ANG = 0.529177210903


def parse_final_structure(text: str):
    """Parse the relaxed structure from a pw.x relax/vc-relax output.

    Returns an ase.Atoms parsed from the ``Begin final coordinates`` ...
    ``End final coordinates`` block, or None if absent/unparseable.
    """
    import numpy as np
    from ase import Atoms

    m = re.search(r"Begin final coordinates(.*?)End final coordinates", text, re.S)
    if not m:
        return None
    block = m.group(1)

    cell = None
    cm = re.search(r"CELL_PARAMETERS\s*[\(\{](angstrom|bohr|alat)[\)\}]",
                   block, re.I)
    if cm:
        unit = cm.group(1).lower()
        rows = []
        for line in block[cm.end():].strip().splitlines()[:3]:
            parts = line.split()
            if len(parts) < 3:
                break
            try:
                rows.append([float(x) for x in parts[:3]])
            except ValueError:
                break
        if len(rows) == 3:
            cell = np.array(rows, dtype=float)
            if unit == "bohr":
                cell = cell * _BOHR_TO_ANG

    pm = re.search(r"ATOMIC_POSITIONS\s*[\(\{](angstrom|crystal|bohr|alat)[\)\}]",
                   block, re.I)
    if not pm or cell is None:
        return None
    unit = pm.group(1).lower()
    syms = []
    coords = []
    for line in block[pm.end():].strip().splitlines():
        parts = line.split()
        if len(parts) < 4:
            break
        try:
            coords.append([float(x) for x in parts[1:4]])
        except ValueError:
            break
        syms.append(parts[0])
    if not syms:
        return None
    coords = np.array(coords, dtype=float)
    if unit == "crystal":
        return Atoms(symbols=syms, scaled_positions=coords, cell=cell, pbc=True)
    if unit == "bohr":
        coords = coords * _BOHR_TO_ANG
    return Atoms(symbols=syms, positions=coords, cell=cell, pbc=True)


def relax_metrics(text: str) -> dict:
    """Extract residual force / stress from a relax output (best effort)."""
    out: dict = {"final_force": None, "final_stress": None}
    forces = re.findall(r"Total force =\s*([-\d.eE+]+)", text)
    if forces:
        try:
            out["final_force"] = float(forces[-1])
        except ValueError:
            pass
    stresses = re.findall(r"Total stress \(Ry/bohr\*\*3\)\s*=?\s*([-\d.eE+]+)", text)
    if stresses:
        try:
            out["final_stress"] = float(stresses[-1])
        except ValueError:
            pass
    return out


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
