"""Parsers for pw.x output: energies, forces, stress, band gaps, bands."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

RY_TO_EV = 13.605693122994
BOHR_TO_ANG = 0.529177210903
KBAR_TO_EV_A3 = 0.0006241509


def voigt_to_matrix(v) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    return np.array([[v[0], v[5], v[4]], [v[5], v[1], v[3]], [v[4], v[3], v[2]]])


def last_total_energy_ry(outfile: str | Path) -> float | None:
    text = Path(outfile).read_text(errors="ignore")
    energies = re.findall(r"!\s+total energy\s+=\s+([-\d.]+)\s+Ry", text)
    return float(energies[-1]) if energies else None


def parse_scf_text(outfile: str | Path, nat: int | None = None) -> dict:
    """Parse total energy / forces / stress directly from pw.x text output."""
    text = Path(outfile).read_text(errors="ignore")
    energies = [float(x) for x in re.findall(r"!\s+total energy\s+=\s+([-\d.]+)\s+Ry", text)]
    if not energies:
        energies = [float(x) for x in re.findall(r"total energy\s+=\s+([-\d.]+)\s+Ry", text)]
    energy = energies[-1] * RY_TO_EV if energies else None

    forces = np.empty((0, 3))
    for line in text.splitlines():
        if "Forces acting on atoms" in line:
            forces = np.empty((0, 3))
            continue
        match = re.match(
            r"\s*atom\s+\d+\s+type\s+\d+\s+force\s*=\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)",
            line,
        )
        if match is not None and (nat is None or len(forces) < nat):
            forces = np.vstack([forces, [float(x) for x in match.groups()]])
    if len(forces):
        forces = forces * RY_TO_EV / BOHR_TO_ANG

    stress = None
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if re.search(r"total\s+stress\s+\(Ry/bohr\*\*3\)", line):
            vals = []
            for follow in lines[i + 1 : i + 8]:
                parts = follow.split()
                if len(parts) >= 3:
                    try:
                        vals.append([float(x) for x in parts[:3]])
                    except ValueError:
                        pass
                if len(vals) == 3:
                    break
            if len(vals) == 3:
                stress = np.asarray(vals, dtype=float) * RY_TO_EV / BOHR_TO_ANG**3
                break
    return {"energy_eV": energy, "forces": forces, "stress_eV_A3": stress}


def scf_results(outfile: str | Path) -> dict:
    """Parse a finished SCF output with ASE's espresso-out reader."""
    from ase.io import read
    import numpy as np
    from qe_bench_core.structure import structure_from_ase

    atoms = read(str(outfile), format="espresso-out")
    return {
        "atoms": atoms,
        "structure": structure_from_ase(atoms),
        "energy_eV": float(atoms.get_potential_energy()),
        "forces": np.asarray(atoms.get_forces(), dtype=float),
        "stress": np.asarray(atoms.get_stress(), dtype=float),
        "volume_A3": float(atoms.get_volume()),
    }


def gap_from_scf(outfile: str | Path) -> float | None:
    gap = None
    for line in Path(outfile).read_text(errors="ignore").splitlines():
        if "highest occupied, lowest unoccupied level" in line:
            vals = line.split(":")[1].split()
            gap = float(vals[1]) - float(vals[0])
        elif "highest occupied level" in line and "lowest unoccupied" not in line:
            gap = 0.0
    return gap


def parse_bands_blocks(outfile: str | Path) -> list[list[float]]:
    blocks: list[list[float]] = []
    cur: list[float] | None = None
    for line in Path(outfile).read_text(errors="ignore").splitlines():
        if "bands (ev)" in line:
            if cur:
                blocks.append(cur)
            cur = []
            continue
        if cur is None:
            continue
        stripped = line.strip()
        if not stripped:
            continue
        try:
            vals = [float(x) for x in stripped.split()]
        except ValueError:
            if cur:
                blocks.append(cur)
            cur = None
            continue
        cur.extend(vals)
    if cur:
        blocks.append(cur)
    return blocks


def gap_from_bands(outfile: str | Path, nocc: int) -> dict:
    blocks = parse_bands_blocks(outfile)
    empty = {
        "n_kblocks": len(blocks),
        "vbm_eV": None,
        "cbm_eV": None,
        "gap_eV": None,
        "direct_gap_eV": None,
    }
    if not blocks or nocc <= 0:
        return empty
    usable = [b for b in blocks if len(b) > nocc]
    if not usable:
        return empty
    vbm = max(max(b[:nocc]) for b in usable)
    cbm = min(min(b[nocc:]) for b in usable)
    direct = min(min(b[nocc:]) - max(b[:nocc]) for b in usable)
    return {
        "n_kblocks": len(blocks),
        "vbm_eV": float(vbm),
        "cbm_eV": float(cbm),
        "gap_eV": float(max(0.0, cbm - vbm)),
        "direct_gap_eV": float(max(0.0, direct)),
    }
