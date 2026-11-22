"""Quantum ESPRESSO input generation helpers.

Extracted from ``repro/scripts/mp20_prepare.py`` and
``repro/scripts/exp104_common.py`` so every skill writes consistent inputs.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping


def cutoff_for(structure, sssp: Mapping[str, Mapping[str, float]]) -> tuple[float, float]:
    """Return (ecutwfc, ecutrho) as the element-wise SSSP recommendation."""
    elements = [str(e) for e in structure.composition.elements]
    ecutwfc = max(sssp[e]["cutoff_wfc"] for e in elements)
    ecutrho = max(sssp[e]["cutoff_rho"] for e in elements)
    return float(ecutwfc), float(ecutrho)


def kgrid(structure, spacing: float, minimum: int = 1) -> tuple[int, int, int]:
    """Gamma-centred k-grid from a real-space spacing in 1/Angstrom."""
    grid = []
    for length in structure.lattice.abc:
        n = int(round(2.0 * math.pi / (length * spacing)))
        grid.append(max(minimum, n))
    return tuple(grid)


def _pseudo_filename(pseudo_map: Mapping[str, object], element: str) -> str:
    rec = pseudo_map[element]
    if isinstance(rec, Mapping):
        return str(rec["filename"])
    return str(rec)


def write_qe_input(
    path: str | Path,
    structure,
    prefix: str,
    calc: str,
    ecutwfc: float,
    ecutrho: float,
    kpts: tuple[int, int, int],
    pseudo_map: Mapping[str, object],
    pseudo_dir: str | Path = "./pseudo",
    outdir: str = "./tmp",
    occupations: str = "fixed",
    nbnd: int | None = None,
    conv_thr: float = 1.0e-8,
    diagonalization: str | None = None,
    diago_thr_init: float | None = None,
    diago_full_acc: bool = False,
    diago_cg_maxiter: int | None = None,
    input_dft: str | None = None,
    smearing: str = "mv",
    degauss: float = 0.02,
    nspin: int = 1,
    mixing_mode: str | None = None,
    mixing_beta: float = 0.7,
    electron_maxstep: int = 300,
    nstep: int | None = None,
    forc_conv_thr: float | None = None,
    etot_conv_thr: float | None = None,
    k_shift: tuple[int, int, int] = (0, 0, 0),
    nosym: bool = False,
    ion_dynamics: str = "bfgs",
    ion_damping: float | None = None,
    cell_dynamics: str = "bfgs",
    cell_damping: float | None = None,
    press_conv_thr: float = 0.5,
    kpoints_path: list[tuple[float, float, float]] | None = None,
) -> Path:
    """Write a pw.x input file and return its path."""
    from pymatgen.core import Element

    path = Path(path)
    species = sorted({str(s) for s in structure.species}, key=lambda x: Element(x).Z)
    lines = [
        "&CONTROL",
        f"  calculation = '{calc}'",
        f"  prefix = '{prefix}'",
        f"  pseudo_dir = '{pseudo_dir}'",
        f"  outdir = '{outdir}'",
        "  verbosity = 'high'",
        "  tprnfor = .true.",
        "  tstress = .true.",
        "/",
        "&SYSTEM",
        "  ibrav = 0",
        f"  nat = {len(structure)}",
        f"  ntyp = {len(species)}",
        f"  ecutwfc = {ecutwfc:.1f}",
        f"  ecutrho = {ecutrho:.1f}",
        f"  occupations = '{occupations}'",
    ]
    if occupations == "smearing":
        lines += [f"  smearing = '{smearing}'", f"  degauss = {degauss}"]
    if nspin != 1:
        lines.append(f"  nspin = {nspin}")
    if nosym:
        lines.append("  nosym = .true.")
    if input_dft is not None:
        lines.append(f"  input_dft = '{input_dft}'")
    if nbnd is not None:
        lines.append(f"  nbnd = {nbnd}")
    lines += [
        "/",
        "&ELECTRONS",
        f"  conv_thr = {conv_thr:.1e}",
        f"  mixing_beta = {mixing_beta}",
        f"  electron_maxstep = {electron_maxstep}",
    ]
    if mixing_mode is not None:
        lines.append(f"  mixing_mode = '{mixing_mode}'")
    if diagonalization is not None:
        lines.append(f"  diagonalization = '{diagonalization}'")
    if diago_thr_init is not None:
        lines.append(f"  diago_thr_init = {diago_thr_init:.1e}")
    if diago_full_acc:
        lines.append("  diago_full_acc = .true.")
    if diago_cg_maxiter is not None:
        lines.append(f"  diago_cg_maxiter = {diago_cg_maxiter}")
    lines += ["/"]
    control_end = lines.index("/")
    extra_control = []
    if nstep is not None:
        extra_control.append(f"  nstep = {nstep}")
    if forc_conv_thr is not None:
        extra_control.append(f"  forc_conv_thr = {forc_conv_thr:.3e}")
    if etot_conv_thr is not None:
        extra_control.append(f"  etot_conv_thr = {etot_conv_thr:.3e}")
    for line in reversed(extra_control):
        lines.insert(control_end, line)
    if calc in ("relax", "vc-relax"):
        lines += ["&IONS", f"  ion_dynamics = '{ion_dynamics}'"]
        if ion_damping is not None:
            lines.append(f"  ion_damping = {ion_damping}")
        lines += ["/"]
    if calc == "vc-relax":
        lines += ["&CELL", f"  cell_dynamics = '{cell_dynamics}'"]
        if cell_damping is not None:
            lines.append(f"  cell_damping = {cell_damping}")
        lines.append(f"  press_conv_thr = {press_conv_thr}")
        lines += ["/"]
    lines += ["ATOMIC_SPECIES"]
    for el in species:
        lines.append(
            f"  {el:2s} {Element(el).atomic_mass:10.4f} {_pseudo_filename(pseudo_map, el)}"
        )
    lines += ["CELL_PARAMETERS angstrom"]
    for vec in structure.lattice.matrix:
        lines.append("  " + " ".join(f"{x:16.10f}" for x in vec))
    lines += ["ATOMIC_POSITIONS crystal"]
    for site in structure:
        lines.append(
            f"  {str(site.specie):2s} "
            + " ".join(f"{x:16.10f}" for x in site.frac_coords)
        )
    if kpoints_path is not None:
        lines += ["K_POINTS crystal", f"  {len(kpoints_path)}"]
        for k in kpoints_path:
            lines.append(f"  {k[0]:.10f} {k[1]:.10f} {k[2]:.10f} 1.0")
    else:
        lines += [
            "K_POINTS automatic",
            f"  {kpts[0]} {kpts[1]} {kpts[2]} {k_shift[0]} {k_shift[1]} {k_shift[2]}",
            "",
        ]
    path.write_text("\n".join(lines))
    return path
