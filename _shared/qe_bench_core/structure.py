"""Structure input/output helpers (pymatgen primary, ASE fallback)."""

from __future__ import annotations

from pathlib import Path


def read_structure(path: str | Path):
    """Read a structure from CIF/POSCAR/QE input using pymatgen."""
    from pymatgen.core import Structure

    return Structure.from_file(str(path))


def structure_from_cif(text: str):
    from pymatgen.core import Structure

    return Structure.from_str(text, fmt="cif")


def structure_from_ase(atoms):
    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor.get_structure(atoms)


def structure_to_ase(structure):
    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor.get_atoms(structure)


def write_cif(structure, path: str | Path) -> None:
    from pymatgen.io.cif import CifWriter

    Path(path).write_text(str(CifWriter(structure)))


def element_symbols(structure) -> list[str]:
    return sorted({str(s) for s in structure.species})
