"""Materials Project 2020 compatibility corrections for QE formation energies."""

from __future__ import annotations

from pymatgen.core import Composition
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedEntry


def correction_per_atom(comp: Composition, energy_ev: float) -> tuple[float, list[tuple[str, float]]]:
    """Return (correction eV/atom, list of (name, value) adjustments)."""
    entry = ComputedEntry(comp, energy_ev)
    entry.parameters["run_type"] = "GGA"
    entry.parameters["potcar_symbols"] = [str(el) for el in comp.elements]
    compat = MaterialsProject2020Compatibility(check_potcar=False)
    adjustments = compat.get_adjustments(entry)
    total = sum(a.value for a in adjustments)
    return total / comp.num_atoms, [(a.name, round(a.value, 4)) for a in adjustments]
