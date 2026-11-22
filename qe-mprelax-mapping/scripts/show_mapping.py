#!/usr/bin/env python3
"""Print the MPRelaxSet -> QE mapping for a structure."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mprelax_config import mapping_for  # noqa: E402


def _load_sssp(json_path: str | None, elements: list[str]) -> dict:
    if json_path:
        return json.loads(Path(json_path).read_text())
    # Fall back to the shared config, if available.
    for parent in [Path(__file__).resolve(), *Path(__file__).resolve().parents]:
        if (parent / "_shared").is_dir():
            sys.path.insert(0, str(parent / "_shared"))
            try:
                from qe_bench_core import paths

                cfg = paths.load_config()
                return json.loads(paths.pseudo_json(cfg, "efficiency").read_text())
            except Exception:
                break
    # Last resort: a minimal manifest with no overrides.
    return {el: {"cutoff_wfc": 40.0, "cutoff_rho": 320.0} for el in elements}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("structure", help="CIF / POSCAR structure")
    ap.add_argument("--sssp-json", default=None, help="SSSP efficiency JSON manifest")
    ap.add_argument("--json", action="store_true", help="print JSON")
    args = ap.parse_args()

    from pymatgen.core import Structure

    structure = Structure.from_file(args.structure)
    elements = [str(e) for e in structure.composition.elements]
    sssp = _load_sssp(args.sssp_json, elements)
    mapping = mapping_for(structure, sssp)

    if args.json:
        print(json.dumps(mapping, indent=2))
        return

    print(f"system: {mapping['formula']}  natoms: {mapping['natoms']}")
    print(f"ecutwfc (Ry): {mapping['ecutwfc_Ry']}")
    print(f"ecutrho (Ry): {mapping['ecutrho_Ry']}")
    print(f"conv_thr (Ry): {mapping['conv_thr_Ry']:.3e}")
    print(f"smearing: {mapping['smearing']}  degauss (Ry): {mapping['degauss_Ry']:.6g}")
    print(f"k-grid: {mapping['k_grid']} shift={mapping['k_shift']}")
    print(f"electron_maxstep: {mapping['electron_maxstep']}  nstep: {mapping['nstep']}")
    print(f"calculation: {mapping['calculation']} ({mapping['ion_dynamics']}/{mapping['cell_dynamics']})")


if __name__ == "__main__":
    main()
