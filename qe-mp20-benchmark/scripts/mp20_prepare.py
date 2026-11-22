#!/usr/bin/env python3
"""Prepare MP20 inputs from bundled CIFs and the MP reference manifest.

Two configuration modes are supported:

* ``sssp``    - SSSP 1.3.0 PBE efficiency cutoffs + 0.25/0.20/0.15 A^-1 k-grids
* ``mprelax`` - MPRelaxSet-like settings (ENCUT 520 eV floor, reciprocal
                density 64, EDIFF 5e-5 eV/atom, SIGMA 0.05 eV, NSW 99)

The elemental reference phases are prepared alongside the systems so the
formation energies can be referenced to QE-relaxed elemental ground states.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core import mprelax, pseudos, qeio  # noqa: E402
from qe_bench_core.structure import read_structure  # noqa: E402

SYSTEMS = {
    "mp-2172": "AlAs", "mp-2176": "ZnTe", "mp-22925": "AgI",
    "mp-2624": "AlSb", "mp-1039": "MgTe", "mp-2201": "PbSe",
    "mp-1057015": "PbS", "mp-13032": "MgS", "mp-1115": "PtSe2",
    "mp-12779": "CdTe", "mp-1342": "BaO", "mp-2472": "SrO",
    "mp-971": "K2O", "mp-3748": "AlCuO2", "mp-830": "GaN",
    "mp-22205": "InN", "mp-2828": "HfN", "mp-22922": "AgCl",
    "mp-22903": "RbI", "mp-763": "MgB2",
}
SYSTEMS_5 = ["mp-2172", "mp-2176", "mp-22925", "mp-2624", "mp-1039"]
REFERENCE_IDS = {
    "Al": "mp-134", "As": "mp-11", "Zn": "mp-79", "Te": "mp-19",
    "Ag": "mp-989737", "I": "mp-23153", "Sb": "mp-104", "Mg": "mp-1056702",
    "Pb": "mp-20483", "Se": "mp-14", "Pt": "mp-126", "Cd": "mp-1096861",
    "B": "mp-160", "Rb": "mp-1179656", "S": "mp-557376",
}


# --- compatibility helpers (also imported by run_mp20) ----------------------
def load_sssp(json_path: str | Path | None = None) -> dict:
    if json_path is None:
        raise ValueError("SSSP JSON path is required")
    return pseudos.load_sssp(json_path)


def z_valence(filename: str, pseudo_dir: str | Path | None = None) -> float:
    if pseudo_dir is None:
        raise ValueError("pseudo_dir is required to read z_valence")
    return pseudos.z_valence(Path(pseudo_dir) / filename)


def cutoff_for(structure, sssp: dict) -> tuple[float, float]:
    return qeio.cutoff_for(structure, sssp)


def kgrid(structure, spacing: float, minimum: int = 1):
    return qeio.kgrid(structure, spacing, minimum)


def write_qe_input(
    path, structure, prefix, calc, ecutwfc, ecutrho, kpts, sssp,
    pseudo_dir=None, **kwargs,
):
    return qeio.write_qe_input(
        path, structure, prefix, calc, ecutwfc, ecutrho, kpts, sssp,
        pseudo_dir=pseudo_dir or "./pseudo", **kwargs,
    )


def _filename_map(sssp: dict, elements: list[str]) -> dict:
    return {el: {"filename": sssp[el]["filename"]} for el in elements}


def _prepare_system(mp_id, meta, structure, sssp, mode, pseudo_dir, out):
    formula = meta["formula"]
    elements = [str(e) for e in structure.composition.elements]
    d = out / f"{mp_id}_{formula}"
    d.mkdir(parents=True, exist_ok=True)
    pmap = _filename_map(sssp, elements)

    if mode == "sssp":
        ecutwfc, ecutrho = meta["ecutwfc"], meta["ecutrho"]
        k_relax = tuple(meta["k_relax"])
        k_scf = tuple(meta.get("k_scf", k_relax))
        k_nscf = tuple(meta["k_nscf"])
        nbnd = meta["nbnd"]
        occupations = meta.get("occupations", "fixed")
        nscf_occ = occupations
        nscf_diag = meta.get("nscf_diagonalization")
        k_shift = (0, 0, 0)
        settings = dict(
            occupations=occupations, conv_thr=1.0e-8, mixing_beta=0.7,
            electron_maxstep=300, nbnd=nbnd, diagonalization=nscf_diag,
        )
        pre_nscf_scf = False
        nosym = False
    else:
        ecutwfc, ecutrho = mprelax.mp_cutoffs(structure, sssp)
        k_relax, k_shift = mprelax.mp_kpoints(structure)
        k_scf = k_relax
        k_nscf = kgrid(structure, spacing=0.15)
        if formula == "GaN":
            k_nscf = k_relax
        pmap = pseudos.sssp_map(sssp, elements, pseudo_dir)
        nbnd = int(round(pseudos.n_electrons(structure, pmap))) // 2 + 10
        conv_thr, degauss = mprelax.mp_electrons(len(structure))
        metal = bool(meta.get("metal", False))
        occupations = "smearing"
        nscf_occ = "smearing" if metal else "fixed"
        nscf_diag = "cg" if formula == "GaN" else None
        nosym = formula == "GaN"
        settings = dict(
            occupations=occupations, smearing="gaussian", degauss=degauss,
            conv_thr=conv_thr, mixing_beta=mprelax.MP_SETTINGS["mixing_beta"],
            electron_maxstep=mprelax.NELM, nstep=mprelax.NSTEP,
            forc_conv_thr=mprelax.FORC_CONV_THR_RY_BOHR,
            etot_conv_thr=mprelax.ETOT_CONV_THR_RY,
            k_shift=k_shift, nosym=nosym,
        )
        pre_nscf_scf = True

    write_qe_input(
        d / "vc_relax.in", structure, f"{mp_id}_relax", "vc-relax",
        ecutwfc, ecutrho, k_relax, pmap, pseudo_dir=pseudo_dir, **settings,
    )
    record = {
        **{k: meta.get(k) for k in (
            "formula", "nat", "elements", "metal", "mp_formation_energy_per_atom",
            "mp_band_gap", "mp_e_above_hull", "mp_spacegroup", "cif_volume",
            "formation_energy_available",
        )},
        "ecutwfc": ecutwfc, "ecutrho": ecutrho,
        "k_relax": list(k_relax), "k_scf": list(k_scf), "k_nscf": list(k_nscf),
        "k_shift": list(k_shift), "nbnd": nbnd,
        "occupations": occupations, "nscf_occupations": nscf_occ,
        "nscf_diagonalization": nscf_diag, "pre_nscf_scf": pre_nscf_scf,
        "nosym": nosym, "config": mode,
    }
    return record


def _prepare_reference(el, mp_id, structure, sssp, mode, pseudo_dir, out, meta):
    elements = [str(e) for e in structure.composition.elements]
    pmap = _filename_map(sssp, elements)
    if mode == "sssp":
        ecutwfc, ecutrho = meta["ecutwfc"], meta["ecutrho"]
        k_relax = tuple(meta["k_relax"])
        k_shift = (0, 0, 0)
        settings = dict(occupations="smearing", smearing="mv", degauss=0.02)
    else:
        ecutwfc, ecutrho = mprelax.mp_cutoffs(structure, sssp)
        k_relax, k_shift = mprelax.mp_kpoints(structure)
        pmap = pseudos.sssp_map(sssp, elements, pseudo_dir)
        conv_thr, degauss = mprelax.mp_electrons(len(structure))
        settings = dict(
            occupations="smearing", smearing="gaussian", degauss=degauss,
            conv_thr=conv_thr, mixing_beta=mprelax.MP_SETTINGS["mixing_beta"],
            electron_maxstep=mprelax.NELM, nstep=mprelax.NSTEP,
            forc_conv_thr=mprelax.FORC_CONV_THR_RY_BOHR,
            etot_conv_thr=mprelax.ETOT_CONV_THR_RY,
            k_shift=k_shift, nosym=True,
        )
    d = out / "references" / el
    d.mkdir(parents=True, exist_ok=True)
    write_qe_input(
        d / "vc_relax.in", structure, f"ref_{el}_relax", "vc-relax",
        ecutwfc, ecutrho, k_relax, pmap, pseudo_dir=pseudo_dir, **settings,
    )
    return {
        "mp_id": mp_id, "nat": len(structure), "ecutwfc": ecutwfc,
        "ecutrho": ecutrho, "k_relax": list(k_relax), "k_shift": list(k_shift),
        "mp_formation_energy_per_atom": meta.get("mp_formation_energy_per_atom"),
        "cif_volume": meta.get("cif_volume"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--mode", choices=["sssp", "mprelax"], default="sssp")
    ap.add_argument("--systems", choices=["5", "20"], default="20")
    ap.add_argument("--data-dir", default=None, help="directory with bundled CIFs + manifest")
    ap.add_argument("--out", default=None, help="output directory for QE inputs")
    ap.add_argument("--pseudo-json", default=None)
    ap.add_argument("--pseudo-dir", default=None)
    ap.add_argument("--skip-references", action="store_true")
    args = ap.parse_args()

    p = get_paths(args.config)
    data_dir = Path(args.data_dir) if args.data_dir else p.skill_root / "data" / "mp20_20"
    out = Path(args.out) if args.out else p.work / (
        "mp20_20" if args.mode == "sssp" else "mp20_20_mprelax"
    )
    out.mkdir(parents=True, exist_ok=True)

    ref_manifest = json.loads((data_dir / "manifest.json").read_text())
    sssp_json = Path(args.pseudo_json) if args.pseudo_json else p.pseudo_json("efficiency")
    pseudo_dir = Path(args.pseudo_dir) if args.pseudo_dir else p.pseudo_dir("sssp_efficiency")
    sssp = load_sssp(sssp_json)

    ids = SYSTEMS_5 if args.systems == "5" else list(SYSTEMS)
    manifest = {"config": args.mode, "systems": {}, "references": {}}
    for mp_id in ids:
        formula = SYSTEMS[mp_id]
        cif = data_dir / f"{mp_id}_{formula}.cif"
        structure = read_structure(cif)
        manifest["systems"][mp_id] = _prepare_system(
            mp_id, ref_manifest["systems"][mp_id], structure, sssp,
            args.mode, pseudo_dir, out,
        )

    if not args.skip_references:
        for el, mp_id in REFERENCE_IDS.items():
            if args.systems == "5":
                needed = set()
                for mp_id2 in ids:
                    needed.update(ref_manifest["systems"][mp_id2].get("elements", []))
                if el not in needed:
                    continue
            cif = data_dir / f"ref_{el}_{mp_id}.cif"
            if not cif.is_file():
                continue
            structure = read_structure(cif)
            manifest["references"][el] = _prepare_reference(
                el, mp_id, structure, sssp, args.mode, pseudo_dir, out,
                ref_manifest["references"].get(el, {}),
            )

    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"prepared {len(manifest['systems'])} systems and "
          f"{len(manifest['references'])} references ({args.mode}) in {out}")
    for k, v in manifest["systems"].items():
        print(k, v["formula"], "ecut", v["ecutwfc"], v["ecutrho"],
              "k", v["k_relax"], "shift", v["k_shift"], "metal", v["metal"])


if __name__ == "__main__":
    main()
