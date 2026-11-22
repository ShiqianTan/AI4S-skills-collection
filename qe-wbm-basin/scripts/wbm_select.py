#!/usr/bin/env python3
"""Select WBM structures and prepare MPRelax-like QE inputs.

Presets:
  base20   - the original 20 WBM-IS2RE structures (8 were basin-mismatched)
  medium20 - 20 additional moderate-distortion structures (40-structure test)
  high20   - 20 high-distortion structures (stress test)
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core import mprelax, qeio  # noqa: E402

BASE20 = [
    "wbm-1-28396", "wbm-1-2535", "wbm-3-38495", "wbm-1-2639", "wbm-2-26680",
    "wbm-3-39840", "wbm-2-26583", "wbm-3-39116", "wbm-1-47411", "wbm-5-11783",
    "wbm-1-927", "wbm-3-63782", "wbm-4-34304", "wbm-1-2696", "wbm-1-53069",
    "wbm-5-1670", "wbm-3-5763", "wbm-5-12529", "wbm-4-18850", "wbm-3-38323",
]
MEDIUM20 = [
    "wbm-3-22021", "wbm-2-15080", "wbm-1-50662", "wbm-1-10853", "wbm-2-15104",
    "wbm-1-23768", "wbm-1-20073", "wbm-3-73395", "wbm-1-9139", "wbm-1-40281",
    "wbm-1-28651", "wbm-1-10017", "wbm-3-50521", "wbm-3-29299", "wbm-1-50645",
    "wbm-5-11579", "wbm-1-9079", "wbm-4-20306", "wbm-3-13286", "wbm-2-30018",
]
HIGH20 = [
    "wbm-3-22021", "wbm-2-15080", "wbm-1-50662", "wbm-1-10853", "wbm-2-15104",
    "wbm-1-23768", "wbm-1-20073", "wbm-3-73395", "wbm-1-9139", "wbm-1-40281",
    "wbm-1-52993", "wbm-1-26087", "wbm-1-3410", "wbm-1-44555", "wbm-1-52998",
    "wbm-3-46925", "wbm-2-23965", "wbm-1-22407", "wbm-3-14775", "wbm-3-56995",
]
PRESETS = {"base20": BASE20, "medium20": MEDIUM20, "high20": HIGH20}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--preset", choices=list(PRESETS), default="base20")
    ap.add_argument("--wbm", default=None, help="matbench-discovery WBM directory")
    ap.add_argument("--is2re", default=None, help="WBM_IS2RE.aselmdb")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    from ase.db import connect
    from ase.io import write
    from pymatgen.core import Structure

    p = get_paths(args.config)
    wbm = Path(args.wbm) if args.wbm else p.dataset("wbm")
    is2re = Path(args.is2re) if args.is2re else p.dataset("wbm_is2re")
    out = Path(args.out) if args.out else p.work / (
        "wbm" if args.preset == "base20" else f"wbm_{args.preset}"
    )
    out.mkdir(parents=True, exist_ok=True)
    sssp = json.loads(p.pseudo_json("efficiency").read_text())
    ids = PRESETS[args.preset]
    selected = set(ids)

    summary = pd.read_csv(wbm / "2023-12-13-wbm-summary.csv.gz").set_index("material_id")
    vasp = {}
    with gzip.open(wbm / "2022-10-19-wbm-computed-structure-entries.jsonl.gz", "rt") as f:
        for line in f:
            rec = json.loads(line)
            mid = rec["material_id"]
            if mid not in selected:
                continue
            cse = rec["computed_structure_entry"]
            vasp[mid] = {
                "structure": Structure.from_dict(cse["structure"]),
                "energy": float(cse["energy"]),
                "potcar_symbols": cse.get("parameters", {}).get("potcar_symbols", []),
            }
            if len(vasp) == len(selected):
                break

    initial = {}
    db = connect(str(is2re))
    for row in db.select():
        mid = str(row.data.get("material_id", ""))
        if mid in selected and mid not in initial:
            initial[mid] = row.toatoms()
        if len(initial) == len(selected):
            break
    missing = sorted(selected - set(vasp) - set(initial))
    if missing:
        raise RuntimeError(f"missing WBM entries: {missing}")

    manifest = {"config": f"WBM_{args.preset}", "systems": {}}
    for mid in ids:
        row = summary.loc[mid]
        formula = str(row.formula).replace(" ", "")
        d = out / mid
        d.mkdir(parents=True, exist_ok=True)
        write(d / "initial.cif", initial[mid])
        init = Structure.from_str((d / "initial.cif").read_text(), fmt="cif")
        v = vasp[mid]
        v["structure"].to(filename=str(d / "vasp_final.cif"))

        ecutwfc, ecutrho = mprelax.mp_cutoffs(init, sssp)
        kpts, k_shift = mprelax.mp_kpoints(init)
        conv_thr, degauss = mprelax.mp_electrons(len(init))
        qeio.write_qe_input(
            d / "vc_relax.in", init, f"wbm_{mid}_relax", "vc-relax",
            ecutwfc, ecutrho, kpts, sssp,
            occupations="smearing", smearing="gaussian", degauss=degauss,
            conv_thr=conv_thr, mixing_beta=0.7, electron_maxstep=mprelax.NELM,
            nstep=mprelax.NSTEP, forc_conv_thr=mprelax.FORC_CONV_THR_RY_BOHR,
            etot_conv_thr=mprelax.ETOT_CONV_THR_RY, k_shift=k_shift,
        )
        ecutwfc_v, ecutrho_v = mprelax.mp_cutoffs(v["structure"], sssp)
        kpts_v, k_shift_v = mprelax.mp_kpoints(v["structure"])
        conv_thr_v, degauss_v = mprelax.mp_electrons(len(v["structure"]))
        qeio.write_qe_input(
            d / "scf_vasp_final.in", v["structure"], f"wbm_{mid}_vaspfinal", "scf",
            ecutwfc_v, ecutrho_v, kpts_v, sssp,
            occupations="smearing", smearing="gaussian", degauss=degauss_v,
            conv_thr=conv_thr_v, mixing_beta=0.7, electron_maxstep=mprelax.NELM,
            k_shift=k_shift_v,
        )
        manifest["systems"][mid] = {
            "material_id": mid, "formula": formula, "n_sites": int(row.n_sites),
            "dir": str(d.relative_to(out)), "vasp_final_energy_eV": v["energy"],
            "vasp_e_form_mp2020_eV_atom": float(row.e_form_per_atom_mp2020_corrected),
            "vasp_e_above_hull_eV_atom": float(row.e_above_hull_wbm),
            "vasp_bandgap_pbe_eV": float(row.bandgap_pbe),
            "vasp_volume_A3": float(row.volume),
            "initial_volume_A3": float(init.volume),
            "potcar_symbols": v["potcar_symbols"],
            "ecutwfc": ecutwfc, "ecutrho": ecutrho, "kpts": list(kpts),
            "k_shift": list(k_shift), "conv_thr": conv_thr, "degauss": degauss,
            "site_stats_fingerprint_init_final_norm_diff": (
                float(row.site_stats_fingerprint_init_final_norm_diff)
                if "site_stats_fingerprint_init_final_norm_diff" in row else None
            ),
            "config": "MPRelaxSet_c4998d92",
        }
    manifest["systems"] = {k: v for k, v in manifest["systems"].items()}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"prepared {len(manifest['systems'])} WBM systems ({args.preset}) in {out}")
    for mid, m in manifest["systems"].items():
        print(mid, m["formula"], "nat", m["n_sites"],
              "V_init", round(m["initial_volume_A3"], 2),
              "V_vasp", round(m["vasp_volume_A3"], 2), "ecut", m["ecutwfc"], "k", m["kpts"])


if __name__ == "__main__":
    main()
