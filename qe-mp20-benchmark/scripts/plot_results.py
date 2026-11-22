#!/usr/bin/env python3
"""Parity plots for the MP20 benchmark (gap, volume, formation energy)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402
from _bootstrap import add_config_arg, get_paths  # noqa: E402
from qe_bench_core.plotting import parity_plot  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_config_arg(ap)
    ap.add_argument("--tag", default="mp20_20")
    args = ap.parse_args()
    p = get_paths(args.config)
    results = p.results

    mp = pd.read_csv(results / f"{args.tag}_results.csv")
    ok = mp[mp.status == "ok"].copy()
    parity_plot(ok.MP_gap_eV, ok.QE_gap_eV, results / f"{args.tag}_gap_parity.png",
                "MP (VASP) gap", "QE (PBE) gap", f"{args.tag} band gap", "eV")
    parity_plot(ok.MP_volume_A3, ok.QE_volume_A3, results / f"{args.tag}_volume_parity.png",
                "MP volume", "QE volume", f"{args.tag} relaxed volume", "Å^3")
    ef = ok[ok.QE_Ef_raw_eV_atom.notna()]
    parity_plot(ef.MP_Ef_eV_atom, ef.QE_Ef_raw_eV_atom,
                results / f"{args.tag}_formation_energy_raw.png",
                "MP Ef", "QE raw Ef", f"{args.tag} formation energy (raw)", "eV/atom")

    mp2020 = results / f"{args.tag}_mp2020.csv"
    if mp2020.is_file():
        corr = pd.read_csv(mp2020)
        parity_plot(corr.MP_Ef_eV_atom, corr.QE_Ef_MP2020_eV_atom,
                    results / f"{args.tag}_formation_energy_mp2020.png",
                    "MP Ef", "QE + MP2020 Ef",
                    f"{args.tag} formation energy (MP2020)", "eV/atom")
    print("wrote plots to", results)


if __name__ == "__main__":
    main()
