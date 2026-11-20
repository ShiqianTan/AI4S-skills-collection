#!/usr/bin/env python3
"""qe_parse — parse QE powder-Raman outputs (scf.out / ph.out / dynmat.out).

Produces ``raman_modes.csv`` (one row per mode, incl. LO tables) and
``summary.json`` in the run directory, and prints a human-readable summary.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

_OPTICAL_EPS = 1e-4  # activity threshold for "IR/Raman active"


# ---------------------------------------------------------------------------
# low-level parsers
# ---------------------------------------------------------------------------

def parse_scf(txt: str) -> dict:
    d: dict = {}

    def _f(pattern: str, key: str, cast=float):
        m = re.search(pattern, txt)
        if m:
            try:
                d[key] = cast(m.group(1))
            except ValueError:
                d[key] = None

    _f(r"number of electrons\s*=\s*([\d.]+)", "nelec")
    _f(r"number of Kohn-Sham states=\s*(\d+)", "nbnd", int)
    _f(r"number of k points=\s*(\d+)", "nk", int)
    _f(r"kinetic-energy cutoff\s*=\s*([\d.]+)", "ecutwfc")
    _f(r"(\d+)\s+Sym\. Ops\.", "sym_ops", int)
    _f(r"convergence has been achieved in\s+(\d+) iterations", "scf_iter", int)

    m = re.search(r"unit-cell volume\s*=\s*([\d.]+)", txt)
    if m:
        try:
            d["volume_a3"] = round(float(m.group(1)) * 0.148184711, 3)
        except ValueError:
            d["volume_a3"] = None

    m = re.search(r"^!\s+total energy\s+=\s+([-\d.]+)\s+Ry", txt, re.M)
    d["etot_ry"] = float(m.group(1)) if m else None

    m = re.search(r"highest occupied level \(ev\):\s*([-\d.]+)", txt)
    d["homo_ev"] = float(m.group(1)) if m else None

    return d


def parse_dielectric(txt: str) -> list[list[float]] | None:
    m = re.search(
        r"Dielectric constant in cartesian axis\s*\n(.*?\n.*?\n.*?\n)", txt, re.S
    )
    if not m:
        return None
    nums = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", m.group(1))
    try:
        vals = [float(x) for x in nums]
    except ValueError:
        return None
    if len(vals) < 9:
        return None
    return [vals[0:3], vals[3:6], vals[6:9]]


def parse_asr(txt: str) -> float | None:
    m = re.search(
        r"Acoustic Sum Rule:\s*\|\|dyn\(ASR\) - dyn\(orig\)\|\|=\s*([-\d.Ee+]+)", txt
    )
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return None
    return None


def parse_dynmat(txt: str) -> dict:
    d: dict = {}

    m = re.search(r"multiply by\s+([\d.]+)\s+for Clausius-Mossotti", txt)
    d["cm_factor"] = float(m.group(1)) if m else None

    m = re.search(r"multiply Raman by\s+([\d.]+)\s+for Clausius-Mossotti", txt)
    d["raman_cm_factor"] = float(m.group(1)) if m else None

    # Polarizability (A^3 units): header line, then (optionally) a "multiply"
    # line, then 3 lines of 3 floats. The reference analyzer leaked the
    # cm_factor into this matrix; parse it correctly here.
    d["polarizability"] = None
    lines = txt.splitlines()
    for i, ln in enumerate(lines):
        if "Polarizability" in ln and "A^3" in ln:
            row_start = i + 1
            # skip any "multiply ... Clausius-Mossotti" lines
            while row_start < len(lines) and "Clausius-Mossotti" in lines[row_start]:
                row_start += 1
            rows: list[list[float]] = []
            for j in range(row_start, min(row_start + 3, len(lines))):
                nums = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", lines[j])
                if len(nums) < 3:
                    break
                rows.append([float(x) for x in nums[:3]])
            if len(rows) == 3:
                d["polarizability"] = rows
            break

    # mode table
    d["has_raman"] = False
    d["modes"] = []
    started = False
    for ln in lines:
        if ln.lstrip().startswith("# mode"):
            started = True
            d["has_raman"] = "Raman" in ln
            continue
        if not started:
            continue
        s = ln.split()
        if not s:
            continue
        if re.fullmatch(r"-?\d+", s[0]) and len(s) in (4, 5, 6):
            try:
                mode = int(s[0])
                freq = float(s[1])
                thz = float(s[2])
                ir = float(s[3])
                raman = float(s[4]) if len(s) >= 6 else 0.0
                depol = float(s[5]) if len(s) >= 6 else None
            except ValueError:
                continue
            d["modes"].append(
                {"mode": mode, "freq": freq, "thz": thz, "ir": ir,
                 "raman": raman, "depol": depol}
            )
        else:
            # first non-mode line after the table starts -> end of table
            if d["modes"]:
                break
    return d


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def _stats(modes: list[dict]) -> dict:
    if not modes:
        return {
            "total": 0, "acoustic": [], "optical": [], "raman_active": [],
            "ir_active": [], "sum_raman": 0.0, "max_raman": None,
            "max_raman_freq": None, "acoustic_max_abs": 0.0,
            "optic_min": None, "optic_max": None, "soft_modes": [],
        }
    idx = sorted(range(len(modes)), key=lambda i: abs(modes[i]["freq"]))[:3]
    acoustic_idx = set(idx)
    acoustic = [modes[i] for i in idx]
    optical = [m for i, m in enumerate(modes) if i not in acoustic_idx]

    raman_active = [m for m in optical if m["raman"] > _OPTICAL_EPS]
    ir_active = [m for m in optical if m["ir"] > _OPTICAL_EPS]
    soft_modes = [m for m in optical if m["freq"] < 0.0]

    sum_raman = round(sum(m["raman"] for m in optical), 4)
    max_raman = max((m for m in optical if m["raman"] > _OPTICAL_EPS),
                    key=lambda m: m["raman"], default=None)
    optic_freqs = [m["freq"] for m in optical]

    return {
        "total": len(modes),
        "acoustic": [m["mode"] for m in acoustic],
        "optical": [m["mode"] for m in optical],
        "raman_active": [m["mode"] for m in raman_active],
        "ir_active": [m["mode"] for m in ir_active],
        "sum_raman": sum_raman,
        "max_raman": None if max_raman is None else round(max_raman["raman"], 4),
        "max_raman_freq": None if max_raman is None else max_raman["freq"],
        "acoustic_max_abs": round(max(abs(m["freq"]) for m in acoustic), 6),
        "optic_min": None if not optic_freqs else min(optic_freqs),
        "optic_max": None if not optic_freqs else max(optic_freqs),
        "soft_modes": [m["mode"] for m in soft_modes],
    }


# ---------------------------------------------------------------------------
# outdir-level
# ---------------------------------------------------------------------------

def parse_outdir(outdir: str | Path) -> dict:
    outdir = Path(outdir)
    rec: dict = {"outdir": str(outdir)}

    state_path = outdir / "state.json"
    if state_path.exists():
        try:
            rec["run_state"] = json.loads(state_path.read_text())
        except Exception:  # noqa: BLE001
            rec["run_state"] = None
    else:
        rec["run_state"] = None

    scf_txt = (outdir / "scf.out").read_text(errors="replace") if (outdir / "scf.out").exists() else ""
    ph_txt = (outdir / "ph.out").read_text(errors="replace") if (outdir / "ph.out").exists() else ""
    dyn_txt = (outdir / "dynmat.out").read_text(errors="replace") if (outdir / "dynmat.out").exists() else ""

    rec["scf"] = parse_scf(scf_txt)
    rec["ph_converged"] = "Convergence has been achieved" in ph_txt
    rec["ph_job_done"] = "JOB DONE" in ph_txt
    rec["dielectric"] = parse_dielectric(ph_txt)
    rec["asr_metric"] = parse_asr(dyn_txt)

    dyn = parse_dynmat(dyn_txt)
    rec["dynmat"] = {
        "cm_factor": dyn["cm_factor"],
        "raman_cm_factor": dyn["raman_cm_factor"],
        "polarizability": dyn["polarizability"],
        "has_raman": dyn["has_raman"],
    }
    rec["modes"] = dyn["modes"]
    rec["stats"] = _stats(dyn["modes"])

    lo: dict = {}
    for p in sorted(outdir.glob("dynmat_lo*.out")):
        key = p.stem[len("dynmat_"):]  # "lo", "lo_1", "lo_2", ...
        if not key:
            key = p.stem
        lod = parse_dynmat(p.read_text(errors="replace"))
        lo[key] = {"modes": lod["modes"], "has_raman": lod["has_raman"],
                   "stats": _stats(lod["modes"])}
    rec["lo"] = lo

    # status
    if not rec["ph_job_done"] or not dyn["modes"]:
        rec["status"] = "failed"
    elif rec["stats"]["soft_modes"]:
        rec["status"] = "warning"
    else:
        rec["status"] = "success"
    return rec


def write_outputs(rec: dict) -> None:
    outdir = Path(rec["outdir"])
    csv_path = outdir / "raman_modes.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["table", "mode", "freq_cm1", "freq_thz", "ir", "raman", "depol"])
        for m in rec["modes"]:
            w.writerow(["to", m["mode"], m["freq"], m["thz"], m["ir"], m["raman"],
                        "" if m["depol"] is None else m["depol"]])
        for key in sorted(rec["lo"]):
            for m in rec["lo"][key]["modes"]:
                w.writerow([key, m["mode"], m["freq"], m["thz"], m["ir"], m["raman"],
                            "" if m["depol"] is None else m["depol"]])

    summary = {k: v for k, v in rec.items() if k not in ("modes", "run_state")}
    summary["run_state"] = rec["run_state"]
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")


def print_summary(rec: dict) -> None:
    print("parse: %s" % rec["outdir"])
    s = rec["scf"]
    print("  scf: nelec=%s nbnd=%s nk=%s ecutwfc=%s sym_ops=%s iters=%s etot=%s Ry"
          % (s.get("nelec"), s.get("nbnd"), s.get("nk"), s.get("ecutwfc"),
             s.get("sym_ops"), s.get("scf_iter"), s.get("etot_ry")))
    print("  ph : converged=%s job_done=%s asr=%s"
          % (rec["ph_converged"], rec["ph_job_done"], rec["asr_metric"]))
    if rec["dielectric"]:
        d = rec["dielectric"]
        diag = [d[i][i] for i in range(3)]
        print("  dielectric diag: %s" % ", ".join("%.4f" % x for x in diag))
    st = rec["stats"]
    print("  modes: total=%d optical=%d raman_active=%d ir_active=%d sum_raman=%.4f"
          % (st["total"], len(st["optical"]), len(st["raman_active"]),
             len(st["ir_active"]), st["sum_raman"]))
    if st["max_raman"] is not None:
        print("  max raman: %.4f at %.2f cm-1" % (st["max_raman"], st["max_raman_freq"]))
    if st["soft_modes"]:
        print("  soft (imaginary) optical modes: %s cm-1"
              % ", ".join("%.2f" % m["freq"] for m in rec["modes"] if m["mode"] in st["soft_modes"]))
    if rec["lo"]:
        print("  LO tables: %s" % ", ".join(sorted(rec["lo"])))
    print("  status: %s" % rec["status"])


def main(argv: list[str] | None = None) -> int:
    import argparse

    p = argparse.ArgumentParser(prog="qe_parse.py", description="parse QE Raman outputs")
    p.add_argument("--outdir", required=True, help="run directory with scf.out/ph.out/dynmat.out")
    args = p.parse_args(argv)

    rec = parse_outdir(args.outdir)
    print_summary(rec)
    write_outputs(rec)
    print("  wrote raman_modes.csv + summary.json")
    return 0 if rec["status"] != "failed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
