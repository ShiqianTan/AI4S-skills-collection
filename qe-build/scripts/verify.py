#!/usr/bin/env python3
"""Run the Si diamond smoke test for a QE build and check the energy."""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SKILL = Path(__file__).resolve().parent.parent
EXPECTED_PBE_RY = -16.92302587


def _load_config(path: str | None):
    """Optionally use the shared qe_bench_core config resolver."""
    for parent in [Path(__file__).resolve(), *Path(__file__).resolve().parents]:
        if (parent / "_shared").is_dir():
            sys.path.insert(0, str(parent / "_shared"))
            try:
                from qe_bench_core import paths  # type: ignore

                return paths.load_config(path)
            except Exception:
                return None
    return None


def _replace_pseudo_dir(text: str, pseudo_dir: Path) -> str:
    if re.search(r"pseudo_dir\s*=", text):
        return re.sub(r"pseudo_dir\s*=\s*'[^']*'", f"pseudo_dir = '{pseudo_dir}'", text)
    return text.replace("&CONTROL", f"&CONTROL\n  pseudo_dir = '{pseudo_dir}'", 1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pw", default=None, help="path to pw.x")
    ap.add_argument("--mpirun", default=None, help="MPI launcher (default: mpirun)")
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--gpu", default=None, help="CUDA_VISIBLE_DEVICES value")
    ap.add_argument("--r2scan", action="store_true", help="use the R2SCAN input")
    ap.add_argument("--pseudo-dir", default=str(SKILL / "pseudo"))
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--config", default=None)
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    cfg = _load_config(args.config) or {}
    qe = cfg.get("qe", {}) if cfg else {}
    pw = args.pw or qe.get("libxc_bin_dir" if args.r2scan else "bin_dir")
    if pw and not str(pw).endswith("pw.x"):
        pw = str(Path(pw) / "pw.x")
    if not pw:
        ap.error("--pw is required (or configure qe.bin_dir in config.yaml)")

    mpirun = args.mpirun or qe.get("mpirun") or "mpirun"
    template = SKILL / "examples" / "si" / ("si_r2scan.in" if args.r2scan else "si.in")
    if not template.is_file():
        ap.error(f"missing template: {template}")

    workdir = Path(args.workdir) if args.workdir else Path(
        tempfile.mkdtemp(prefix="qe-build-verify-")
    )
    workdir.mkdir(parents=True, exist_ok=True)
    inp = workdir / "si.in"
    inp.write_text(_replace_pseudo_dir(template.read_text(), Path(args.pseudo_dir).resolve()))
    out = workdir / "si.out"

    env = os.environ.copy()
    env["OMPI_ALLOW_RUN_AS_ROOT"] = "1"
    env["OMPI_ALLOW_RUN_AS_ROOT_CONFIRM"] = "1"
    if args.gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    cmd = [mpirun, "-np", str(args.np), str(pw), "-in", inp.name]
    print("running:", " ".join(cmd))
    with out.open("w") as fo, (workdir / "si.err").open("w") as fe:
        proc = subprocess.run(cmd, cwd=workdir, env=env, stdout=fo, stderr=fe)

    text = out.read_text(errors="ignore")
    ok = proc.returncode == 0 and "JOB DONE" in text
    energies = [float(x) for x in re.findall(r"!\s+total energy\s+=\s+([-\d.]+)\s+Ry", text)]
    energy = energies[-1] if energies else None
    gpu_active = "GPU acceleration is ACTIVE" in text
    r2scan_active = "Exchange-correlation = R2SCAN" in text

    print(f"exit={proc.returncode} JOB_DONE={'JOB DONE' in text}")
    print(f"total_energy_Ry={energy}")
    if args.gpu is not None:
        print(f"GPU_banner={gpu_active}")
    if args.r2scan:
        print(f"R2SCAN_banner={r2scan_active}")

    if energy is not None and not args.r2scan:
        if abs(energy - EXPECTED_PBE_RY) <= 1e-6:
            print(f"PASS: energy matches {EXPECTED_PBE_RY} Ry")
        else:
            print(f"FAIL: expected {EXPECTED_PBE_RY} Ry, got {energy}")
            ok = False
    if args.r2scan and not r2scan_active:
        print("FAIL: R2SCAN functional banner missing")
        ok = False
    if args.gpu is not None and not gpu_active:
        print("FAIL: GPU acceleration banner missing")
        ok = False

    if not args.keep and args.workdir is None:
        shutil.rmtree(workdir, ignore_errors=True)
    else:
        print(f"artifacts: {workdir}")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
