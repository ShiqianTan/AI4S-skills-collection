#!/usr/bin/env python3
"""qe_run — execute the QE powder-Raman workflow (pw.x -> ph.x -> dynmat.x).

Writes/updates ``state.json`` in the run directory with per-stage status,
exit code, elapsed time, the exact command, and the log file name.
"""
from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

STATE_FILE = "state.json"


def _launcher(np: int, mpi_command: str, exe: Path, args: list[str]) -> list[str]:
    """Return the full argv for a QE executable (with MPI prefix when np > 1)."""
    if np and int(np) > 1 and mpi_command:
        return [mpi_command, "-np", str(int(np)), str(exe), *args]
    return [str(exe), *args]


def _run_stage(
    state: dict,
    name: str,
    cmd: list[str],
    log: Path,
    cwd: Path,
    skip_existing: bool,
    verbose: bool,
    env: dict | None = None,
) -> bool:
    stages = state["stages"]
    if name not in stages:
        stages[name] = {}
    rec = stages[name]

    if skip_existing and log.exists() and rec.get("status") in ("success", "skipped"):
        rec["status"] = "skipped"
        rec["reason"] = "log exists"
        if verbose:
            print("  [skip] %s (already done: %s)" % (name, log.name))
        return True

    if verbose:
        print("  [run]  %-10s %s" % (name, " ".join(cmd)))

    rec["command"] = cmd
    rec["log"] = log.name
    start = time.time()
    try:
        with open(log, "w") as fh:
            proc = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env)
        rec["exit_code"] = proc.returncode
        rec["seconds"] = round(time.time() - start, 3)
        rec["status"] = "success" if proc.returncode == 0 else "failed"
    except Exception as exc:  # noqa: BLE001
        rec["status"] = "failed"
        rec["exit_code"] = None
        rec["seconds"] = round(time.time() - start, 3)
        rec["error"] = str(exc)

    if rec["status"] == "failed":
        if verbose:
            print("  [FAIL] %s (exit %s) -> %s" % (name, rec.get("exit_code"), log.name))
        return False
    if verbose:
        print("  [done] %-10s %.1fs -> %s" % (name, rec.get("seconds", 0.0), log.name))
    return True


def run_workflow(
    outdir: str | Path,
    bin_dir: str | Path,
    pseudo_dir: str | Path,
    np: int = 4,
    mpi_command: str = "mpirun",
    omp_num_threads: int = 1,
    stage: str = "all",
    skip_existing: bool = False,
    verbose: bool = True,
) -> tuple[dict, bool]:
    """Run the workflow stages in an existing prepared outdir.

    Returns ``(state, ok)``. ``state`` is persisted to ``<outdir>/state.json``.
    """
    outdir = Path(outdir)
    bin_dir = Path(bin_dir)
    pseudo_dir = Path(pseudo_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    state_path = outdir / STATE_FILE
    if state_path.exists():
        try:
            state = json.loads(state_path.read_text())
        except Exception:  # noqa: BLE001
            state = {}
    else:
        state = {}

    state.update(
        {
            "outdir": str(outdir),
            "bin_dir": str(bin_dir),
            "pseudo_dir": str(pseudo_dir),
            "np": int(np),
            "mpi_command": mpi_command,
            "omp_num_threads": int(omp_num_threads),
        }
    )
    state.setdefault("stages", {})

    pw = bin_dir / "pw.x"
    ph = bin_dir / "ph.x"
    dyn = bin_dir / "dynmat.x"

    def launcher(exe: Path, args: list[str]) -> list[str]:
        return _launcher(np, mpi_command, exe, args)

    stages_to_run: list[str] = []
    if stage == "all":
        stages_to_run = ["scf", "ph", "dynmat", "lo"]
    else:
        stages_to_run = [stage]

    env_run = {}
    if omp_num_threads is not None:
        env_run["OMP_NUM_THREADS"] = str(int(omp_num_threads))
    env_run["ESPRESSO_PSEUDO"] = str(pseudo_dir)
    # subprocess inherits os.environ plus these overrides.
    import os

    env = dict(os.environ)
    env.update(env_run)

    if verbose:
        print("run: %s" % outdir)
        print("  binaries=%s pseudo=%s np=%s" % (bin_dir, pseudo_dir, np))

    ok = True
    for name in stages_to_run:
        if name == "scf":
            if not (outdir / "scf.in").exists():
                print("  [FAIL] scf.in not found in %s" % outdir)
                ok = False
                break
            cmd = launcher(pw, ["-in", "scf.in"])
            if not _run_stage(state, "scf", cmd, outdir / "scf.out", outdir, skip_existing, verbose, env):
                ok = False
                break
        elif name == "ph":
            if not (outdir / "ph.in").exists():
                print("  [FAIL] ph.in not found in %s" % outdir)
                ok = False
                break
            cmd = launcher(ph, ["-in", "ph.in"])
            if not _run_stage(state, "ph", cmd, outdir / "ph.out", outdir, skip_existing, verbose, env):
                ok = False
                break
        elif name == "dynmat":
            if not (outdir / "dynmat.in").exists():
                print("  [FAIL] dynmat.in not found in %s" % outdir)
                ok = False
                break
            cmd = [str(dyn), "-in", "dynmat.in"]
            if not _run_stage(state, "dynmat", cmd, outdir / "dynmat.out", outdir, skip_existing, verbose, env):
                ok = False
                break
        elif name == "lo":
            lo_inputs = sorted(outdir.glob("dynmat_lo*.in"))
            if not lo_inputs:
                if verbose:
                    print("  [skip] no LO inputs (dynmat_lo*.in) in %s" % outdir)
                continue
            lo_stages = state["stages"].setdefault("lo", {})
            for loin in lo_inputs:
                key = loin.stem  # e.g. dynmat_lo_1
                tag = key[len("dynmat_"):] if key.startswith("dynmat_") else key
                loout = outdir / (loin.stem + ".out")
                cmd = [str(dyn), "-in", loin.name]
                # reuse _run_stage but with nested dict
                rec = lo_stages.setdefault(key, {})
                if skip_existing and loout.exists() and rec.get("status") in ("success", "skipped"):
                    rec["status"] = "skipped"
                    rec["reason"] = "log exists"
                    if verbose:
                        print("  [skip] %s (already done)" % key)
                    continue
                if verbose:
                    print("  [run]  %-10s %s" % (key, " ".join(cmd)))
                rec["command"] = cmd
                rec["log"] = loout.name
                start = time.time()
                try:
                    with open(loout, "w") as fh:
                        proc = subprocess.run(cmd, cwd=outdir, stdout=fh, stderr=subprocess.STDOUT, env=env)
                    rec["exit_code"] = proc.returncode
                    rec["seconds"] = round(time.time() - start, 3)
                    rec["status"] = "success" if proc.returncode == 0 else "failed"
                except Exception as exc:  # noqa: BLE001
                    rec["status"] = "failed"
                    rec["error"] = str(exc)
                if rec["status"] == "failed":
                    if verbose:
                        print("  [FAIL] %s -> %s" % (key, loout.name))
                    ok = False
                elif verbose:
                    print("  [done] %-10s %.1fs -> %s" % (key, rec.get("seconds", 0.0), loout.name))
                # keep per-direction eigenvector artifacts (dynmat.mold/axsf are
                # overwritten by every LO run)
                for src_name, ext in (("dynmat.mold", ".mold"), ("dynmat.axsf", ".axsf")):
                    src = outdir / src_name
                    if src.exists():
                        (outdir / ("dynmat_%s%s" % (tag, ext))).write_bytes(src.read_bytes())

    state["overall_status"] = "success" if ok else "failed"
    tmp_path = state_path.with_suffix(".json.tmp")
    tmp_path.write_text(json.dumps(state, indent=2, default=str) + "\n")
    tmp_path.replace(state_path)
    return state, ok


def main(argv: list[str] | None = None) -> int:
    import argparse

    p = argparse.ArgumentParser(prog="qe_run.py", description="run QE powder-Raman workflow")
    p.add_argument("--outdir", required=True, help="prepared run directory")
    p.add_argument("--bin-dir", default=None, help="QE binary directory")
    p.add_argument("--pseudo-dir", default=None, help="pseudopotential directory")
    p.add_argument("--np", type=int, default=4, help="MPI processes")
    p.add_argument("--mpi-command", default="mpirun", help="MPI launcher")
    p.add_argument("--omp-num-threads", type=int, default=1, help="OMP_NUM_THREADS")
    p.add_argument("--stage", default="all", choices=["all", "scf", "ph", "dynmat", "lo"])
    p.add_argument("--skip-existing", action="store_true")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)

    # Resolve defaults the same way the driver does.
    from qe_common import resolve_bin_dir, resolve_pseudo_dir

    bin_dir = resolve_bin_dir(args.bin_dir)
    pseudo_dir = resolve_pseudo_dir(args.pseudo_dir)

    _, ok = run_workflow(
        args.outdir, bin_dir, pseudo_dir,
        np=args.np, mpi_command=args.mpi_command, omp_num_threads=args.omp_num_threads,
        stage=args.stage, skip_existing=args.skip_existing, verbose=not args.quiet,
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
