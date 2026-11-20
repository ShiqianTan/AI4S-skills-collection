#!/usr/bin/env python3
"""qe_run — execute the QE powder-Raman workflow (pw.x -> ph.x -> dynmat.x).

Writes/updates ``state.json`` in the run directory with per-stage status,
exit code, elapsed time, the exact command, and the log file name.
"""
from __future__ import annotations

import json
import re
import subprocess
import time
from pathlib import Path

from qe_parse import assess_electronic_structure, parse_scf

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


def _regenerate_from_relax(outdir: Path, pseudo_dir: Path, bin_dir: Path,
                           verbose: bool) -> bool:
    """After a successful relax, rebuild scf/ph/dynmat inputs from the relaxed
    structure and record residual force/stress/symmetry in relax_manifest.json."""
    import json

    from qe_struct import parse_final_structure, relax_metrics
    from qe_inputs import generate_inputs

    relax_out = outdir / "relax.out"
    if not relax_out.exists():
        return False
    text = relax_out.read_text(errors="replace")
    atoms = parse_final_structure(text)
    if atoms is None:
        if verbose:
            print("  [FAIL] could not parse final coordinates from relax.out")
        return False

    eff = outdir / "effective_config.json"
    cfg = json.loads(eff.read_text()) if eff.exists() else {}
    im = outdir / "inputs_manifest.json"
    man = json.loads(im.read_text()) if im.exists() else {}
    prefix = man.get("prefix") or cfg.get("prefix", "raman")
    lo_dirs = [tuple(q) for q in (man.get("lo_dirs") or [])]

    sg = None
    try:
        import spglib

        ds = spglib.get_symmetry_dataset(
            (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers()),
            symprec=1e-3)
        if ds is not None:
            sg = {"number": ds.number, "symbol": ds.international}
    except Exception:  # noqa: BLE001
        sg = None

    metrics = relax_metrics(text)
    metrics.update({"n_atoms": len(atoms), "spacegroup": sg})
    (outdir / "relax_manifest.json").write_text(json.dumps(metrics, indent=2) + "\n")
    if verbose:
        print("  [relax] inputs regenerated from relaxed structure "
              "(%d atoms, residual force=%s)" % (len(atoms), metrics.get("final_force")))

    generate_inputs(atoms, cfg, outdir, pseudo_dir, prefix, lo_dirs=lo_dirs,
                    with_relax=False, bin_dir=bin_dir)
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
    with_relax: bool = False,
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
        if with_relax:
            stages_to_run = ["relax"] + stages_to_run
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
        if name == "relax":
            if not (outdir / "relax.in").exists():
                print("  [FAIL] relax.in not found in %s (regenerate with --relax)" % outdir)
                ok = False
                break
            cmd = launcher(pw, ["-in", "relax.in"])
            if not _run_stage(state, "relax", cmd, outdir / "relax.out", outdir,
                              skip_existing, verbose, env):
                ok = False
                break
            if not _regenerate_from_relax(outdir, pseudo_dir, bin_dir, verbose):
                ok = False
                break
        elif name == "scf":
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
            scf_log = outdir / "scf.out"
            if scf_log.is_file():
                electronic = assess_electronic_structure(
                    parse_scf(scf_log.read_text(errors="replace"))
                )
                state["electronic_structure"] = electronic
                if electronic["classification"] == "metallic":
                    state["stages"]["ph"] = {"status": "blocked", "reason": electronic["reason"]}
                    print("  [FAIL] ph: %s" % electronic["reason"])
                    ok = False
                    break
                if verbose and electronic["classification"] == "undeterminable":
                    print("  [WARN] ph: %s" % electronic["reason"])
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
                # Older inputs used QE's default filenames for both TO and LO.
                if rec["status"] == "success":
                    lo_input = loin.read_text(errors="replace")
                    for variable, ext in (("filmol", ".mold"), ("filxsf", ".axsf")):
                        src = outdir / ("dynmat" + ext)
                        if not re.search(rf"\b{variable}\s*=", lo_input, re.I) and src.is_file():
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
