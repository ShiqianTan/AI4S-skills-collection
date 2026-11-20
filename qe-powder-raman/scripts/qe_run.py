#!/usr/bin/env python3
"""qe_run — execute the QE powder-Raman workflow (pw.x -> ph.x -> dynmat.x).

Writes/updates ``state.json`` in the run directory with per-stage status,
exit code, elapsed time, the exact command, and the log file name.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import tempfile
import time
from pathlib import Path

from qe_parse import assess_electronic_structure, parse_scf

STATE_FILE = "state.json"
DOWNSTREAM = {
    "relax": ("scf", "ph", "dynmat", "lo"),
    "scf": ("ph", "dynmat", "lo"),
    "ph": ("dynmat", "lo"),
}


def _save_state(cwd: Path, state: dict) -> None:
    handle = tempfile.NamedTemporaryFile(mode="w", dir=cwd, prefix=".state-",
                                         suffix=".tmp", delete=False)
    tmp = Path(handle.name)
    try:
        with handle:
            json.dump(state, handle, indent=2, default=str)
            handle.write("\n")
        tmp.replace(cwd / STATE_FILE)
    finally:
        tmp.unlink(missing_ok=True)


def _artifacts(name: str, input_text: str, cwd: Path) -> dict | None:
    def setting(key):
        match = re.search(rf"\b{key}\s*=\s*(['\"])(.*?)\1", input_text, re.I)
        return match.group(2) if match else None

    files = []
    if name == "scf" and (prefix := setting("prefix")) and (scratch := setting("outdir")):
        save = cwd / scratch / (prefix + ".save")
        files = [save / "data-file-schema.xml"]
        # QE may be built with or without HDF5: accept either storage suffix.
        charge = (sorted(save.glob("charge-density.dat"))
                  + sorted(save.glob("charge-density.hdf5")))
        if not charge:
            return None
        files.append(charge[0])
        files.extend(sorted(set(save.glob("wfc*.dat"))
                            | set(save.glob("wfc*.hdf5"))))
    elif name == "ph" and (matrix := setting("fildyn")):
        files = [cwd / matrix]
    elif name.startswith("dynmat"):
        files = [cwd / value for key in ("filout", "filmol", "filxsf")
                 if (value := setting(key))]
    try:
        sizes = {str(path): path.stat().st_size for path in files if path.is_file()}
        if len(sizes) != len(files) or not all(sizes.values()):
            return None
    except OSError:
        return None
    if name == "ph" and files:
        return {str(files[0]): hashlib.sha256(files[0].read_bytes()).hexdigest()}
    return sizes


def _upstream_issue(stages: dict, name: str, cwd: Path) -> str | None:
    for upstream in ("scf", "ph"):
        if name not in DOWNSTREAM[upstream]:
            continue
        record = stages.get(upstream, {})
        if record.get("status") in ("stale", "running", "failed", "blocked"):
            return f"{upstream} is {record['status']}; rerun it first"
        signature = record.get("signature")
        if signature:
            try:
                text = (cwd / signature["command"][-1]).read_bytes()
                log = (cwd / record["log"]).read_bytes()
            except OSError:
                return f"{upstream} input or log is missing; rerun it first"
            artifacts = _artifacts(upstream, text.decode(errors="replace"), cwd)
            if (hashlib.sha256(text).hexdigest() != signature["input_sha256"]
                    or hashlib.sha256(log).hexdigest() != record.get("log_sha256")
                    or artifacts is None or artifacts != record.get("artifacts")):
                return f"{upstream} inputs or artifacts changed; rerun it first"
    return None


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
    records = stages.setdefault("lo", {}) if name.startswith("dynmat_lo") else stages
    rec = records.setdefault(name, {})
    input_bytes = (cwd / cmd[-1]).read_bytes()
    input_text = input_bytes.decode(errors="replace")
    signature = {
        "command": cmd,
        "input_sha256": hashlib.sha256(input_bytes).hexdigest(),
        "environment": {key: (env or {}).get(key) for key in
                        ("OMP_NUM_THREADS", "ESPRESSO_PSEUDO")},
    }
    if (skip_existing and rec.get("status") in ("success", "skipped")
            and rec.get("signature") == signature and log.is_file()
            and rec.get("log_sha256") == hashlib.sha256(log.read_bytes()).hexdigest()
            and (artifacts := _artifacts(name, input_text, cwd)) is not None
            and rec.get("artifacts") == artifacts):
        rec["status"] = "skipped"
        rec["reason"] = "matching inputs, log and artifacts"
        if verbose:
            print("  [skip] %s (already done: %s)" % (name, log.name))
        return True

    if verbose:
        print("  [run]  %-10s %s" % (name, " ".join(cmd)))

    for downstream in DOWNSTREAM.get(name, ()):
        dependent = stages.get(downstream)
        if dependent:
            for entry in dependent.values() if downstream == "lo" else (dependent,):
                entry.update(status="stale", reason=f"{name} restarted")
    rec.clear()
    rec.update(command=cmd, log=log.name, signature=signature, status="running")
    state["overall_status"] = "running"
    _save_state(cwd, state)
    start = time.time()
    try:
        with open(log, "w") as fh:
            proc = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env)
        rec["exit_code"] = proc.returncode
        rec["seconds"] = round(time.time() - start, 3)
        rec["status"] = "failed"
        output = log.read_bytes()
        if proc.returncode == 0:
            if b"JOB DONE." in output:
                rec.update(status="success", log_sha256=hashlib.sha256(output).hexdigest(),
                           artifacts=_artifacts(name, input_text, cwd))
            else:
                rec["error"] = "QE completion marker missing"
        elif name == "ph" and b"set_irr_sym_new" in output and b"wrong representation" in output:
            rec["error"] = ("QE phonon symmetry representation failed; retry in a fresh "
                            "directory with ph_search_sym: false if the structure is valid")
    except Exception as exc:  # noqa: BLE001
        rec["status"] = "failed"
        rec["exit_code"] = None
        rec["seconds"] = round(time.time() - start, 3)
        rec["error"] = str(exc)

    _save_state(cwd, state)
    if rec["status"] == "failed":
        if verbose:
            print("  [FAIL] %s (exit %s) -> %s%s" % (
                name, rec.get("exit_code"), log.name,
                ": " + rec["error"] if rec.get("error") else ""))
        return False
    if verbose:
        print("  [done] %-10s %.1fs -> %s" % (name, rec.get("seconds", 0.0), log.name))
    return True


def _regenerate_from_relax(outdir: Path, pseudo_dir: Path, bin_dir: Path,
                           verbose: bool) -> str | None:
    """After a successful relax, rebuild scf/ph/dynmat inputs from the relaxed
    structure and record residual force/stress/symmetry in relax_manifest.json."""
    import json

    from qe_struct import parse_final_structure, relax_metrics
    from qe_inputs import generate_inputs

    relax_out = outdir / "relax.out"
    if not relax_out.exists():
        return "relax.out missing"
    text = relax_out.read_text(errors="replace")
    if "maximum number of steps has been reached" in text.lower():
        return "vc-relax reached the maximum number of steps before convergence"
    if not re.search(r"\bbfgs converged in\b", text, re.I):
        return "vc-relax has no BFGS convergence marker"
    try:
        atoms = parse_final_structure(text)
    except (KeyError, ValueError):
        atoms = None
    if atoms is None:
        return "could not parse final coordinates from relax.out"
    initial = (outdir / "relax.in").read_text(errors="replace")
    nat = re.search(r"\bnat\s*=\s*(\d+)", initial, re.I)
    if nat and len(atoms) != int(nat.group(1)):
        return f"vc-relax returned {len(atoms)} atoms, expected {nat.group(1)}"

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
    return None


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
        issue = _upstream_issue(state["stages"], name, outdir)
        if issue:
            print(f"  [FAIL] {name}: {issue}")
            ok = False
            break
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
            error = _regenerate_from_relax(outdir, pseudo_dir, bin_dir, verbose)
            if error:
                state["stages"]["relax"].update(status="failed", error=error)
                print(f"  [FAIL] relax: {error}")
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
            for loin in lo_inputs:
                key = loin.stem  # e.g. dynmat_lo_1
                tag = key[len("dynmat_"):] if key.startswith("dynmat_") else key
                loout = outdir / (loin.stem + ".out")
                cmd = [str(dyn), "-in", loin.name]
                if not _run_stage(state, key, cmd, loout, outdir, skip_existing, verbose, env):
                    ok = False
                    continue
                # Older inputs used QE's default filenames for both TO and LO.
                if state["stages"]["lo"][key]["status"] == "success":
                    lo_input = loin.read_text(errors="replace")
                    for variable, ext in (("filmol", ".mold"), ("filxsf", ".axsf")):
                        src = outdir / ("dynmat" + ext)
                        if not re.search(rf"\b{variable}\s*=", lo_input, re.I) and src.is_file():
                            (outdir / ("dynmat_%s%s" % (tag, ext))).write_bytes(src.read_bytes())

    state["overall_status"] = "success" if ok else "failed"
    _save_state(outdir, state)
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
    p.add_argument("--skip-existing", action="store_true", help="reuse completed stages with matching inputs and outputs")
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
