#!/usr/bin/env python3
"""
qe-powder-raman — command-line driver for the Quantum ESPRESSO powder-Raman
workflow (pw.x SCF -> ph.x epsil+lraman at Gamma -> dynmat.x -> spectrum).

Subcommands (implemented incrementally, one step at a time):

    preflight      check environment + structure + pseudopotentials  [Step 0]
    prepare        read structure, map pseudopotentials, write inputs [Step 1-2]
    run            execute pw.x / ph.x / dynmat.x                     [Step 3]
    parse          parse outputs into raman_modes.csv / summary.json  [Step 3]
    spectrum       broaden modes into a powder-Raman spectrum          [Step 4]
    report         assemble report.md                                 [Step 4]
    all            prepare [-> run ->] parse -> spectrum -> report  (Step 5)
    validate-input validate a generated QE input file
    list-templates list available input templates

All subcommands implemented (Steps 0-5).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# Allow running as `python scripts/qe_raman.py` from anywhere.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qe_common import (  # noqa: E402
    DEFAULT_BIN_DIR,
    QE_BINARIES,
    PSEUDO_DIR,
    load_defaults,
    resolve_bin_dir,
    resolve_pseudo_dir,
    validate_config,
)
from qe_pp import check_pseudo, element_pp_map, build_pseudopotential_manifest  # noqa: E402
from qe_struct import (  # noqa: E402
    check_structure,
    format_check_result,
    read_structure,
    refine_structure,
    write_validated_cif,
)
from qe_inputs import generate_inputs  # noqa: E402
from qe_run import run_workflow  # noqa: E402
from qe_parse import parse_outdir, print_summary, write_outputs  # noqa: E402
from qe_spectrum import generate_spectra, load_modes, select_peaks  # noqa: E402
from qe_report import render_report, write_report  # noqa: E402

_NOT_IMPLEMENTED = (
    "`{cmd}` is not implemented yet (arrives in a later step of the skill "
    "build). The available command is `preflight`."
)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="qe_raman.py",
        description="QE powder-Raman workflow driver",
    )
    sub = p.add_subparsers(dest="command", required=True)

    pf = sub.add_parser("preflight", help="check env + structure + pseudos")
    pf.add_argument("structure", help="CIF / POSCAR / PW input file")
    pf.add_argument("--bin-dir", default=None, help="QE binary directory")
    pf.add_argument("--pseudo-dir", default=None, help="pseudopotential directory")
    pf.add_argument("--skip-structure", action="store_true", help="skip structure checks")

    vi = sub.add_parser("validate-input", help="validate a generated QE input")
    vi.add_argument("file", help="QE input file to validate")

    lt = sub.add_parser("list-templates", help="list input templates")

    st = sub.add_parser(
        "structure", help="standardize structure + map pseudopotentials [Step 1]"
    )
    st.add_argument("structure", help="CIF / POSCAR / PW input file")
    st.add_argument("--outdir", default="runs", help="output directory (validated.cif + manifest)")
    st.add_argument("--pseudo-dir", default=None, help="pseudopotential directory")
    st.add_argument("--symprec", type=float, default=1e-3, help="spglib symmetry tolerance")
    st.add_argument("--no-refine", action="store_true", help="skip spglib.refine_cell")

    pr = sub.add_parser(
        "prepare", help="standardize structure + map pseudos + write inputs [Step 2]"
    )
    pr.add_argument("structure", help="CIF / POSCAR / PW input file")
    pr.add_argument("--config", default=None, help="YAML config file")
    pr.add_argument("--outdir", default="runs", help="output directory for inputs")
    pr.add_argument("--pseudo-dir", default=None, help="pseudopotential directory")
    pr.add_argument("--prefix", default=None, help="QE prefix (default from config)")
    pr.add_argument("--symprec", type=float, default=1e-3, help="spglib symmetry tolerance")
    pr.add_argument("--no-refine", action="store_true", help="skip spglib.refine_cell")
    pr.add_argument("--relax", action="store_true", help="also write relax.in (vc-relax)")
    pr.add_argument("--lo", nargs="*", type=float, default=None,
                    help="LO directions 'x y z [x y z ...]'")

    rn = sub.add_parser("run", help="execute pw.x/ph.x/dynmat.x in a prepared outdir [Step 3]")
    rn.add_argument("--outdir", default=None, help="prepared run directory")
    rn.add_argument("--config", default=None, help="YAML config file")
    rn.add_argument("--bin-dir", default=None, help="QE binary directory")
    rn.add_argument("--pseudo-dir", default=None, help="pseudopotential directory")
    rn.add_argument("--np", type=int, default=None, help="MPI processes")
    rn.add_argument("--mpi-command", default=None, help="MPI launcher")
    rn.add_argument("--omp-num-threads", type=int, default=None)
    rn.add_argument("--stage", default="all", choices=["all", "scf", "ph", "dynmat", "lo"])
    rn.add_argument("--skip-existing", action="store_true", help="skip stages whose log exists")

    ps = sub.add_parser("parse", help="parse outputs into raman_modes.csv + summary.json [Step 3]")
    ps.add_argument("--outdir", default=None, help="run directory")

    spc = sub.add_parser("spectrum", help="broaden modes into a powder-Raman spectrum [Step 4]")
    spc.add_argument("--outdir", default=None, help="run directory with raman_modes.csv")
    spc.add_argument("--config", default=None, help="YAML config file")
    spc.add_argument("--type", default=None, choices=["activity", "stokes", "both"],
                     help="spectrum type (default activity)")
    spc.add_argument("--broadening", default=None, choices=["lorentzian", "gaussian"])
    spc.add_argument("--fwhm", type=float, default=None, help="broadening FWHM (cm-1)")
    spc.add_argument("--temperature", type=float, default=None, help="Kelvin (Stokes)")
    spc.add_argument("--laser-wavelength", type=float, default=None, help="nm (Stokes)")
    spc.add_argument("--threshold", type=float, default=None, help="Raman activity cutoff")
    spc.add_argument("--formats", default=None, help="comma list: png,svg,csv")

    rpt = sub.add_parser("report", help="assemble report.md [Step 4]")
    rpt.add_argument("--outdir", default=None, help="run directory with summary.json")
    rpt.add_argument("--config", default=None, help="YAML config file")

    al = sub.add_parser(
        "all", help="prepare [-> run ->] parse -> spectrum -> report [Step 5]"
    )
    al.add_argument("structure", help="CIF / POSCAR / PW input file")
    al.add_argument("--config", default=None, help="YAML config file")
    al.add_argument("--outdir", default=None, help="run directory")
    al.add_argument("--pseudo-dir", default=None, help="pseudopotential directory")
    al.add_argument("--prefix", default=None, help="QE prefix")
    al.add_argument("--symprec", type=float, default=1e-3, help="spglib symmetry tolerance")
    al.add_argument("--no-refine", action="store_true", help="skip spglib.refine_cell")
    al.add_argument("--relax", action="store_true", help="also write relax.in (vc-relax)")
    al.add_argument("--bin-dir", default=None, help="QE binary directory")
    al.add_argument("--np", type=int, default=None, help="MPI processes")
    al.add_argument("--mpi-command", default=None, help="MPI launcher")
    al.add_argument("--omp-num-threads", type=int, default=None)
    al.add_argument("--stage", default="all", choices=["all", "scf", "ph", "dynmat", "lo"])
    al.add_argument("--skip-existing", action="store_true", help="skip stages whose log exists")
    al.add_argument("--execute", action="store_true",
                    help="actually run pw.x/ph.x/dynmat.x (default: generate only)")
    al.add_argument("--type", default=None, choices=["activity", "stokes", "both"],
                    help="spectrum type")
    al.add_argument("--broadening", default=None, choices=["lorentzian", "gaussian"])
    al.add_argument("--fwhm", type=float, default=None)
    al.add_argument("--temperature", type=float, default=None)
    al.add_argument("--laser-wavelength", type=float, default=None)
    al.add_argument("--threshold", type=float, default=None)
    al.add_argument("--formats", default=None, help="comma list: png,svg,csv")
    return p


# ---------------------------------------------------------------------------
# preflight
# ---------------------------------------------------------------------------


def _check_binaries(bin_dir: Path) -> list[tuple[str, str, str]]:
    msgs: list[tuple[str, str, str]] = []
    if not bin_dir.is_dir():
        return [("FAIL", "bin_dir", "QE binary directory not found: %s" % bin_dir)]
    version = None
    for name in QE_BINARIES:
        exe = bin_dir / name
        if exe.exists() and (exe.stat().st_mode & 0o111):
            msgs.append(("PASS", name, str(exe)))
        else:
            msgs.append(("FAIL", name, "missing/not executable: %s" % exe))
    pw = bin_dir / "pw.x"
    if pw.exists():
        try:
            out = subprocess.run(
                [str(pw), "--version"], capture_output=True, text=True, timeout=20
            )
            first = (out.stdout or out.stderr).splitlines()
            version = next((ln.strip() for ln in first if "Program PWSCF" in ln), None)
        except Exception:  # noqa: BLE001
            version = None
    if version:
        msgs.append(("PASS", "qe_version", version))
    else:
        msgs.append(("WARN", "qe_version", "could not read QE version from pw.x --version"))
    return msgs


def _check_mpi() -> list[tuple[str, str, str]]:
    mpirun = shutil.which("mpirun") or shutil.which("mpiexec")
    if mpirun:
        return [("PASS", "mpi", "mpirun/mpiexec found: %s" % mpirun)]
    return [("WARN", "mpi", "no mpirun/mpiexec on PATH; single-process runs only")]


def _check_python_deps() -> list[tuple[str, str, str]]:
    msgs: list[tuple[str, str, str]] = []
    required = {"numpy": "numpy", "yaml": "PyYAML", "ase": "ASE", "spglib": "spglib"}
    optional = {"matplotlib": "matplotlib"}
    for mod, label in required.items():
        try:
            __import__(mod)
            msgs.append(("PASS", "dep", "%s available" % label))
        except ImportError:
            msgs.append(("FAIL", "dep", "%s missing (pip install %s)" % (label, label.lower())))
    for mod, label in optional.items():
        try:
            __import__(mod)
            msgs.append(("PASS", "dep", "%s available (spectrum)" % label))
        except ImportError:
            msgs.append(("WARN", "dep", "%s missing (spectrum/report disabled)" % label))
    return msgs


def _check_pseudos(pseudo_dir: Path, elements: list[str]) -> list[tuple[str, str, str]]:
    msgs: list[tuple[str, str, str]] = []
    if not pseudo_dir.is_dir():
        return [("FAIL", "pseudo_dir", "pseudopotential directory not found: %s" % pseudo_dir)]
    mapping = element_pp_map(pseudo_dir, elements)
    for el in sorted(set(elements)):
        cands = mapping.get(el) or []
        if not cands:
            msgs.append(("FAIL", "pseudo", "no UPF for element %s in %s" % (el, pseudo_dir)))
        elif len(cands) > 1:
            msgs.append(("FAIL", "pseudo", "ambiguous UPF for %s: %d candidates (%s)"
                         % (el, len(cands), ", ".join(cands))))
        else:
            msgs.extend(check_pseudo(cands[0], expected_element=el))
    return msgs


def _run_preflight(args: argparse.Namespace) -> int:
    bin_dir = resolve_bin_dir(args.bin_dir)
    pseudo_dir = resolve_pseudo_dir(args.pseudo_dir)

    all_msgs: list[tuple[str, str, str]] = []
    all_msgs.append(("INFO", "section", "--- environment ---"))
    all_msgs.extend(_check_binaries(bin_dir))
    all_msgs.extend(_check_mpi())
    all_msgs.extend(_check_python_deps())

    elements: list[str] = []
    if not args.skip_structure:
        res = check_structure(args.structure)
        all_msgs.append(("INFO", "section", "--- structure ---"))
        all_msgs.append(
            ("INFO", "summary", "atoms=%d elements=%s" % (res.n_atoms, ",".join(res.elements)))
        )
        if res.angles is not None:
            all_msgs.append(
                ("INFO", "angles", "cell angles (deg): %s"
                 % ", ".join("%.6f" % a for a in res.angles))
            )
        all_msgs.extend(res.messages)
        elements = res.elements

    all_msgs.append(("INFO", "section", "--- pseudopotentials ---"))
    all_msgs.extend(_check_pseudos(pseudo_dir, elements))

    n_fail = 0
    n_warn = 0
    for level, key, text in all_msgs:
        print("[%s] %-14s %s" % (level, key, text))
        if level == "FAIL":
            n_fail += 1
        elif level == "WARN":
            n_warn += 1

    print()
    if n_fail:
        print("preflight: FAIL (%d failure(s), %d warning(s))" % (n_fail, n_warn))
        return 1
    print("preflight: OK (%d warning(s))" % n_warn)
    return 0


def _run_structure(args: argparse.Namespace) -> int:
    pseudo_dir = resolve_pseudo_dir(args.pseudo_dir)
    outdir = Path(args.outdir)

    atoms = read_structure(args.structure)
    res = check_structure(args.structure)
    print("--- structure checks ---")
    print(format_check_result(res))
    if res.has_fail:
        print("\nstructure: FAIL (fix the issues above before continuing)")
        return 1

    print("\n--- symmetry refinement ---")
    if args.no_refine:
        refined, info = atoms, {
            "standardized": False,
            "original": {
                "spacegroup_number": res.spacegroup_number,
                "spacegroup_symbol": res.spacegroup_symbol,
                "n_sym_ops": res.n_sym_ops,
            },
            "refined": {
                "spacegroup_number": res.spacegroup_number,
                "spacegroup_symbol": res.spacegroup_symbol,
                "n_sym_ops": res.n_sym_ops,
            },
            "cellpar": [round(float(x), 8) for x in atoms.cell.cellpar()],
            "n_atoms": len(atoms),
            "elements": list(atoms.get_chemical_symbols()),
        }
    else:
        refined, info = refine_structure(atoms, symprec=args.symprec)

    o = info["original"]
    r = info["refined"]
    print(
        "original:  sg #%s %s (%s ops)"
        % (o.get("spacegroup_number"), o.get("spacegroup_symbol"), o.get("n_sym_ops"))
    )
    print(
        "refined :  sg #%s %s (%s ops)"
        % (r.get("spacegroup_number"), r.get("spacegroup_symbol"), r.get("n_sym_ops"))
    )
    print("cellpar (a,b,c,alpha,beta,gamma): %s" % ", ".join("%.6f" % x for x in info["cellpar"]))
    if info["standardized"]:
        print("standardized: cell was numerically cleaned (angle/length noise removed)")
    else:
        print("standardized: no numerical cleanup needed")

    cif_path = write_validated_cif(refined, outdir)
    print("\nwrote %s" % cif_path)

    manifest = build_pseudopotential_manifest(info["elements"], pseudo_dir, outdir)
    print("wrote %s" % (outdir / "pseudopotential_manifest.json"))
    print("\n--- pseudopotential manifest ---")
    for el in manifest["elements"]:
        e = manifest["pseudopotentials"].get(el)
        if e:
            print(
                "  %-3s %-16s type=%s xc=%s z=%s"
                % (el, e["file"], e["pseudo_type"], e["functional"], e["z_valence"])
            )
    for w in manifest["warnings"]:
        print("  [WARN] %s" % w)
    for err in manifest["errors"]:
        print("  [FAIL] %s" % err)
    print("  total valence electrons = %s" % manifest["n_valence_electrons"])

    if manifest["errors"]:
        print("\nstructure: FAIL (pseudopotential problems)")
        return 1
    print("\nstructure: OK")
    return 0


# ---------------------------------------------------------------------------
# stubs
# ---------------------------------------------------------------------------


def _load_config(path: str | None) -> dict:
    cfg = load_defaults()
    if path:
        import yaml

        user = yaml.safe_load(Path(path).read_text()) or {}
        cfg.update({k: v for k, v in user.items()})
    return cfg


def _run_prepare(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    for level, key, text in validate_config(cfg):
        print("[%s] %-14s %s" % (level, key, text))
        if level == "FAIL":
            print("\nprepare: FAIL (invalid config; see above)")
            return 1
    pseudo_dir = resolve_pseudo_dir(args.pseudo_dir or cfg.get("pseudo_dir"))
    bin_dir = resolve_bin_dir(getattr(args, "bin_dir", None) or cfg.get("bin_dir"))
    outdir = Path(args.outdir or cfg.get("outdir", "runs"))
    prefix = args.prefix or cfg.get("prefix", "raman")
    lo_vals = getattr(args, "lo", None)
    if lo_vals is not None:
        lo_dirs = [tuple(lo_vals[i:i + 3]) for i in range(0, len(lo_vals), 3)]
    else:
        lo_dirs = [tuple(q) for q in (cfg.get("lo_to_directions") or [])]

    atoms = read_structure(args.structure)
    res = check_structure(args.structure)
    print("--- structure checks ---")
    print(format_check_result(res))
    if res.has_fail:
        print("\nprepare: FAIL (structure issues; fix before continuing)")
        return 1

    print("\n--- symmetry refinement ---")
    if args.no_refine:
        refined, info = atoms, {
            "standardized": False,
            "cellpar": [round(float(x), 8) for x in atoms.cell.cellpar()],
            "elements": list(atoms.get_chemical_symbols()),
        }
    else:
        refined, info = refine_structure(atoms, symprec=args.symprec)
    print("cellpar (a,b,c,alpha,beta,gamma): %s"
          % ", ".join("%.6f" % x for x in info["cellpar"]))
    print("standardized: %s" % ("cleaned" if info["standardized"] else "no cleanup needed"))

    cif_path = write_validated_cif(refined, outdir)
    print("wrote %s" % cif_path)

    manifest = build_pseudopotential_manifest(info["elements"], pseudo_dir, outdir)
    if manifest["errors"]:
        for err in manifest["errors"]:
            print("  [FAIL] %s" % err)
        print("\nprepare: FAIL (pseudopotential problems)")
        return 1

    print("\n--- inputs ---")
    summary = generate_inputs(
        refined, cfg, outdir, pseudo_dir, prefix,
        lo_dirs=lo_dirs, with_relax=args.relax, bin_dir=bin_dir,
    )
    for f in summary["files"]:
        print("  " + f)
    print("kpoints %s, species %s, lo %s"
          % (summary["kpoints"], ",".join(summary["species"]), summary["lo_dirs"]))
    print("\nprepare: OK -> %s" % outdir)
    return 0


def _stub(command: str) -> int:
    print(_NOT_IMPLEMENTED.format(cmd=command), file=sys.stderr)
    return 2


def _run_run(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    for level, key, text in validate_config(cfg):
        print("[%s] %-14s %s" % (level, key, text))
        if level == "FAIL":
            print("\nrun: FAIL (invalid config; see above)")
            return 1
    outdir = Path(args.outdir or cfg.get("outdir", "runs"))
    if not outdir.is_dir() or not (outdir / "scf.in").exists():
        print("run: FAIL (no prepared inputs in %s; run `prepare` first)" % outdir,
              file=sys.stderr)
        return 1

    bin_dir = resolve_bin_dir(args.bin_dir or cfg.get("bin_dir"))
    pseudo_dir = resolve_pseudo_dir(args.pseudo_dir or cfg.get("pseudo_dir"))
    np = args.np if args.np is not None else int(cfg.get("np", 4))
    mpi_command = args.mpi_command or cfg.get("mpi_command", "mpirun")
    omp = args.omp_num_threads if args.omp_num_threads is not None else int(cfg.get("omp_num_threads", 1))

    _, ok = run_workflow(
        outdir, bin_dir, pseudo_dir,
        np=np, mpi_command=mpi_command, omp_num_threads=omp,
        stage=args.stage, skip_existing=args.skip_existing, verbose=True,
    )
    print("\nrun: %s" % ("OK" if ok else "FAIL"))
    return 0 if ok else 1


def _run_parse(args: argparse.Namespace) -> int:
    cfg = _load_config(getattr(args, "config", None))
    outdir = Path(args.outdir or cfg.get("outdir", "runs"))
    if not outdir.is_dir():
        print("parse: FAIL (directory not found: %s)" % outdir, file=sys.stderr)
        return 1
    rec = parse_outdir(outdir)
    print_summary(rec)
    write_outputs(rec)
    print("  wrote raman_modes.csv + summary.json")
    return 0 if rec["status"] != "failed" else 1


def _run_spectrum(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    outdir = Path(args.outdir or cfg.get("outdir", "runs"))
    if not (outdir / "raman_modes.csv").exists():
        print("spectrum: FAIL (no raman_modes.csv in %s; run `parse` first)" % outdir,
              file=sys.stderr)
        return 1

    spectrum_type = args.type or cfg.get("spectrum_type", "activity")
    broadening = args.broadening or cfg.get("spectrum_broadening", "lorentzian")
    fwhm = args.fwhm if args.fwhm is not None else float(cfg.get("fwhm", 8.0))
    temperature = args.temperature if args.temperature is not None else float(cfg.get("temperature", 300.0))
    laser = args.laser_wavelength if args.laser_wavelength is not None else float(cfg.get("laser_wavelength", 532.0))
    threshold = args.threshold if args.threshold is not None else 1e-4
    formats = tuple(f.strip() for f in (args.formats or "png,csv").split(",") if f.strip())

    peaks = select_peaks(load_modes(outdir).get("to", []), threshold=threshold)
    if not peaks:
        print("spectrum: no Raman-active optical modes (blank spectrum)")
    else:
        for f, a in sorted(peaks, key=lambda t: t[0]):
            print("  peak %8.2f cm-1  activity %.6f" % (f, a))

    written = generate_spectra(
        outdir, spectrum_type=spectrum_type, broadening=broadening, fwhm=fwhm,
        temperature=temperature, laser_wavelength=laser, threshold=threshold,
        formats=formats,
    )
    for w in written:
        print("  wrote %s" % w)
    return 0


def _run_report(args: argparse.Namespace) -> int:
    cfg = _load_config(getattr(args, "config", None))
    outdir = Path(args.outdir or cfg.get("outdir", "runs"))
    if not outdir.is_dir():
        print("report: FAIL (directory not found: %s)" % outdir, file=sys.stderr)
        return 1
    if not (outdir / "summary.json").exists() and not (outdir / "raman_modes.csv").exists():
        print("report: FAIL (no summary.json / raman_modes.csv in %s; run `parse` first)"
              % outdir, file=sys.stderr)
        return 1
    path = write_report(outdir, render_report(outdir))
    print("  wrote %s" % path)
    return 0


def _run_all(args: argparse.Namespace) -> int:
    rc = _run_prepare(args)
    if rc != 0:
        return rc

    if not getattr(args, "execute", False):
        outdir = Path(args.outdir or _load_config(args.config).get("outdir", "runs"))
        print("\nall: inputs prepared in %s (re-run with --execute to compute)"
              % outdir)
        return 0

    for step in (_run_run, _run_parse, _run_spectrum, _run_report):
        rc = step(args)
        if rc != 0:
            print("\nall: FAIL (a step returned %d)" % rc)
            return rc
    print("\nall: OK")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "preflight":
        return _run_preflight(args)
    if args.command == "structure":
        return _run_structure(args)
    if args.command == "prepare":
        return _run_prepare(args)
    if args.command == "run":
        return _run_run(args)
    if args.command == "parse":
        return _run_parse(args)
    if args.command == "spectrum":
        return _run_spectrum(args)
    if args.command == "report":
        return _run_report(args)
    if args.command == "all":
        return _run_all(args)
    if args.command == "validate-input":
        from validate_input import main as vi_main

        return vi_main([args.file])
    if args.command == "list-templates":
        from list_templates import main as lt_main

        return lt_main([])
    return _stub(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
