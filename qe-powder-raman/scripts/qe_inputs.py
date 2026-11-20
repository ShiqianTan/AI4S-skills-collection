#!/usr/bin/env python3
"""
qe_inputs — render Quantum ESPRESSO input files (pw.x / ph.x / dynmat.x) for the
solid-state powder-Raman workflow.

Generated inputs implement the validated recipe:

    pw.x scf (LDA + NC, occupations='fixed')
      -> ph.x (epsil=.true. trans=.true. lraman=.true. asr=.true., q=Gamma)
      -> dynmat.x (asr='crystal') -> IR/Raman activities + depolarization factors

Only norm-conserving LDA pseudopotentials are supported (lraman forbids US/PAW).
Optional non-analytic LO-TO runs are written as separate dynmat_lo*.in inputs
using the `q` (3-vector) namelist variable (NOT q(3)).
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
from ase.data import atomic_masses, chemical_symbols

from qe_common import DEFAULT_BIN_DIR, TEMPLATE_DIR, load_defaults, resolve_pseudo_dir


# ---------------------------------------------------------------------------
# rendering helpers
# ---------------------------------------------------------------------------


def render(text: str, variables: dict) -> str:
    """Substitute {key} placeholders (leaves bash ${key} and other braces untouched)."""
    def repl(match: re.Match) -> str:
        key = match.group(1)
        val = variables.get(key)
        return "" if val is None else str(val)

    # Negative lookbehind for '$' keeps bash "${VAR}" references untouched.
    return re.sub(r"(?<!\$)\{([A-Za-z_][A-Za-z0-9_]*)\}", repl, text)


def fmt_num(val) -> str:
    """Format a Python number in a QE-namelist-friendly way."""
    if isinstance(val, bool):
        return ".true." if val else ".false."
    if isinstance(val, (int, np.integer)):
        return str(int(val))
    f = float(val)
    if f.is_integer() and abs(f) < 1e12:
        return str(int(f))
    if abs(f) >= 1e-3:
        return "%g" % f
    return "%.1e" % f


def species_order(atoms) -> list[str]:
    """Unique element symbols in order of first appearance."""
    seen: list[str] = []
    for s in atoms.get_chemical_symbols():
        if s not in seen:
            seen.append(s)
    return seen


def kpoint_grid(cell, density: float, shift=(0, 0, 0)) -> tuple[int, int, int]:
    """Monkhorst-Pack grid from a real-space length target.

    ``k_i = max(1, round(density / |a_i|))`` with ``|a_i|`` the i-th lattice
    vector length in Angstrom.  ``density`` is therefore a real-space sampling
    target (Angstrom of cell per k-point along each axis), NOT a
    reciprocal-Angstrom spacing.  For an explicit grid use ``kpt_grid``.
    """
    lengths = np.linalg.norm(cell, axis=1)
    grid = [max(1, int(round(density / float(L)))) for L in lengths]
    return (grid[0], grid[1], grid[2])


def resolve_kpoint_grid(cfg, cell) -> tuple[int, int, int]:
    """Return the k-grid: an explicit ``kpt_grid`` wins over ``kpt_density``."""
    kg = cfg.get("kpt_grid")
    if kg:
        vals = [max(1, int(x)) for x in kg]
        return (vals[0], vals[1], vals[2])
    return kpoint_grid(cell, float(cfg.get("kpt_density", 30.0)))


# ---------------------------------------------------------------------------
# namelist block builders
# ---------------------------------------------------------------------------


def atomic_species_block(species, pp_map) -> str:
    lines = ["ATOMIC_SPECIES"]
    for s in species:
        z = chemical_symbols.index(s)
        fname = Path(pp_map[s]).name
        lines.append("  %-3s %10.4f  %s" % (s, atomic_masses[z], fname))
    return "\n".join(lines)


def cell_block(atoms) -> str:
    lines = ["CELL_PARAMETERS angstrom"]
    for v in atoms.get_cell():
        lines.append("  %16.10f %16.10f %16.10f" % tuple(float(x) for x in v))
    return "\n".join(lines)


def positions_block(atoms) -> str:
    lines = ["ATOMIC_POSITIONS crystal"]
    for sym, pos in zip(atoms.get_chemical_symbols(), atoms.get_scaled_positions()):
        lines.append("  %-3s %16.10f %16.10f %16.10f" % (sym, pos[0], pos[1], pos[2]))
    return "\n".join(lines)


def kgrid_block(grid, shift) -> str:
    return "K_POINTS automatic\n  %d %d %d %d %d %d" % (
        grid[0], grid[1], grid[2], shift[0], shift[1], shift[2],
    )


def amass_block(species) -> str:
    lines = []
    for i, s in enumerate(species, start=1):
        z = chemical_symbols.index(s)
        lines.append("  amass(%d) = %.4f" % (i, atomic_masses[z]))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# input renderers
# ---------------------------------------------------------------------------


def _scf_vars(cfg, atoms, species, pp_map, prefix, pseudo_dir):
    ecutwfc = float(cfg["ecutwfc"])
    ecutrho = max(float(cfg.get("ecutrho", 4 * ecutwfc)), 4 * ecutwfc)
    shift = tuple(cfg.get("kpt_shift", [0, 0, 0]))
    grid = resolve_kpoint_grid(cfg, atoms.get_cell())
    return {
        "prefix": prefix,
        "outdir": cfg.get("qe_outdir", "./tmp"),
        "pseudo_dir": str(pseudo_dir),
        "verbosity": cfg.get("verbosity", "high"),
        "nat": len(atoms),
        "ntyp": len(species),
        "ecutwfc": fmt_num(ecutwfc),
        "ecutrho": fmt_num(ecutrho),
        "occupations": cfg.get("occupations", "fixed"),
        "conv_thr": fmt_num(cfg.get("conv_thr", 1.0e-10)),
        "mixing_beta": fmt_num(cfg.get("mixing_beta", 0.5)),
        "atomic_species": atomic_species_block(species, pp_map),
        "cell": cell_block(atoms),
        "positions": positions_block(atoms),
        "kgrid": kgrid_block(grid, shift),
    }


def generate_scf(cfg, atoms, species, pp_map, prefix, pseudo_dir) -> str:
    from qe_pp import parse_upf

    tpl = (TEMPLATE_DIR / "scf.in.tmpl").read_text()
    v = _scf_vars(cfg, atoms, species, pp_map, prefix, pseudo_dir)
    v["calculation"] = cfg.get("calculation", "scf")
    valence = {s: parse_upf(pp_map[s]).z_valence for s in species}
    v["bands"] = ""
    v["empty_band_accuracy"] = ""
    if all(z is not None and np.isfinite(z) for z in valence.values()):
        occupied = sum(valence[s] for s in atoms.get_chemical_symbols()) / 2
        if occupied > 0 and np.isclose(occupied, round(occupied), rtol=0, atol=1e-8):
            v["bands"] = f"  nbnd         = {round(occupied) + 1}\n"
            v["empty_band_accuracy"] = "  diago_full_acc = .true.\n"
    return render(tpl, v)


def generate_relax(cfg, atoms, species, pp_map, prefix, pseudo_dir) -> str:
    tpl = (TEMPLATE_DIR / "relax.in.tmpl").read_text()
    return render(tpl, _scf_vars(cfg, atoms, species, pp_map, prefix, pseudo_dir))


def generate_ph(cfg, species, prefix) -> str:
    tpl = (TEMPLATE_DIR / "ph.in.tmpl").read_text()
    return render(tpl, {
        "prefix": prefix,
        "outdir": cfg.get("qe_outdir", "./tmp"),
        "amass": amass_block(species),
        "tr2_ph": fmt_num(cfg.get("tr2_ph", 1.0e-14)),
        "fildyn": prefix + ".dyn",
    })


def generate_dynmat(prefix, filout="dynmat_table.dat") -> str:
    tpl = (TEMPLATE_DIR / "dynmat.in.tmpl").read_text()
    return render(tpl, {"fildyn": prefix + ".dyn", "filout": filout})


def generate_dynmat_lo(prefix, q, filout) -> str:
    tpl = (TEMPLATE_DIR / "dynmat_lo.in.tmpl").read_text()
    return render(tpl, {
        "fildyn": prefix + ".dyn",
        "qx": fmt_num(q[0]), "qy": fmt_num(q[1]), "qz": fmt_num(q[2]),
        "filout": filout,
    })


def generate_runsh(cfg, bin_dir, pseudo_dir) -> str:
    tpl = (TEMPLATE_DIR / "run.sh.tmpl").read_text()
    return render(tpl, {
        "np": int(cfg.get("np", 4)),
        "bin_dir": str(bin_dir),
        "mpi_command": cfg.get("mpi_command", "mpirun"),
        "omp": int(cfg.get("omp_num_threads", 1)),
        "pseudo_dir": str(pseudo_dir),
    })


# ---------------------------------------------------------------------------
# top-level generator
# ---------------------------------------------------------------------------


def generate_inputs(atoms, cfg, outdir, pseudo_dir, prefix, lo_dirs=None,
                    with_relax=False, bin_dir=None) -> dict:
    """Write all inputs + run.sh into outdir. Returns a summary dict."""
    from qe_pp import element_pp_map

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    species = species_order(atoms)
    pp_map = element_pp_map(pseudo_dir, species)
    pp_single: dict[str, str] = {}
    for s in species:
        cands = pp_map.get(s) or []
        if not cands:
            raise FileNotFoundError(
                "no UPF for element %s in %s" % (s, pseudo_dir)
            )
        if len(cands) > 1:
            raise ValueError(
                "ambiguous UPF for %s: %d candidates (%s)"
                % (s, len(cands), ", ".join(cands))
            )
        pp_single[s] = cands[0]

    lo_dirs = lo_dirs or []
    clean_lo = []
    for q in lo_dirs:
        q = tuple(float(x) for x in q)
        if len(q) != 3 or not np.isfinite(q).all():
            raise ValueError("LO direction must be 3 finite numbers: %r" % (q,))
        if float(np.linalg.norm(q)) == 0.0:
            raise ValueError("LO direction must be non-zero: %r" % (q,))
        clean_lo.append(q)
    lo_dirs = clean_lo
    shift = tuple(cfg.get("kpt_shift", [0, 0, 0]))
    grid = resolve_kpoint_grid(cfg, atoms.get_cell())

    # clean stale LO artifacts from a previous prepare (P0-05)
    for stale in outdir.glob("dynmat_lo*.in"):
        stale.unlink()
    for stale in outdir.glob("dynmat_lo*.dat"):
        stale.unlink()
    for stale in outdir.glob("dynmat_lo*.out"):
        stale.unlink()
    dt = outdir / "dynmat_table.dat"
    if dt.exists():
        dt.unlink()

    written: list[str] = []

    def _write(name, content):
        (outdir / name).write_text(content)
        written.append(name)

    if with_relax:
        _write("relax.in", generate_relax(cfg, atoms, species, pp_single, prefix, pseudo_dir))
    _write("scf.in", generate_scf(cfg, atoms, species, pp_single, prefix, pseudo_dir))
    _write("ph.in", generate_ph(cfg, species, prefix))
    _write("dynmat.in", generate_dynmat(prefix, "dynmat_table.dat"))

    for i, q in enumerate(lo_dirs):
        base = "dynmat_lo" if len(lo_dirs) == 1 else "dynmat_lo_%d" % (i + 1)
        _write(base + ".in", generate_dynmat_lo(prefix, q, base + ".dat"))

    _write("run.sh", generate_runsh(cfg, bin_dir, pseudo_dir))
    (outdir / "run.sh").chmod(0o755)

    manifest = {
        "generator": "qe-powder-raman",
        "prefix": prefix,
        "species": species,
        "kpoints": list(grid),
        "qe_outdir": cfg.get("qe_outdir", "./tmp"),
        "lo_dirs": [list(q) for q in lo_dirs],
        "files": written,
    }
    (outdir / "inputs_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    return {
        "prefix": prefix,
        "species": species,
        "kpoints": grid,
        "lo_dirs": list(lo_dirs),
        "files": written,
    }


# ---------------------------------------------------------------------------
# standalone CLI (smoke-testing; the authoritative path is `qe_raman.py prepare`)
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    import argparse

    from qe_struct import read_structure, refine_structure

    p = argparse.ArgumentParser(prog="qe_inputs.py", description="render QE powder-Raman inputs")
    p.add_argument("structure")
    p.add_argument("--outdir", default="runs")
    p.add_argument("--config", default=None)
    p.add_argument("--pseudo-dir", default=None)
    p.add_argument("--prefix", default="raman")
    p.add_argument("--relax", action="store_true")
    p.add_argument("--lo", nargs="*", type=float, default=[],
                   help="LO directions as 'x y z [x y z ...]'")
    args = p.parse_args(argv)

    cfg = load_defaults()
    if args.config:
        import yaml

        user = yaml.safe_load(Path(args.config).read_text()) or {}
        cfg.update(user)

    pseudo_dir = resolve_pseudo_dir(args.pseudo_dir)
    atoms = read_structure(args.structure)
    atoms, _info = refine_structure(atoms)

    lo_dirs = [tuple(args.lo[i:i + 3]) for i in range(0, len(args.lo), 3)]
    summary = generate_inputs(atoms, cfg, args.outdir, pseudo_dir, args.prefix,
                              lo_dirs=lo_dirs, with_relax=args.relax,
                              bin_dir=DEFAULT_BIN_DIR)
    print("wrote %s:" % args.outdir)
    for f in summary["files"]:
        print("  " + f)
    print("kpoints %s, lo %s" % (summary["kpoints"], summary["lo_dirs"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
