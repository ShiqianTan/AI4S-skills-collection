#!/usr/bin/env python3
"""qe_converge — parameter convergence scans for the powder-Raman workflow.

Each scan varies one parameter at a time (``ecutwfc``, ``kpt_density``,
``tr2_ph``) and runs the full ``pw.x -> ph.x -> dynmat.x`` chain in a separate
subdirectory.  The ``fwhm`` group is post-processing only: it re-broadens the
reference mode table and does not need a QE run.

Outputs (in the scan root):

* ``convergence.json`` — per-variant metrics, deviations and verdicts
* ``convergence.csv``  — flat table for spreadsheets
* ``convergence.md``   — human-readable Chinese report with criteria

Criteria (deliberately matching ``examples/README.md`` tolerances):

* optical frequency: ``max |f - f_ref| <= 2 cm-1``
* Raman activity:    ``max relative change <= 5 %``
* total energy:      ``|dE| <= 1e-4 Ry/atom``
* dielectric:        ``relative change <= 1 %``
"""
from __future__ import annotations

import copy
import csv
import json
from pathlib import Path

import numpy as np

from qe_common import (
    load_defaults,
    resolve_bin_dir,
    resolve_pseudo_dir,
    validate_config,
)
from qe_inputs import generate_inputs
from qe_parse import parse_outdir
from qe_pp import build_pseudopotential_manifest, raman_xc_issues
from qe_run import run_workflow
from qe_spectrum import build_spectrum, load_modes, select_peaks
from qe_struct import (
    check_structure,
    format_check_result,
    read_structure,
    refine_structure,
    write_validated_cif,
)

CRITERIA = {
    "freq_cm1": 2.0,
    "activity_rel": 0.05,
    "etot_ry_per_atom": 1.0e-4,
    "dielectric_rel": 0.01,
}
ACTIVITY_FLOOR = 1.0e-4

DEFAULT_SWEEPS = {
    "ecutwfc": [50.0, 60.0, 70.0],
    "kpt_density": [10.0, 20.0, 30.0],
    "tr2_ph": [1.0e-10, 1.0e-12, 1.0e-14],
    "fwhm": [4.0, 6.0, 8.0, 10.0],
}
QE_GROUPS = ("ecutwfc", "kpt_density", "tr2_ph")
ALL_GROUPS = QE_GROUPS + ("fwhm",)

GROUP_LABELS = {
    "ecutwfc": "截断能 ecutwfc (Ry)",
    "kpt_density": "k 网格密度 kpt_density",
    "tr2_ph": "声子收敛阈值 tr2_ph",
    "fwhm": "展宽 FWHM (cm⁻¹)",
}


def _fmt(value) -> str:
    value = float(value)
    if value == int(value) and abs(value) < 1e6:
        return "%d" % int(value)
    return "%g" % value


def variant_tag(group: str, value) -> str:
    """Filesystem-safe tag, e.g. ``ecutwfc_60`` or ``tr2_1e-10``."""
    prefix = {"ecutwfc": "ecutwfc", "kpt_density": "kpt",
              "tr2_ph": "tr2", "fwhm": "fwhm"}[group]
    return "%s_%s" % (prefix, _fmt(value))


def variant_config(base: dict, group: str, value) -> dict:
    """Return a copy of *base* with one parameter overridden.

    ``ecutwfc`` keeps the configured ``ecutrho/ecutwfc`` ratio (default 4).
    ``kpt_density`` removes an explicit ``kpt_grid`` so that the density
    actually controls the grid.
    """
    cfg = copy.deepcopy(base)
    value = float(value)
    if group == "ecutwfc":
        base_ecut = float(base.get("ecutwfc", 60.0)) or 60.0
        eco = float(base.get("ecutrho", 4.0 * base_ecut)) or 4.0 * base_ecut
        ratio = eco / base_ecut
        if ratio < 1.0:
            ratio = 4.0
        cfg["ecutwfc"] = value
        cfg["ecutrho"] = value * ratio
    elif group == "kpt_density":
        cfg["kpt_density"] = value
        cfg.pop("kpt_grid", None)
    elif group == "tr2_ph":
        cfg["tr2_ph"] = value
    elif group == "fwhm":
        cfg["fwhm"] = value
    else:
        raise ValueError("unknown convergence group: %s" % group)
    return cfg


def sweep_values(group: str, explicit=None) -> list[float]:
    if group not in DEFAULT_SWEEPS:
        raise ValueError("unknown convergence group: %s" % group)
    if explicit:
        values = sorted({float(v) for v in explicit})
    else:
        values = list(DEFAULT_SWEEPS[group])
    if not values:
        raise ValueError("empty sweep for %s" % group)
    if group in ("ecutwfc", "kpt_density", "fwhm") and any(v <= 0 for v in values):
        raise ValueError("%s values must be positive" % group)
    if group == "tr2_ph" and any(not 0 < v < 1 for v in values):
        raise ValueError("tr2_ph values must lie in (0, 1)")
    return values


def _load_config(path: str | None) -> dict:
    cfg = load_defaults()
    if path:
        import yaml

        user = yaml.safe_load(Path(path).read_text()) or {}
        cfg.update(user)
    return cfg


def _optical_numbers(rec: dict) -> list[int]:
    return list((rec.get("stats") or {}).get("optical") or [])


def metrics_from_rec(rec: dict, n_atoms: int | None = None) -> dict:
    """Flatten the quantities used by the convergence criteria."""
    modes = rec.get("modes") or []
    stats = rec.get("stats") or {}
    scf = rec.get("scf") or {}
    dielectric = rec.get("dielectric")
    eps_avg = None
    if dielectric and len(dielectric) == 3:
        eps_avg = float(np.trace(np.asarray(dielectric, dtype=float))) / 3.0
    etot = scf.get("etot_ry")
    if etot is not None and n_atoms:
        etot = float(etot)
    return {
        "status": rec.get("status"),
        "n_modes": len(modes),
        "n_atoms": n_atoms,
        "etot_ry": etot,
        "etot_ry_per_atom": (etot / n_atoms) if etot is not None and n_atoms else None,
        "gap_ev": (rec.get("electronic") or {}).get("gap_ev"),
        "eps_avg": eps_avg,
        "max_raman": stats.get("max_raman"),
        "max_raman_freq": stats.get("max_raman_freq"),
        "sum_raman": stats.get("sum_raman"),
        "acoustic_max_abs": stats.get("acoustic_max_abs"),
        "soft_modes": stats.get("soft_modes") or [],
        "optical": _optical_numbers(rec),
        "freq_by_mode": {int(m["mode"]): float(m["freq"]) for m in modes},
        "raman_by_mode": {int(m["mode"]): float(m["raman"]) for m in modes},
        "warnings": rec.get("warnings") or [],
    }


def compare_metrics(metrics: dict, reference: dict, criteria=None) -> dict:
    """Deviations of one variant from the reference variant."""
    criteria = criteria or CRITERIA
    out = {
        "available": bool(metrics.get("n_modes")),
        "max_freq_dev": None,
        "max_freq_mode": None,
        "max_activity_rel": None,
        "max_activity_mode": None,
        "etot_dev": None,
        "eps_rel": None,
        "converged": False,
        "checks": {},
    }
    if not out["available"]:
        return out

    freq_dev = []
    for mode in reference.get("optical", []):
        f = metrics["freq_by_mode"].get(mode)
        fref = reference["freq_by_mode"].get(mode)
        if f is None or fref is None:
            continue
        freq_dev.append((abs(f - fref), mode))
    if freq_dev:
        dev, mode = max(freq_dev)
        out["max_freq_dev"] = round(dev, 4)
        out["max_freq_mode"] = mode

    act_rel = []
    for mode in reference.get("optical", []):
        a = metrics["raman_by_mode"].get(mode)
        aref = reference["raman_by_mode"].get(mode)
        if a is None or aref is None or aref <= ACTIVITY_FLOOR:
            continue
        act_rel.append((abs(a - aref) / aref, mode))
    if act_rel:
        dev, mode = max(act_rel)
        out["max_activity_rel"] = round(dev, 6)
        out["max_activity_mode"] = mode

    if metrics.get("etot_ry_per_atom") is not None and reference.get("etot_ry_per_atom") is not None:
        out["etot_dev"] = round(abs(metrics["etot_ry_per_atom"]
                                    - reference["etot_ry_per_atom"]), 8)
    if metrics.get("eps_avg") is not None and reference.get("eps_avg"):
        out["eps_rel"] = round(abs(metrics["eps_avg"] - reference["eps_avg"])
                               / abs(reference["eps_avg"]), 6)

    checks = {
        "freq": out["max_freq_dev"] is None
        or out["max_freq_dev"] <= criteria["freq_cm1"],
        "activity": out["max_activity_rel"] is None
        or out["max_activity_rel"] <= criteria["activity_rel"],
        "etot": out["etot_dev"] is None
        or out["etot_dev"] <= criteria["etot_ry_per_atom"],
        "dielectric": out["eps_rel"] is None
        or out["eps_rel"] <= criteria["dielectric_rel"],
    }
    out["checks"] = checks
    out["converged"] = all(checks.values())
    return out


def _reference_index(group: str, values: list[float]) -> int:
    return int(np.argmax(values)) if group != "tr2_ph" else int(np.argmin(values))


def evaluate_group(group: str, values: list[float], metrics: list[dict]) -> dict:
    """Compare all variants of one QE group against its reference."""
    ref_idx = _reference_index(group, values)
    reference = metrics[ref_idx]
    variants = []
    for value, m in zip(values, metrics):
        entry = {"value": float(value), "tag": variant_tag(group, value),
                 "metrics": m, "is_reference": False}
        if m.get("n_modes"):
            entry.update(compare_metrics(m, reference))
        else:
            entry["available"] = False
            entry["converged"] = False
        entry["is_reference"] = len(variants) == ref_idx
        variants.append(entry)

    usable = [v for v in variants if v.get("available")]
    recommended = None
    if usable and all(v["converged"] for v in usable):
        if group == "tr2_ph":
            recommended = max(v["value"] for v in usable)
        else:
            recommended = min(v["value"] for v in usable)
    return {
        "group": group,
        "label": GROUP_LABELS[group],
        "values": [float(v) for v in values],
        "reference_value": float(values[ref_idx]),
        "variants": variants,
        "converged": bool(usable) and all(v["converged"] for v in usable),
        "recommended": recommended,
        "criteria": dict(CRITERIA),
    }


def evaluate_fwhm(values: list[float], reference_dir: Path) -> dict:
    """Post-processing scan: peak height / integrated area vs FWHM."""
    tables = load_modes(reference_dir)
    peaks = select_peaks(tables.get("to", []))
    variants = []
    for value in values:
        entry = {"value": float(value), "tag": variant_tag("fwhm", value)}
        if not peaks:
            entry.update({"available": False, "max_intensity": None,
                          "integral": None})
            variants.append(entry)
            continue
        x, y, _ = build_spectrum(peaks, fwhm=float(value), npoints=4000)
        integrate = getattr(np, "trapezoid", None) or np.trapz
        entry.update({
            "available": True,
            "max_intensity": float(np.max(y)),
            "integral": float(integrate(y, x)),
            "max_freq": float(x[int(np.argmax(y))]),
            "note": "后处理参数：改变峰宽与重叠（孤立峰峰值不变，积分随 FWHM 增大）",
        })
        variants.append(entry)
    return {
        "group": "fwhm",
        "label": GROUP_LABELS["fwhm"],
        "values": [float(v) for v in values],
        "reference_dir": str(reference_dir),
        "variants": variants,
        "converged": None,
        "recommended": None,
        "note": "FWHM 只影响展宽后的峰高/积分，没有 DFT 收敛判据。",
    }


def _prepare_variant(refined, cfg: dict, outdir: Path, pseudo_dir: Path,
                     prefix: str, lo_dirs, bin_dir, structure_info: dict) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    write_validated_cif(refined, outdir)
    (outdir / "structure_manifest.json").write_text(
        json.dumps(structure_info, indent=2, ensure_ascii=False) + "\n")
    (outdir / "effective_config.json").write_text(
        json.dumps(cfg, indent=2, ensure_ascii=False, default=str) + "\n")
    manifest = build_pseudopotential_manifest(structure_info["elements"],
                                              pseudo_dir, outdir)
    if manifest["errors"]:
        raise RuntimeError("pseudopotential problems: " + "; ".join(manifest["errors"]))
    issues = [m for m in raman_xc_issues(manifest, cfg.get("raman", True))
              if m[0] == "FAIL"]
    if issues:
        raise RuntimeError(issues[0][2])
    generate_inputs(refined, cfg, outdir, pseudo_dir, prefix,
                    lo_dirs=lo_dirs, bin_dir=bin_dir)


def _write_outputs(root: Path, result: dict) -> list[Path]:
    written = []
    root.mkdir(parents=True, exist_ok=True)
    (root / "convergence.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, default=str) + "\n")
    written.append(root / "convergence.json")

    header = ["group", "tag", "value", "status", "n_modes", "etot_ry_per_atom",
              "gap_ev", "eps_avg", "max_raman", "max_raman_freq", "sum_raman",
              "max_freq_dev", "max_activity_rel", "etot_dev", "eps_rel",
              "converged", "note"]
    with (root / "convergence.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for group in result["groups"]:
            for v in group["variants"]:
                m = v.get("metrics") or {}
                writer.writerow([
                    group["group"], v.get("tag"), v.get("value"),
                    m.get("status", "prepared" if group["group"] in QE_GROUPS else ""),
                    m.get("n_modes", ""), m.get("etot_ry_per_atom", ""),
                    m.get("gap_ev", ""), m.get("eps_avg", ""),
                    m.get("max_raman", ""), m.get("max_raman_freq", ""),
                    m.get("sum_raman", ""),
                    v.get("max_freq_dev", ""), v.get("max_activity_rel", ""),
                    v.get("etot_dev", ""), v.get("eps_rel", ""),
                    v.get("converged", ""),
                    v.get("note") or ("参考" if v.get("is_reference") else ""),
                ])
    written.append(root / "convergence.csv")

    (root / "convergence.md").write_text(render_markdown(result))
    written.append(root / "convergence.md")
    return written


def _num(value, nd=4) -> str:
    return "—" if value is None else ("%.*f" % (nd, float(value)))


def _mark(flag) -> str:
    return "✅" if flag else "❌"


def render_markdown(result: dict) -> str:
    c = result["criteria"]
    lines = [
        "# 收敛扫描报告",
        "",
        "- 结构：`%s`" % result["structure"],
        "- 扫描根目录：`%s`" % result["outdir"],
        "- 模式：%s" % ("已执行 QE（可直接判读）"
                        if result["execute"] else "仅生成输入/解析现有日志"),
        "- 判据：峰位 |Δf| ≤ %.1f cm⁻¹；拉曼活性相对变化 ≤ %.0f%%；"
        "总能 |ΔE| ≤ %.1e Ry/atom；介电常数相对变化 ≤ %.0f%%"
        % (c["freq_cm1"], 100 * c["activity_rel"],
           c["etot_ry_per_atom"], 100 * c["dielectric_rel"]),
        "",
    ]

    for group in result["groups"]:
        lines.append("## %s" % group["label"])
        lines.append("")
        if group["group"] == "fwhm":
            lines.append("> %s" % group["note"])
            lines.append("")
            lines.append("| FWHM | 峰值位置 (cm⁻¹) | 峰值强度 | 积分面积 |")
            lines.append("|---|---|---|---|")
            for v in group["variants"]:
                lines.append("| %s | %s | %s | %s |" % (
                    _fmt(v["value"]), _num(v.get("max_freq"), 2),
                    _num(v.get("max_intensity"), 1), _num(v.get("integral"), 1)))
            lines.append("")
            continue

        lines.append("| 参数 | 状态 | 总能/atom (Ry) | 带隙 (eV) | ε_avg | "
                     "最强 Raman (cm⁻¹) | ΣRaman | Δf_max | Δ活性 | ΔE | Δε | 判定 |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for v in group["variants"]:
            m = v.get("metrics") or {}
            mark = "参考" if v.get("is_reference") else _mark(v.get("converged"))
            if not v.get("available"):
                mark = ("待运行" if m.get("status") == "prepared"
                        else "失败")
            lines.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                _fmt(v["value"]), m.get("status") or "prepared",
                _num(m.get("etot_ry_per_atom"), 6), _num(m.get("gap_ev"), 3),
                _num(m.get("eps_avg"), 4),
                ("%s (%s)" % (_num(m.get("max_raman"), 2), _num(m.get("max_raman_freq"), 2))
                 if m.get("max_raman") is not None else "—"),
                _num(m.get("sum_raman"), 3),
                _num(v.get("max_freq_dev"), 2),
                ("%.2f%%" % (100 * v["max_activity_rel"])
                 if v.get("max_activity_rel") is not None else "—"),
                _num(v.get("etot_dev"), 6),
                ("%.2f%%" % (100 * v["eps_rel"])
                 if v.get("eps_rel") is not None else "—"),
                mark,
            ))
        lines.append("")
        if group["recommended"] is not None:
            lines.append("**推荐最小可用值：%s**（所有已运行变体均在容差内）"
                         % _fmt(group["recommended"]))
        elif any(v.get("available") for v in group["variants"]):
            lines.append("**尚未收敛**：至少一个变体超出容差；请加密参数后重跑。")
        else:
            lines.append("尚未执行 QE；使用 `--execute` 运行后重新生成报告。")
        failures = [v for v in group["variants"] if v.get("metrics")
                    and (v["metrics"].get("status") == "failed")]
        if failures:
            lines.append("")
            lines.append("失败变体：" + ", ".join(
                "%s (%s)" % (v["tag"], (v["metrics"].get("status_reason") or "failed"))
                for v in failures))
        lines.append("")

    lines.append("## 备注")
    lines.append("")
    lines.append("- 参考变体取最严格/最密的一侧（ecutwfc 与 kpt_density 取最大值，"
                 "tr2_ph 取最小值）。")
    lines.append("- 峰位按模式编号与参考变体配对；结构对称性不变时模式顺序稳定。")
    lines.append("- 本扫描只验证数值收敛，不替代实验对比；LDA/HGH 与 PBE/ONCV 的"
                 "系统误差不在判据内。")
    lines.append("")
    return "\n".join(lines)


def _find_reference_dir(root: Path, groups: list[dict], fallback: Path | None) -> Path | None:
    for name in ("ecutwfc", "kpt_density", "tr2_ph"):
        group = next((g for g in groups if g["group"] == name), None)
        if not group:
            continue
        for v in group["variants"]:
            if v.get("available"):
                return root / name / v["tag"]
    if fallback and (fallback / "raman_modes.csv").exists():
        return fallback
    return None


def run(args) -> int:
    from qe_raman import _apply_cli_overrides  # local import avoids a cycle

    cfg = _apply_cli_overrides(_load_config(args.config), args)
    for level, key, text in validate_config(cfg):
        print("[%s] %-14s %s" % (level, key, text))
        if level == "FAIL":
            print("\nconverge: FAIL (invalid config; see above)")
            return 1

    groups = [g.strip() for g in (args.groups or "").split(",") if g.strip()]
    unknown = [g for g in groups if g not in ALL_GROUPS]
    if unknown:
        print("converge: FAIL (unknown groups: %s)" % ", ".join(unknown))
        return 1
    if not groups:
        groups = list(ALL_GROUPS)

    if args.outdir:
        root = Path(args.outdir)
    else:
        root = Path(cfg.get("outdir", "runs")) / "convergence"
    root.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or cfg.get("prefix", "raman")
    pseudo_dir = resolve_pseudo_dir(args.pseudo_dir or cfg.get("pseudo_dir"))
    bin_dir = resolve_bin_dir(getattr(args, "bin_dir", None) or cfg.get("bin_dir"))
    lo_dirs = [tuple(q) for q in (cfg.get("lo_to_directions") or [])]

    atoms = read_structure(args.structure)
    res = check_structure(args.structure)
    print("--- structure checks ---")
    print(format_check_result(res))
    if res.has_fail:
        print("\nconverge: FAIL (structure issues)")
        return 1
    if args.no_refine:
        refined, info = atoms, {
            "standardized": False,
            "cellpar": [round(float(x), 8) for x in atoms.cell.cellpar()],
            "elements": list(atoms.get_chemical_symbols()),
        }
    else:
        refined, info = refine_structure(atoms, symprec=args.symprec)
    info["source"] = str(args.structure)
    print("reference cellpar: %s"
          % ", ".join("%.6f" % x for x in info["cellpar"]))
    n_atoms = len(refined)

    explicit = {
        "ecutwfc": getattr(args, "ecutwfc", None),
        "kpt_density": getattr(args, "kpt_values", None),
        "tr2_ph": getattr(args, "tr2_ph", None),
        "fwhm": getattr(args, "fwhm", None),
    }

    result = {
        "structure": str(args.structure),
        "outdir": str(root),
        "execute": bool(args.execute),
        "groups": [],
        "criteria": dict(CRITERIA),
    }

    np_ = args.np if args.np is not None else int(cfg.get("np", 4))
    mpi = args.mpi_command or cfg.get("mpi_command", "mpirun")
    omp = (args.omp_num_threads if args.omp_num_threads is not None
           else int(cfg.get("omp_num_threads", 1)))

    for group in groups:
        if group not in QE_GROUPS:
            continue
        values = sweep_values(group, explicit.get(group))
        print("\n=== %s: %s ===" % (group, ", ".join(_fmt(v) for v in values)))
        metrics = []
        group_dir = root / group
        for value in values:
            tag = variant_tag(group, value)
            vdir = group_dir / tag
            vcfg = variant_config(cfg, group, value)
            prepared = (vdir / "scf.in").exists()
            if not prepared or args.refresh:
                try:
                    _prepare_variant(refined, vcfg, vdir, pseudo_dir, prefix,
                                     lo_dirs, bin_dir, info)
                    print("  [%s] inputs -> %s" % (tag, vdir))
                except (RuntimeError, FileNotFoundError, ValueError) as exc:
                    print("  [%s] PREPARE FAIL: %s" % (tag, exc))
                    metrics.append({"status": "failed", "n_modes": 0,
                                    "error": str(exc)})
                    continue
            else:
                print("  [%s] reusing inputs" % tag)

            if args.execute and (args.refresh or not (vdir / "scf.out").exists()):
                _, ok = run_workflow(
                    vdir, bin_dir, pseudo_dir, np=np_, mpi_command=mpi,
                    omp_num_threads=omp, stage="all", skip_existing=True,
                    verbose=args.verbose,
                )
                if not ok:
                    print("  [%s] RUN FAIL" % tag)
            if (vdir / "scf.out").exists():
                rec = parse_outdir(vdir)
                m = metrics_from_rec(rec, n_atoms=n_atoms)
                m["status"] = rec.get("status")
                m["status_reason"] = rec.get("status_reason")
                metrics.append(m)
                print("  [%s] parse: status=%s, modes=%d, eps=%s"
                      % (tag, m["status"], m["n_modes"],
                         "—" if m.get("eps_avg") is None else "%.4f" % m["eps_avg"]))
            else:
                metrics.append({"status": "prepared", "n_modes": 0})
        result["groups"].append(evaluate_group(group, values, metrics))

    fwhm_group = None
    if "fwhm" in groups:
        values = sweep_values("fwhm", explicit["fwhm"])
        ref_dir = _find_reference_dir(root, result["groups"], root.parent)
        if ref_dir is None:
            print("\n=== fwhm: no reference modes.csv; skipped ===")
        else:
            print("\n=== fwhm: %s (reference %s) ==="
                  % (", ".join(_fmt(v) for v in values), ref_dir))
            fwhm_group = evaluate_fwhm(values, ref_dir)
            result["groups"].append(fwhm_group)

    for path in _write_outputs(root, result):
        print("  wrote %s" % path)

    print("\n=== summary ===")
    okay = True
    any_available = False
    for group in result["groups"]:
        if group["group"] == "fwhm":
            print("  fwhm: post-processing only (%d values)" % len(group["variants"]))
            continue
        available = [v for v in group["variants"] if v.get("available")]
        any_available = any_available or bool(available)
        if not available:
            print("  %s: inputs prepared (not run)" % group["group"])
        elif group["converged"]:
            print("  %s: converged (recommended %s)"
                  % (group["group"], _fmt(group["recommended"])))
        else:
            print("  %s: NOT converged" % group["group"])
            okay = False
    if not args.execute and any(g["group"] in QE_GROUPS for g in result["groups"]):
        print("  inputs generated only; re-run with --execute to evaluate convergence")
        okay = okay and any(v.get("status") != "failed"
                            for g in result["groups"] for v in g["variants"])
    return 0 if okay else 1


def add_arguments(parser) -> None:
    parser.add_argument("structure", help="CIF / POSCAR / PW input file")
    parser.add_argument("--config", default=None, help="YAML config file")
    parser.add_argument("--outdir", default=None,
                        help="scan root (default: <outdir>/convergence)")
    parser.add_argument("--pseudo-dir", default=None, help="pseudopotential directory")
    parser.add_argument("--prefix", default=None, help="QE prefix")
    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="spglib symmetry tolerance")
    parser.add_argument("--no-refine", action="store_true",
                        help="skip spglib.refine_cell")
    parser.add_argument("--no-lraman", action="store_true",
                        help="phonon/IR only (needed for GGA pseudopotentials)")
    parser.add_argument("--groups", default="ecutwfc,kpt_density,tr2_ph,fwhm",
                        help="comma-separated groups to scan")
    parser.add_argument("--ecutwfc", nargs="*", type=float, default=None,
                        help="explicit ecutwfc values (Ry)")
    parser.add_argument("--kpt-values", "--kpt-density-values", nargs="*",
                        type=float, default=None, dest="kpt_values",
                        help="explicit kpt_density values to scan")
    parser.add_argument("--kpt-density", type=float, default=None,
                        help="base kpt_density for non-kpt scans")
    parser.add_argument("--kpt-grid", nargs=3, type=int, default=None,
                        metavar=("NX", "NY", "NZ"),
                        help="explicit base k-grid for non-kpt scans")
    parser.add_argument("--tr2-ph", nargs="*", type=float, default=None,
                        help="explicit tr2_ph values")
    parser.add_argument("--fwhm", nargs="*", type=float, default=None,
                        help="explicit FWHM values (cm-1)")
    parser.add_argument("--execute", action="store_true",
                        help="run pw.x/ph.x/dynmat.x for each variant")
    parser.add_argument("--refresh", action="store_true",
                        help="regenerate inputs and rerun even if logs exist")
    parser.add_argument("--np", type=int, default=None, help="MPI processes")
    parser.add_argument("--bin-dir", default=None, help="QE binary directory")
    parser.add_argument("--mpi-command", default=None, help="MPI launcher")
    parser.add_argument("--omp-num-threads", type=int, default=None)
    parser.add_argument("--verbose", action="store_true", default=True,
                        help="stream QE stage progress")
