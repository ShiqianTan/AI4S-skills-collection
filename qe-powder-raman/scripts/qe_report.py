#!/usr/bin/env python3
"""qe_report — assemble a human-readable ``report.md`` from a parsed run.

Reads ``summary.json`` (SCF/dielectric/ASR/stats/LO) and ``raman_modes.csv``
(full mode tables) in the run directory and writes ``report.md``.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


def _load_summary(outdir: Path) -> dict:
    path = outdir / "summary.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:  # noqa: BLE001
        return {}


def _load_mode_tables(outdir: Path) -> dict[str, list[dict]]:
    path = outdir / "raman_modes.csv"
    tables: dict[str, list[dict]] = {}
    if not path.exists():
        return tables
    with open(path) as fh:
        for r in csv.DictReader(fh):
            t = r["table"]
            m = {
                "mode": int(r["mode"]),
                "freq": float(r["freq_cm1"]),
                "thz": float(r["freq_thz"]),
                "ir": float(r["ir"]),
                "raman": float(r["raman"]),
                "depol": float(r["depol"]) if r["depol"] != "" else None,
            }
            tables.setdefault(t, []).append(m)
    return tables


def _fmt(v, nd: int = 4) -> str:
    if v is None:
        return "—"
    return ("%." + str(nd) + "f") % v


def _mode_table(title: str, modes: list[dict]) -> list[str]:
    lines = ["### %s" % title, "",
             "| 模 | cm⁻¹ | THz | IR | Raman | 退偏 |",
             "|----|------|-----|----|-------|------|"]
    for m in modes:
        dep = "—" if m["depol"] is None else _fmt(m["depol"])
        lines.append("| %d | %s | %s | %s | %s | %s |"
                     % (m["mode"], _fmt(m["freq"], 2), _fmt(m["thz"], 4),
                        _fmt(m["ir"]), _fmt(m["raman"]), dep))
    return lines + [""]


def render_report(outdir: Path) -> str:
    s = _load_summary(outdir)
    tables = _load_mode_tables(outdir)

    lines: list[str] = []
    lines.append("# QE 粉末拉曼计算报告")
    lines.append("")
    lines.append("- 运行目录：`%s`" % outdir)
    lines.append("- 生成时间：%s" % datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"))
    status = s.get("status", "unknown")
    lines.append("- 状态：**%s**" % status)
    lines.append("")

    scf = s.get("scf") or {}
    lines.append("## 1. SCF 参数")
    lines.append("")
    lines.append("| 项目 | 值 |")
    lines.append("|------|----|")
    rows = [
        ("电子数 nelec", scf.get("nelec")),
        ("能带数 nbnd", scf.get("nbnd")),
        ("k 点数 nk", scf.get("nk")),
        ("截断能 ecutwfc (Ry)", scf.get("ecutwfc")),
        ("对称操作数", scf.get("sym_ops")),
        ("SCF 迭代数", scf.get("scf_iter")),
        ("总能量 (Ry)", scf.get("etot_ry")),
        ("体积 (Å³)", scf.get("volume_a3")),
    ]
    for k, v in rows:
        lines.append("| %s | %s |" % (k, _fmt(v, 6) if isinstance(v, float) else (v if v is not None else "—")))
    lines.append("")

    lines.append("## 2. 介电性质与声子")
    lines.append("")
    ph = [
        ("ph 收敛", "是" if s.get("ph_converged") else "否"),
        ("ph JOB DONE", "是" if s.get("ph_job_done") else "否"),
        ("声学求和规则 ASR", _fmt(s.get("asr_metric"))),
    ]
    for k, v in ph:
        lines.append("- %s：%s" % (k, v))
    diel = s.get("dielectric")
    if diel:
        diag = [_fmt(diel[i][i]) for i in range(3)]
        lines.append("- 介电常数对角元：%s" % ", ".join(diag))
    dyn = s.get("dynmat") or {}
    if dyn.get("cm_factor") is not None:
        lines.append("- Clausius-Mossotti 因子：%s" % _fmt(dyn["cm_factor"]))
    if dyn.get("raman_cm_factor") is not None:
        lines.append("- Raman CM 因子：%s" % _fmt(dyn["raman_cm_factor"]))
    if dyn.get("polarizability"):
        p = dyn["polarizability"]
        lines.append("- 极化率 (Å³) 对角元：%s" % ", ".join(_fmt(p[i][i]) for i in range(3)))
    if not dyn.get("has_raman", True):
        lines.append("- 注意：dynmat 未输出 Raman 列（该体系无一级拉曼活性）。")
    lines.append("")

    st = s.get("stats") or {}
    lines.append("## 3. 模统计")
    lines.append("")
    lines.append("| 项目 | 值 |")
    lines.append("|------|----|")
    stats_rows = [
        ("总模数", st.get("total")),
        ("光学模数", len(st.get("optical", []))),
        ("拉曼活性模数", len(st.get("raman_active", []))),
        ("红外活性模数", len(st.get("ir_active", []))),
        ("Raman 活性和", _fmt(st.get("sum_raman"))),
        ("最强 Raman 活性", _fmt(st.get("max_raman"))),
        ("最强 Raman 频率 (cm⁻¹)", _fmt(st.get("max_raman_freq"), 2)),
    ]
    for k, v in stats_rows:
        lines.append("| %s | %s |" % (k, v if v is not None else "—"))
    if st.get("soft_modes"):
        lines.append("")
        lines.append("**警告：存在虚频（软模）**，结构未充分弛豫。")
    lines.append("")

    to_modes = tables.get("to", [])
    lines.append("## 4. 模表（横向光学 TO，Γ 点）")
    lines.append("")
    if to_modes:
        lines.extend(_mode_table("TO 模", to_modes))
    else:
        lines.append("（无）")
        lines.append("")

    lo_keys = sorted(k for k in tables if k != "to")
    if lo_keys:
        lines.append("## 5. LO-TO 分裂（纵向光学 LO）")
        lines.append("")
        for k in lo_keys:
            lines.extend(_mode_table("LO 表 `%s`" % k, tables[k]))

    files = [
        "raman_modes.csv", "summary.json",
        "spectrum_activity.png", "spectrum_stokes.png",
    ]
    lines.append("## 6. 输出文件")
    lines.append("")
    for f in files:
        mark = "（存在）" if (outdir / f).exists() else "（未生成）"
        lines.append("- `%s` %s" % (f, mark))
    lines.append("")

    lines.append("> 本报告由 `qe_report.py` 自动生成。")
    return "\n".join(lines) + "\n"


def write_report(outdir: Path, text: str) -> Path:
    path = Path(outdir) / "report.md"
    path.write_text(text)
    return path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="qe_report.py", description="assemble report.md")
    p.add_argument("--outdir", required=True, help="run directory with summary.json")
    args = p.parse_args(argv)

    outdir = Path(args.outdir)
    if not (outdir / "summary.json").exists() and not (outdir / "raman_modes.csv").exists():
        print("report: FAIL (no summary.json / raman_modes.csv in %s; run `parse` first)"
              % outdir, file=__import__("sys").stderr)
        return 1

    path = write_report(outdir, render_report(outdir))
    print("wrote %s" % path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
