#!/usr/bin/env python3
"""qe_report — assemble a human-readable ``report.md`` from a parsed run.

Reads ``summary.json`` (SCF/dielectric/ASR/stats/LO) and ``raman_modes.csv``
(full mode tables) in the run directory and writes ``report.md``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
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


def _load_json(outdir: Path, name: str) -> dict:
    path = outdir / name
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:  # noqa: BLE001
        return {}


def _load_convergence(outdir: Path) -> dict:
    path = outdir / "convergence" / "convergence.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:  # noqa: BLE001
        return {}


def _load_assignment(outdir: Path) -> list[dict]:
    path = outdir / "mode_assignment.csv"
    if not path.exists():
        return []
    try:
        with open(path) as fh:
            return list(csv.DictReader(fh))
    except OSError:
        return []


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
    ppman = _load_json(outdir, "pseudopotential_manifest.json")
    struct_manifest = _load_json(outdir, "structure_manifest.json")

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

    el = s.get("electronic") or {}
    lines.append("### 电子结构适用性")
    lines.append("")
    lines.append("| 项目 | 值 |")
    lines.append("|------|----|")
    lines.append("| 最高占据能级 (eV) | %s |" % _fmt(el.get("homo_ev")))
    lines.append("| 最低未占据能级 (eV) | %s |" % _fmt(el.get("lumo_ev")))
    lines.append("| 带隙 (eV) | %s |" % _fmt(el.get("gap_ev")))
    lines.append("| 判定 | %s |" % (el.get("classification") or "unknown"))
    if el.get("reason"):
        lines.append("| 依据 | %s |" % el["reason"])
    lines.append("")

    relax = s.get("relax")
    relax_stage = ((s.get("run_state") or {}).get("stages") or {}).get("relax") or {}
    lines.append("### 结构弛豫（vc-relax）")
    lines.append("")
    if relax:
        lines.append("| 项目 | 值 |")
        lines.append("|------|----|")
        if relax_stage.get("status"):
            lines.append("| 弛豫阶段状态 | %s |" % relax_stage["status"])
        lines.append("| 残余力 final_force | %s |" % _fmt(relax.get("final_force")))
        lines.append("| 残余应力 final_stress | %s |" % _fmt(relax.get("final_stress")))
        lines.append("| 原子数 | %s |" % (relax.get("n_atoms") if relax.get("n_atoms") is not None else "—"))
        if relax.get("spacegroup") is not None:
            lines.append("| 弛豫后空间群 | %s |" % relax["spacegroup"])
    else:
        lines.append("- 未执行弛豫（未启用 `--relax` 或缺少 `relax_manifest.json`）。")
    lines.append("")

    layers = struct_manifest.get("layers") or {}
    if layers.get("is_slab") or layers.get("is_layered"):
        lines.append("### 层状/低维提示")
        lines.append("")
        lines.append("| 项目 | 值 |")
        lines.append("|------|----|")
        lines.append("| 类型 | %s |" % ("slab/准二维（真空层）" if layers.get("is_slab") else "层状/范德华间隙"))
        lines.append("| 方向 | %s |" % layers.get("axis", "—"))
        lines.append("| 最大间隙 (Å) | %s |" % _fmt(layers.get("vacuum_thickness"), 2))
        lines.append("| 层数 | %s |" % layers.get("layer_count", "—"))
        if layers.get("gap_elements"):
            lines.append("| 间隙两侧元素 | %s |"
                         % "-".join(str(e) for e in layers["gap_elements"]))
        lines.append("")
        if (s.get("stats") or {}).get("soft_modes"):
            lines.append("- 检测到软模：层间/面外模式可能为虚频，建议先 vc-relax 并考虑色散修正。")
        else:
            lines.append("- 建议：层间/面外模式仍可能偏软，生产计算前建议 vc-relax；"
                         "Gamma 点流程不含层间耦合修正。")
        lines.append("")

    lines.append("## 2. 介电性质与声子")
    lines.append("")
    ph = [
        ("ph 收敛", "是" if s.get("ph_converged") else "否"),
        ("ph JOB DONE", "是" if s.get("ph_job_done") else "否"),
        ("Raman 计算 (lraman)", "开启" if s.get("lraman", True) else "关闭（仅声子/IR）"),
        ("Raman 张量生成", "是" if s.get("lraman_tensor") else "否"),
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
        lines.append("  - 备注：表中 Raman 列为 QE 原始值，未乘该因子；若启用局域场"
                     "修正，需将 Raman 乘以该因子的平方，且仅应用一次。")
    if dyn.get("polarizability"):
        p = dyn["polarizability"]
        lines.append("- 极化率 (Å³) 对角元：%s" % ", ".join(_fmt(p[i][i]) for i in range(3)))
    fams = ppman.get("xc_families") or []
    if fams:
        lines.append("- 赝势 XC 家族：%s" % ", ".join(fams))
    if ppman.get("pseudo_source"):
        lines.append("- 赝势来源：%s" % "; ".join(ppman["pseudo_source"]))
    if fams and "GGA" in fams:
        lines.append("- 注意：QE `lraman` 不支持 GGA/PBE（`third order derivatives "
                     "not implemented with GGA`），PBE NC 仅用于声子/IR（`--no-lraman`）。")
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
        lines.append("**警告：存在显著虚频（软模）**，结构可能未充分弛豫。")
    if st.get("noise_modes"):
        lines.append("")
        lines.append("- 小型数值虚频（视为噪声，不入谱）：%s"
                     % ", ".join(str(x) for x in st["noise_modes"]))
    if st.get("acoustic_clean") is False:
        lines.append("")
        lines.append("- 注意：声学模偏离零（max |freq| = %s cm⁻¹），建议检查收敛。"
                     % _fmt(st.get("acoustic_max_abs"), 2))
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
        lines.append("默认 TO 谱来自 Γ 点、未加非解析修正；下列 LO 表按指定方向单独计算。")
        lines.append("")
        lo_dirs = _load_json(outdir, "inputs_manifest.json").get("lo_dirs") or []
        for k in lo_keys:
            label = k
            idx = None
            if lo_dirs:
                if k == "lo":
                    idx = 0
                elif k.startswith("lo_"):
                    try:
                        idx = int(k.split("_")[-1]) - 1
                    except ValueError:
                        idx = None
            if idx is not None and 0 <= idx < len(lo_dirs):
                label = "%s  q=(%s)" % (k, ", ".join(_fmt(x, 3) for x in lo_dirs[idx]))
            lines.extend(_mode_table("LO 表 `%s`" % label, tables[k]))

    files = [
        "raman_modes.csv", "summary.json", "relax_manifest.json",
    ]
    spectrum = _load_json(outdir, "spectrum_manifest.json")
    spectrum = spectrum if isinstance(spectrum, dict) else {}
    mode_file = outdir / "raman_modes.csv"
    settings = spectrum.get("settings")
    spectrum_files = spectrum.get("files")
    spectrum_current = (s.get("status") != "failed" and spectrum.get("status") == "complete"
                        and isinstance(settings, dict) and isinstance(spectrum_files, list)
                        and all(isinstance(name, str) for name in spectrum_files)
                        and mode_file.is_file() and spectrum.get("modes_sha256")
                        == hashlib.sha256(mode_file.read_bytes()).hexdigest())
    if spectrum_current:
        files.extend(["spectrum_manifest.json", *spectrum_files])
    if (outdir / "convergence" / "convergence.json").exists():
        files.extend(["convergence/convergence.json",
                      "convergence/convergence.csv",
                      "convergence/convergence.md"])
    if (outdir / "mode_assignment.csv").exists():
        files.extend(["mode_assignment.csv", "assignment.md",
                      "assignment_overlay.png"])
    lines.append("## 6. 输出文件")
    lines.append("")
    for f in files:
        mark = "（存在）" if (outdir / f).exists() else "（未生成）"
        lines.append("- `%s` %s" % (f, mark))
    if not spectrum_current and any(outdir.glob("spectrum_*")):
        lines.append("- 当前模式表没有可核验的谱图生成记录，请重新运行 `spectrum`。")
    lines.append("")

    struct = _load_json(outdir, "structure_manifest.json")
    effcfg = _load_json(outdir, "effective_config.json")
    run_state = s.get("run_state") or {}

    lines.append("## 7. 可复现信息")
    lines.append("")
    lines.append("| 项目 | 值 |")
    lines.append("|------|----|")
    lines.append("| QE 版本 | %s |" % (scf.get("qe_version") or "未知"))
    if struct:
        orig = struct.get("original") or {}
        ref = struct.get("refined") or {}
        lines.append("| 原始空间群 | #%s %s |" % (orig.get("spacegroup_number"), orig.get("spacegroup_symbol")))
        lines.append("| 精化空间群 | #%s %s |" % (ref.get("spacegroup_number"), ref.get("spacegroup_symbol")))
        lines.append("| 是否标准化 | %s |" % ("是" if struct.get("standardized") else "否"))
        if struct.get("fallback_reason"):
            lines.append("| 精化回退原因 | %s |" % struct["fallback_reason"])
        cp = struct.get("cellpar")
        if cp and len(cp) >= 6:
            lines.append("| 晶胞参数 | a=%.4f b=%.4f c=%.4f α=%.3f β=%.3f γ=%.3f |"
                         % (cp[0], cp[1], cp[2], cp[3], cp[4], cp[5]))
    for k in ("ecutwfc", "ecutrho", "conv_thr", "tr2_ph", "kpt_density"):
        if k in effcfg:
            lines.append("| %s | %s |" % (k, effcfg[k]))
    if spectrum_current:
        for k in ("spectrum_type", "spectrum_broadening", "fwhm", "temperature",
                  "laser_wavelength", "threshold"):
            if k in settings:
                lines.append("| %s | %s |" % (k, settings[k]))
    lines.append("| Raman (lraman) | %s |"
                 % ("开启" if effcfg.get("raman", s.get("lraman", True)) else "关闭（仅声子/IR）"))
    if ppman.get("xc_families"):
        lines.append("| 赝势 XC 家族 | %s |" % ", ".join(ppman["xc_families"]))
    if ppman.get("pseudo_source"):
        lines.append("| 赝势来源 | %s |" % "; ".join(ppman["pseudo_source"]))
    pps = ppman.get("pseudopotentials") or {}
    for elx in sorted(pps):
        e = pps[elx]
        lines.append("| UPF %s | %s (sha256 %s…) |"
                     % (elx, e.get("file"), (e.get("sha256") or "?")[:12]))
    stages = run_state.get("stages") or {}
    cmds = []
    for nm, rc in stages.items():
        if isinstance(rc, dict) and rc.get("command"):
            cmds.append("%s: %s" % (nm, " ".join(map(str, rc["command"]))))
    if cmds:
        lines.append("")
        lines.append("命令：")
        for c in cmds:
            lines.append("- `%s`" % c)
    lines.append("")

    conv = _load_convergence(outdir)
    if conv:
        lines.append("## 8. 收敛扫描")
        lines.append("")
        lines.append("来源：`convergence/convergence.json`（详见 "
                     "`convergence/convergence.md`）。")
        lines.append("")
        lines.append("| 参数组 | 参考值 | 推荐值 | 是否收敛 | 最大峰位偏差 (cm⁻¹) | 最大活性相对偏差 |")
        lines.append("|---|---|---|---|---|---|")
        for g in conv.get("groups") or []:
            if g.get("group") == "fwhm":
                lines.append("| FWHM | — | — | 后处理参数 | — | — |")
                continue
            fdev = act = None
            for v in g.get("variants") or []:
                if v.get("max_freq_dev") is not None:
                    fdev = max(fdev or 0.0, v["max_freq_dev"])
                if v.get("max_activity_rel") is not None:
                    act = max(act or 0.0, v["max_activity_rel"])
            lines.append("| %s | %s | %s | %s | %s | %s |" % (
                g.get("label") or g.get("group"),
                _fmt(g.get("reference_value")),
                _fmt(g.get("recommended")),
                "是" if g.get("converged") else "否",
                _fmt(fdev, 2),
                "—" if act is None else "%.2f%%" % (100 * act),
            ))
        lines.append("")

    assign = _load_assignment(outdir)
    if assign:
        counts: dict[str, int] = {}
        to_rows = [r for r in assign if (r.get("table") or "to") == "to"]
        for row in to_rows:
            counts[row.get("activity", "?")] = counts.get(row.get("activity", "?"), 0) + 1
        lines.append("## 9. 模式指认")
        lines.append("")
        lines.append("来源：`mode_assignment.csv` / `assignment.md`"
                     "（由 `assign` 子命令生成）。")
        lines.append("")
        lines.append("| 活性类型 | 模数 |")
        lines.append("|---|---|")
        labels = {"both": "Raman+IR", "raman": "仅 Raman", "ir": "仅 IR",
                  "silent": "静默", "soft": "软模（虚频）"}
        for key in ("both", "raman", "ir", "silent", "soft"):
            if counts.get(key):
                lines.append("| %s | %d |" % (labels[key], counts[key]))
        lines.append("")
        lines.append("| 模 | 频率 (cm⁻¹) | 活性 | Raman | IR | 退偏 | 简并 |")
        lines.append("|---|---|---|---|---|---|---|")
        for row in to_rows[:20]:
            depol = row.get("depol") or "—"
            lines.append("| %s | %s | %s | %s | %s | %s | %s |" % (
                row.get("mode", "—"), _fmt(float(row["freq_cm1"]), 2)
                if row.get("freq_cm1") else "—",
                labels.get(row.get("activity", ""), row.get("activity", "—")),
                _fmt(float(row["raman"]), 4) if row.get("raman") else "—",
                _fmt(float(row["ir"]), 4) if row.get("ir") else "—",
                depol if depol == "—" else _fmt(float(depol), 4),
                row.get("degeneracy", "1"),
            ))
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
