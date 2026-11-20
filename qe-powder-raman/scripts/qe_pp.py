#!/usr/bin/env python3
"""
Pseudopotential reading and validation for the qe-powder-raman skill.

Step 0 scope:
  - parse a UPF header (element, pseudo_type, functional, z_valence)
  - check that the file is norm-conserving (NC) and LDA (PZ family)

Step 1 will extend this with a full element -> UPF manifest and the
bundled `pseudo/` directory mapping.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from qe_common import ACCEPTED_XC_LDA


class UPFInfo:
    def __init__(self) -> None:
        self.path: Path | None = None
        self.element: str | None = None
        self.pseudo_type: str | None = None
        self.functional: str | None = None
        self.z_valence: float | None = None
        self.is_ultrasoft: bool | None = None
        self.is_paw: bool | None = None


def parse_upf(path: str | Path) -> UPFInfo:
    """Parse a UPF file header into a UPFInfo. Raises on unreadable file."""
    info = UPFInfo()
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"pseudopotential not found: {p}")
    text = p.read_text(errors="replace")
    info.path = p

    def _attr(name: str) -> str | None:
        m = re.search(r'%s\s*=\s*"([^"]*)"' % re.escape(name), text)
        return m.group(1).strip() if m else None

    info.element = _attr("element")
    info.pseudo_type = _attr("pseudo_type")
    info.functional = _attr("functional")
    info.is_ultrasoft = (_attr("is_ultrasoft") or "").strip().upper() == "T"
    info.is_paw = (_attr("is_paw") or "").strip().upper() == "T"

    zv = _attr("z_valence")
    if zv is not None:
        try:
            info.z_valence = float(zv)
        except ValueError:
            info.z_valence = None

    # Fallback: the human-readable <PP_INFO> block sometimes carries these.
    if info.pseudo_type is None:
        m = re.search(r"Pseudopotential type:\s*(\S+)", text)
        if m:
            info.pseudo_type = m.group(1)
    if info.functional is None:
        m = re.search(r"Functional:\s*(\S+)", text)
        if m:
            info.functional = m.group(1)
    if info.element is None:
        m = re.search(r"Element:\s*(\S+)", text)
        if m:
            info.element = m.group(1)
    return info


def check_pseudo(path: str | Path) -> list[tuple[str, str, str]]:
    """Validate a single UPF for Raman compatibility. Returns (level,key,text)."""
    msgs: list[tuple[str, str, str]] = []
    try:
        upf = parse_upf(path)
    except Exception as exc:  # noqa: BLE001
        return [("FAIL", "read", "cannot read pseudopotential %s: %s" % (path, exc))]

    tag = "%s (%s)" % (upf.element or "?", Path(path).name)

    if upf.is_paw:
        msgs.append(("FAIL", "paw", "%s is PAW; lraman requires norm-conserving (NC)" % tag))
    if upf.is_ultrasoft or (upf.pseudo_type or "").upper() in ("US", "USPP"):
        msgs.append(("FAIL", "uspp", "%s is ultrasoft; lraman requires NC" % tag))

    pt = (upf.pseudo_type or "").upper()
    if pt == "NC" and not upf.is_ultrasoft and not upf.is_paw:
        msgs.append(("PASS", "type", "%s is norm-conserving (NC)" % tag))
    elif pt not in ("NC", "US", "USPP", "PAW"):
        msgs.append(("WARN", "type", "%s pseudo_type=%r unrecognized" % (tag, upf.pseudo_type)))

    func = upf.functional or ""
    if func in ACCEPTED_XC_LDA:
        msgs.append(("PASS", "xc", "%s functional=%s (LDA)" % (tag, func)))
    elif "PZ" in func or "LDA" in func.upper():
        msgs.append(("PASS", "xc", "%s functional=%s (LDA)" % (tag, func)))
    else:
        msgs.append(
            ("FAIL", "xc", "%s functional=%s is not LDA; Raman skill requires LDA+NC" % (tag, func))
        )

    if upf.z_valence is None:
        msgs.append(("WARN", "z_valence", "%s has no z_valence in header" % tag))
    else:
        msgs.append(("PASS", "z_valence", "%s z_valence=%.6g" % (tag, upf.z_valence)))
    return msgs


def element_pp_map(pseudo_dir: str | Path, elements: list[str]) -> dict[str, str]:
    """Map element symbols to UPF files in pseudo_dir (name <El>.*.UPF)."""
    d = Path(pseudo_dir)
    mapping: dict[str, str] = {}
    if not d.is_dir():
        return mapping
    files = sorted(p.name for p in d.glob("*.UPF")) + sorted(p.name for p in d.glob("*.upf"))
    for el in elements:
        cand = [f for f in files if f.split(".")[0].lower() == el.lower()]
        if cand:
            mapping[el] = str(d / cand[0])
    return mapping


def build_pseudopotential_manifest(
    elements: list[str],
    pseudo_dir: str | Path,
    outdir: str | Path | None = None,
) -> dict:
    """Map a structure's elements to validated UPF files.

    Writes `pseudopotential_manifest.json` (if outdir is given) recording each
    element's UPF header fields (element / pseudo_type / functional / z_valence)
    plus any warnings or hard errors (PAW / ultrasoft / non-LDA / missing).
    """
    manifest: dict = {
        "generator": "qe-powder-raman",
        "pseudo_dir": str(Path(pseudo_dir)),
        "elements": sorted(set(elements)),
        "n_atoms": len(elements),
        "pseudopotentials": {},
        "warnings": [],
        "errors": [],
    }
    mapping = element_pp_map(pseudo_dir, elements)
    zsum = 0.0
    for el in sorted(set(elements)):
        count = elements.count(el)
        if el not in mapping:
            manifest["errors"].append("no UPF for element %s in %s" % (el, pseudo_dir))
            continue
        path = mapping[el]
        upf = parse_upf(path)
        entry = {
            "file": Path(path).name,
            "element": upf.element or el,
            "pseudo_type": upf.pseudo_type,
            "functional": upf.functional,
            "z_valence": upf.z_valence,
            "is_paw": upf.is_paw,
            "is_ultrasoft": upf.is_ultrasoft,
            "count": count,
        }
        manifest["pseudopotentials"][el] = entry
        if upf.z_valence is not None:
            zsum += float(upf.z_valence) * count
        for level, key, text in check_pseudo(path):
            if level == "FAIL":
                manifest["errors"].append("%s: %s" % (el, text))
            elif level == "WARN":
                manifest["warnings"].append("%s: %s" % (el, text))
    manifest["n_valence_electrons"] = round(zsum, 6)

    if outdir is not None:
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "pseudopotential_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
        )
    return manifest
