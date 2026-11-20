#!/usr/bin/env python3
"""
Pseudopotential reading and validation for the qe-powder-raman skill.

  - parse a UPF header (element, pseudo_type, functional, z_valence)
  - check that the file is norm-conserving (NC) with LDA or GGA/PBE XC
  - build the element -> UPF manifest with provenance and SHA-256

Note: QE's ``lraman`` only supports LDA (``PHonon/PH/phq_setup.f90`` errors on
GGA). GGA/PBE NC pseudopotentials are therefore accepted for phonon/IR runs
(``raman: false`` / ``--no-lraman``) but rejected when Raman is requested.
"""

from __future__ import annotations

import fnmatch
import hashlib
import json
import re
from pathlib import Path

from qe_common import DATA_DIR, xc_family


def _sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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


def check_pseudo(
    path: str | Path, expected_element: str | None = None
) -> list[tuple[str, str, str]]:
    """Validate a single UPF for Raman compatibility. Returns (level,key,text).

    If `expected_element` is given, the UPF header element must match it
    (case-insensitive) or a FAIL ``element_mismatch`` is emitted.
    """
    msgs: list[tuple[str, str, str]] = []
    try:
        upf = parse_upf(path)
    except Exception as exc:  # noqa: BLE001
        return [("FAIL", "read", "cannot read pseudopotential %s: %s" % (path, exc))]

    tag = "%s (%s)" % (upf.element or "?", Path(path).name)

    if (
        expected_element
        and upf.element
        and upf.element.lower() != expected_element.lower()
    ):
        msgs.append(
            (
                "FAIL",
                "element_mismatch",
                "%s: UPF header element=%s does not match requested %s"
                % (Path(path).name, upf.element, expected_element),
            )
        )

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
    fam = xc_family(func)
    if fam == "LDA":
        msgs.append(("PASS", "xc", "%s functional=%s (LDA)" % (tag, func)))
    elif fam == "GGA":
        msgs.append(
            ("PASS", "xc",
             "%s functional=%s (GGA; usable for phonon/IR with lraman disabled)"
             % (tag, func))
        )
    else:
        msgs.append(
            ("FAIL", "xc",
             "%s functional=%s is not LDA/GGA; this skill requires NC pseudopotentials "
             "with LDA or PBE" % (tag, func))
        )

    if upf.z_valence is None:
        msgs.append(("WARN", "z_valence", "%s has no z_valence in header" % tag))
    else:
        msgs.append(("PASS", "z_valence", "%s z_valence=%.6g" % (tag, upf.z_valence)))
    return msgs


def raman_xc_issues(
    manifest: dict, raman_enabled: bool = True
) -> list[tuple[str, str, str]]:
    """Check whether a pseudopotential manifest can be used for the Raman step.

    QE's ``lraman`` implements only LDA third-order derivatives
    (``PHonon/PH/phq_setup.f90``: "third order derivatives not implemented with
    GGA").  GGA/PBE pseudopotentials are therefore accepted by
    :func:`check_pseudo` but must run with ``lraman=.false.`` (phonon/IR only).
    """
    msgs: list[tuple[str, str, str]] = []
    if not raman_enabled:
        msgs.append(("INFO", "lraman", "Raman disabled (lraman=.false.): phonon/IR only"))
        return msgs
    entries = (manifest.get("pseudopotentials") or {}).items()
    for el, entry in sorted(entries):
        func = entry.get("functional")
        fam = xc_family(func)
        if fam == "GGA":
            msgs.append(
                ("FAIL", "raman_gga",
                 "%s: functional=%s is GGA; QE lraman rejects GGA "
                 "('third order derivatives not implemented with GGA', "
                 "PHonon/PH/phq_setup.f90). Use LDA pseudopotentials or "
                 "--no-lraman for phonon/IR." % (el, func))
            )
        elif fam is None:
            msgs.append(
                ("FAIL", "raman_xc",
                 "%s: functional=%s is not recognized as LDA/GGA" % (el, func))
            )
    return msgs


def _pseudo_sources() -> dict[str, str]:
    """Filename-pattern -> provenance strings from data/pseudopotential_map.json."""
    path = DATA_DIR / "pseudopotential_map.json"
    try:
        data = json.loads(path.read_text())
    except Exception:  # noqa: BLE001
        return {}
    sources = data.get("sources") or {}
    return {str(k): str(v) for k, v in sources.items()}


def _source_for(file_name: str, sources: dict[str, str]) -> str:
    for pattern, source in sources.items():
        if fnmatch.fnmatch(file_name.lower(), pattern.lower()):
            return source
    return "user-supplied"


def element_pp_map(pseudo_dir: str | Path, elements: list[str]) -> dict[str, list[str]]:
    """Map requested elements to candidate UPF paths by UPF *content* (header
    ``element``), not filename.

    Returns ``{element: [path, ...]}``: an empty list means the element has no
    UPF, and more than one candidate means the mapping is ambiguous (callers
    must surface an actionable error rather than silently pick one).
    """
    d = Path(pseudo_dir)
    mapping: dict[str, list[str]] = {el: [] for el in elements}
    if not d.is_dir():
        return mapping
    indexed: dict[str, list[str]] = {}
    for p in sorted(list(d.glob("*.UPF")) + list(d.glob("*.upf"))):
        try:
            info = parse_upf(p)
        except Exception:  # noqa: BLE001 - skip unreadable files; caller surfaces missing
            continue
        el = (info.element or "").strip().lower()
        if el:
            indexed.setdefault(el, []).append(str(p))
    for el in elements:
        mapping[el] = indexed.get(el.lower(), [])
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
    sources = _pseudo_sources()
    zsum = 0.0
    for el in sorted(set(elements)):
        count = elements.count(el)
        cands = mapping.get(el, [])
        if not cands:
            manifest["errors"].append("no UPF for element %s in %s" % (el, pseudo_dir))
            continue
        if len(cands) > 1:
            manifest["errors"].append(
                "ambiguous UPF for element %s: %d candidates (%s)"
                % (el, len(cands), ", ".join(sorted(cands)))
            )
            continue
        path = cands[0]
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
            "sha256": _sha256(path),
            "source": _source_for(Path(path).name, sources),
        }
        manifest["pseudopotentials"][el] = entry
        if upf.z_valence is not None:
            zsum += float(upf.z_valence) * count
        for level, key, text in check_pseudo(path, expected_element=el):
            if level == "FAIL":
                manifest["errors"].append("%s: %s" % (el, text))
            elif level == "WARN":
                manifest["warnings"].append("%s: %s" % (el, text))
    manifest["n_valence_electrons"] = round(zsum, 6)
    families = sorted({
        fam for fam in (
            xc_family(e.get("functional"))
            for e in manifest["pseudopotentials"].values()
        ) if fam
    })
    manifest["xc_families"] = families
    manifest["raman_supported"] = bool(families) and "GGA" not in families
    manifest["pseudo_source"] = sorted({
        e.get("source") for e in manifest["pseudopotentials"].values()
        if e.get("source")
    })

    if outdir is not None:
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "pseudopotential_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
        )
    return manifest
