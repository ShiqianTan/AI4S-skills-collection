#!/usr/bin/env python3
"""
Light structural/syntactic validation of a generated QE input file.

Detects pw.x / ph.x / dynmat.x inputs and checks:
  - balanced namelists (& ... /)
  - required keywords present per input type
  - pw.x: ATOMIC_SPECIES / CELL_PARAMETERS / ATOMIC_POSITIONS / K_POINTS present,
    ntyp == number of species lines, nat == number of position lines,
    occupations must be 'fixed' for the lraman MVP
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def _detect(text: str) -> str:
    if re.search(r"&inputph", text, re.I):
        return "ph"
    if re.search(r"&\s*input\b", text, re.I):
        return "dynmat"
    if re.search(r"&control|&system|&electrons", text, re.I):
        return "pw"
    return "unknown"


def _section_lines(text: str, start_marker: str, end_markers: tuple[str, ...]) -> list[str]:
    lines = text.splitlines()
    out: list[str] = []
    started = False
    for ln in lines:
        s = ln.strip()
        if not started:
            if s.upper().startswith(start_marker):
                started = True
            continue
        if s.startswith("&") or s == "/":
            break
        if any(s.upper().startswith(m) for m in end_markers):
            break
        if s:
            out.append(s)
    return out


def _validate_pw(text: str) -> list[tuple[str, str, str]]:
    msgs: list[tuple[str, str, str]] = []

    def req(present, key, label):
        msgs.append(("PASS" if present else "FAIL", key, label))

    n_open = len(re.findall(r"&(?:control|system|electrons|ions|cell)\b", text, re.I))
    n_close = len(re.findall(r"^\s*/\s*$", text, re.M))
    if n_open != n_close:
        msgs.append(("FAIL", "namelist", "unbalanced namelists: %d open, %d close" % (n_open, n_close)))

    req(bool(re.search(r"\bnat\b\s*=", text, re.I)), "nat", "nat present")
    req(bool(re.search(r"\bntyp\b\s*=", text, re.I)), "ntyp", "ntyp present")
    req(bool(re.search(r"\becutwfc\b\s*=", text, re.I)), "ecutwfc", "ecutwfc present")
    req(bool(re.search(r"\boccupations\b\s*=", text, re.I)), "occupations", "occupations present")
    req(bool(re.search(r"ATOMIC_SPECIES", text)), "atomic_species", "ATOMIC_SPECIES present")
    req(bool(re.search(r"CELL_PARAMETERS", text)), "cell", "CELL_PARAMETERS present")
    req(bool(re.search(r"ATOMIC_POSITIONS", text)), "positions", "ATOMIC_POSITIONS present")
    req(bool(re.search(r"K_POINTS", text)), "kpoints", "K_POINTS present")

    m = re.search(r"occupations\s*=\s*['\"]?(\w+)", text, re.I)
    if m and m.group(1).lower() != "fixed":
        msgs.append(("FAIL", "occupations",
                     "occupations must be 'fixed' for the Raman MVP (got '%s')" % m.group(1)))

    ntyp = re.search(r"\bntyp\b\s*=\s*(\d+)", text, re.I)
    nat = re.search(r"\bnat\b\s*=\s*(\d+)", text, re.I)
    species = _section_lines(text, "ATOMIC_SPECIES", ("CELL_PARAMETERS",))
    positions = _section_lines(text, "ATOMIC_POSITIONS", ("K_POINTS",))
    if ntyp and len(species) != int(ntyp.group(1)):
        msgs.append(("FAIL", "ntyp", "ntyp=%s but %d ATOMIC_SPECIES lines" % (ntyp.group(1), len(species))))
    if nat and len(positions) != int(nat.group(1)):
        msgs.append(("FAIL", "nat", "nat=%s but %d ATOMIC_POSITIONS lines" % (nat.group(1), len(positions))))
    return msgs


def _validate_ph(text: str) -> list[tuple[str, str, str]]:
    msgs: list[tuple[str, str, str]] = []

    def req(present, key, label):
        msgs.append(("PASS" if present else "FAIL", key, label))

    req(bool(re.search(r"&inputph", text, re.I)), "namelist", "&inputph present")
    req(bool(re.search(r"\btr2_ph\b\s*=", text, re.I)), "tr2_ph", "tr2_ph present")
    req(bool(re.search(r"\bepsil\b\s*=\s*\.true\.", text, re.I)), "epsil", "epsil=.true.")
    req(bool(re.search(r"\btrans\b\s*=\s*\.true\.", text, re.I)), "trans", "trans=.true.")
    req(bool(re.search(r"\blraman\b\s*=\s*\.true\.", text, re.I)), "lraman", "lraman=.true.")
    req(bool(re.search(r"\basr\b\s*=\s*\.true\.", text, re.I)), "asr", "asr=.true.")
    req(bool(re.search(r"\bfildyn\b\s*=", text, re.I)), "fildyn", "fildyn present")
    qline = re.search(r"/\s*\n\s*([-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+)", text)
    if qline:
        msgs.append(("PASS", "q", "q-point line present: %s" % qline.group(1).strip()))
    else:
        msgs.append(("FAIL", "q", "missing q-point line after &inputph"))
    return msgs


def _validate_dynmat(text: str) -> list[tuple[str, str, str]]:
    msgs: list[tuple[str, str, str]] = []

    def req(present, key, label):
        msgs.append(("PASS" if present else "FAIL", key, label))

    req(bool(re.search(r"&\s*input\b", text, re.I)), "namelist", "&input present")
    req(bool(re.search(r"\bfildyn\b\s*=", text, re.I)), "fildyn", "fildyn present")
    req(bool(re.search(r"\basr\b\s*=\s*['\"]crystal['\"]", text, re.I)), "asr", "asr='crystal'")
    return msgs


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("usage: validate_input.py FILE", file=sys.stderr)
        return 2
    path = Path(args[0])
    if not path.exists():
        print("validate_input: not found: %s" % path, file=sys.stderr)
        return 1

    text = path.read_text()
    typ = _detect(text)
    print("validate_input: %s (type: %s)" % (path.name, typ))

    if typ == "pw":
        msgs = _validate_pw(text)
    elif typ == "ph":
        msgs = _validate_ph(text)
    elif typ == "dynmat":
        msgs = _validate_dynmat(text)
    else:
        msgs = [("FAIL", "type", "could not detect input type (expected pw.x/ph.x/dynmat.x)")]

    n_fail = 0
    for level, key, label in msgs:
        print("  [%s] %-14s %s" % (level, key, label))
        if level == "FAIL":
            n_fail += 1

    if n_fail:
        print("validate_input: FAIL (%d issue(s))" % n_fail)
        return 1
    print("validate_input: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
