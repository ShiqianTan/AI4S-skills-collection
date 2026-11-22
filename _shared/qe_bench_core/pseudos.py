"""Pseudopotential bookkeeping: SSSP manifests, PseudoDojo, valences, unpacking."""

from __future__ import annotations

import json
import re
import tarfile
from pathlib import Path


def load_sssp(json_path: str | Path) -> dict:
    return json.loads(Path(json_path).read_text())


def load_pseudodojo(json_path: str | Path | None = None) -> dict | None:
    if json_path is None:
        return None
    return json.loads(Path(json_path).read_text())


def z_valence(upf_path: str | Path) -> float:
    text = Path(upf_path).read_text(errors="ignore")
    match = re.search(r'z_valence\s*=\s*"\s*([\d.eE+-]+)\s*"', text)
    if match is None:
        match = re.search(r"([\d.eE+-]+)\s+Z valence", text)
    if match is None:
        raise RuntimeError(f"cannot read z_valence from {upf_path}")
    return float(match.group(1))


def sssp_map(sssp: dict, elements: list[str], upf_dir: str | Path) -> dict:
    upf_dir = Path(upf_dir)
    mapping = {}
    for el in elements:
        rec = sssp[el]
        mapping[el] = {
            "filename": rec["filename"],
            "cutoff_wfc": float(rec["cutoff_wfc"]),
            "cutoff_rho": float(rec["cutoff_rho"]),
            "z_valence": z_valence(upf_dir / rec["filename"]),
        }
    return mapping


def pseudodojo_map(elements: list[str], upf_dir: str | Path) -> dict:
    upf_dir = Path(upf_dir)
    mapping = {}
    for el in elements:
        path = upf_dir / f"{el}.upf"
        mapping[el] = {
            "filename": path.name,
            "cutoff_wfc": 60.0,
            "cutoff_rho": 480.0,
            "z_valence": z_valence(path),
        }
    return mapping


def n_electrons(structure, pmap: dict) -> float:
    return float(
        sum(
            pmap[str(el)]["z_valence"] * structure.composition[el]
            for el in structure.composition.elements
        )
    )


def extract_tarball(tar_path: str | Path, dest_dir: str | Path) -> Path:
    """Extract an SSSP/PseudoDojo tarball into ``dest_dir``."""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:*") as tar:
        tar.extractall(dest_dir)
    return dest_dir
