#!/usr/bin/env python3
"""Extract the bundled SSSP tarballs into the directories the skills expect.

The SSSP 1.3.0 efficiency/precision tarballs and JSON manifests are committed
here (compressed) to keep the repository manageable.  Run this once after
cloning to create ``sssp/upf/`` and ``sssp/upf_precision/``.

PseudoDojo NC-SR PBE UPF files are committed uncompressed in ``pseudodojo/``.
"""
from __future__ import annotations

import argparse
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SSSP = HERE / "sssp"
JOBS = [
    (SSSP / "SSSP_1.3.0_PBE_efficiency.tar.gz", SSSP / "upf"),
    (SSSP / "SSSP_1.3.0_PBE_precision.tar.gz", SSSP / "upf_precision"),
]


def _count_upf(dest: Path) -> int:
    return len([p for p in dest.glob("*") if p.suffix.lower() == ".upf"])


def _extract(tar_path: Path, dest: Path, force: bool) -> int:
    if not tar_path.is_file():
        raise FileNotFoundError(tar_path)
    if dest.is_dir() and _count_upf(dest) and not force:
        return _count_upf(dest)
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:*") as tar:
        tar.extractall(dest)
    return _count_upf(dest)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-extract even if present")
    args = ap.parse_args()
    for tar_path, dest in JOBS:
        n = _extract(tar_path, dest, args.force)
        print(f"{dest}: {n} UPF files")
    print(f"pseudodojo: {_count_upf(HERE / 'pseudodojo')} UPF files")
    print("done")


if __name__ == "__main__":
    main()
