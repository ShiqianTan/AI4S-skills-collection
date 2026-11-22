#!/usr/bin/env python3
"""Apply the known Quantum ESPRESSO make.inc fixes for this build recipe."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

# Profiles correspond to the documented CPU / GPU / libxc build failures.
CPU_FFLAGS = "-fallow-argument-mismatch -fallow-invalid-boz"


def _find_var(text: str, var: str) -> re.Match | None:
    return re.search(rf"^(\s*{re.escape(var)}\s*=\s*)(.*)$", text, flags=re.MULTILINE)


def set_var(text: str, var: str, value: str) -> str:
    match = _find_var(text, var)
    if match:
        return text[: match.start()] + match.group(1) + value + text[match.end():]
    return text + f"\n{var} = {value}\n"


def append_var(text: str, var: str, value: str) -> str:
    match = _find_var(text, var)
    if not match:
        return text + f"\n{var} = {value}\n"
    current = match.group(2).strip()
    if value in current:
        return text
    new = f"{current} {value}".strip()
    return text[: match.start()] + match.group(1) + new + text[match.end():]


def apply_profile(text: str, args: argparse.Namespace) -> str:
    if args.profile in ("cpu", "cpu-libxc"):
        text = append_var(text, "FFLAGS", CPU_FFLAGS)
        text = set_var(text, "LD", args.mpif90)
        text = append_var(text, "LDFLAGS", "-fopenmp")
        if args.conda_prefix:
            text = set_var(
                text, "FFT_LIBS", f"-L{args.conda_prefix}/lib -lfftw3_omp -lfftw3"
            )
        if args.profile == "cpu-libxc":
            text = append_var(text, "DFLAGS", "-D__LIBXC")
    elif args.profile == "gpu":
        text = set_var(text, "LD", args.mpif90)
        if args.conda_prefix:
            text = set_var(
                text,
                "FFT_LIBS",
                f"-L{args.conda_prefix}/lib -lfftw3_threads -lfftw3",
            )
        if args.fft_include:
            text = set_var(text, "FFTW_INCLUDE", args.fft_include)
    else:  # pragma: no cover
        raise ValueError(args.profile)
    return text


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--make-inc", default="make.inc")
    ap.add_argument("--profile", choices=["cpu", "gpu", "cpu-libxc"], required=True)
    ap.add_argument("--mpif90", default="mpif90")
    ap.add_argument("--conda-prefix", default="")
    ap.add_argument("--fft-include", default="")
    args = ap.parse_args()

    path = Path(args.make_inc)
    text = path.read_text()
    updated = apply_profile(text, args)
    path.write_text(updated)
    print(f"patched {path} for profile '{args.profile}'")


if __name__ == "__main__":
    main()
