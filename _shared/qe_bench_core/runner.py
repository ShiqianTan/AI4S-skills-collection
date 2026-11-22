"""pw.x execution helpers (CPU/GPU, MPI, job-status detection)."""

from __future__ import annotations

import shutil
import subprocess
import time
from pathlib import Path
from typing import Mapping


def prepare_tmp(directory: str | Path, name: str = "tmp") -> Path:
    """Recreate a clean QE scratch directory."""
    d = Path(directory) / name
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True, exist_ok=True)
    return d


def is_done(outfile: str | Path) -> bool:
    p = Path(outfile)
    return p.exists() and "JOB DONE" in p.read_text(errors="ignore")


def last_energy_ry(outfile: str | Path) -> float | None:
    energy = None
    for line in Path(outfile).read_text(errors="ignore").splitlines():
        if "!    total energy" in line and "Ry" in line:
            try:
                energy = float(line.split("=")[1].split("Ry")[0])
            except ValueError:
                pass
    return energy


def run_pw(
    inp: str | Path,
    out: str | Path,
    *,
    pw: str | Path,
    mpirun: str = "mpirun",
    np: int = 1,
    cwd: str | Path | None = None,
    env: Mapping[str, str] | None = None,
    gpu: int | str | None = None,
    timeout: int = 6 * 3600,
    bind_to_none: bool = True,
) -> dict:
    """Run ``pw.x`` and return exit code, wall time and completion flag."""
    inp = Path(inp)
    out = Path(out)
    cwd = Path(cwd) if cwd is not None else inp.parent
    cmd = [str(mpirun), "-np", str(np)]
    if bind_to_none:
        cmd += ["--bind-to", "none"]
    cmd += [str(pw), "-in", inp.name]
    run_env = dict(env or {})
    if gpu is not None:
        run_env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    t0 = time.time()
    with out.open("w") as fo, out.with_suffix(".err").open("w") as fe:
        proc = subprocess.run(
            cmd, cwd=cwd, env=run_env or None, stdout=fo, stderr=fe, timeout=timeout
        )
    return {
        "exit": proc.returncode,
        "wall_s": round(time.time() - t0, 1),
        "done": is_done(out),
        "output": str(out),
    }
