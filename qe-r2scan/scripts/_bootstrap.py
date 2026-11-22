"""Locate the shared ``qe_bench_core`` package and resolve skill paths."""
from __future__ import annotations

import sys
from pathlib import Path


def _find_root(start: Path) -> Path | None:
    for parent in [start, *start.parents]:
        if (parent / "_shared" / "qe_bench_core").is_dir():
            return parent
    return None


_REPO_ROOT = _find_root(Path(__file__).resolve())
if _REPO_ROOT is not None:
    sys.path.insert(0, str(_REPO_ROOT / "_shared"))

from qe_bench_core import paths as _paths  # noqa: E402


class SkillPaths:
    def __init__(self, config_path: str | None = None):
        self.cfg = _paths.load_config(config_path)
        self.repo_root = Path(self.cfg["_base_dir"])
        self.work = _paths.work_dir(self.cfg)
        self.results = _paths.results_dir(self.cfg)
        self.skill_root = Path(__file__).resolve().parent.parent

    def dataset(self, kind: str) -> Path:
        return _paths.dataset_path(self.cfg, kind)

    def pseudo_dir(self, kind: str = "sssp_efficiency") -> Path:
        return _paths.pseudo_dir(self.cfg, kind)

    def pseudo_json(self, kind: str = "efficiency") -> Path:
        return _paths.pseudo_json(self.cfg, kind)

    def qe_bin(self, kind: str = "cpu") -> Path:
        return _paths.qe_bin_dir(self.cfg, kind)

    def pw(self, kind: str = "cpu") -> Path:
        return self.qe_bin(kind) / "pw.x"

    def mpirun(self) -> str:
        return _paths.mpirun(self.cfg)

    def env(self, gpu: int | str | None = None) -> dict[str, str]:
        return _paths.build_env(self.cfg, gpu)


def get_paths(config_path: str | None = None) -> SkillPaths:
    return SkillPaths(config_path)


def add_config_arg(parser) -> None:
    parser.add_argument("--config", default=None, help="config.yaml (shared qe_bench_core)")
