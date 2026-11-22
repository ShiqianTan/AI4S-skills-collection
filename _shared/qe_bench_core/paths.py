"""Configuration and path resolution for the QE benchmark skills.

Resolution priority (highest first):

1. explicit ``--config`` file passed to :func:`load_config`
2. ``QE_BENCH_CONFIG`` environment variable
3. environment variables listed in :data:`ENV_MAP`
4. a ``config.yaml`` next to the current working directory
5. built-in defaults

The original reproduction scripts hard-coded absolute paths.  Every packaged
script calls :func:`get_config` and reads paths from the returned mapping so the
same code runs on another machine by editing a single ``config.yaml``.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Mapping

try:  # PyYAML is optional; a JSON config also works.
    import yaml
except Exception:  # pragma: no cover - yaml is in every skill requirement
    yaml = None

DEFAULTS: dict[str, Any] = {
    "root": ".",
    "qe": {
        "bin_dir": "",
        "gpu_bin_dir": "",
        "libxc_bin_dir": "",
        "mpirun": "mpirun",
        "nvhpc_root": "",
        "conda_lib": "",
        "mpi_np": 1,
        "omp_num_threads": 8,
    },
    "pseudo": {
        "sssp_efficiency": "pseudos/sssp/upf",
        "sssp_precision": "pseudos/sssp/upf_precision",
        "sssp_efficiency_json": "pseudos/sssp/SSSP_1.3.0_PBE_efficiency.json",
        "sssp_precision_json": "pseudos/sssp/SSSP_1.3.0_PBE_precision.json",
        "pseudodojo": "pseudos/pseudodojo",
    },
    "datasets": {
        "mp20": "",
        "mptraj": "",
        "salex": "",
        "matpes": "",
        "wbm": "",
        "wbm_is2re": "",
    },
    "output": {
        "work_dir": "runs",
        "results_dir": "runs/results",
    },
}

# environment variable -> (section, key)
ENV_MAP: dict[str, tuple[str, str]] = {
    "QE_BIN_DIR": ("qe", "bin_dir"),
    "QE_GPU_BIN_DIR": ("qe", "gpu_bin_dir"),
    "QE_LIBXC_BIN_DIR": ("qe", "libxc_bin_dir"),
    "QE_MPIRUN": ("qe", "mpirun"),
    "QE_NVHPC_ROOT": ("qe", "nvhpc_root"),
    "QE_CONDA_LIB": ("qe", "conda_lib"),
    "QE_DATASETS_MP20": ("datasets", "mp20"),
    "QE_DATASETS_MPTRAJ": ("datasets", "mptraj"),
    "QE_DATASETS_SALEX": ("datasets", "salex"),
    "QE_DATASETS_MATPES": ("datasets", "matpes"),
    "QE_DATASETS_WBM": ("datasets", "wbm"),
    "QE_DATASETS_WBM_IS2RE": ("datasets", "wbm_is2re"),
    "QE_PSEUDO_SSSP": ("pseudo", "sssp_efficiency"),
    "QE_PSEUDO_SSSP_PRECISION": ("pseudo", "sssp_precision"),
    "QE_PSEUDO_PSEUDODOJO": ("pseudo", "pseudodojo"),
    "QE_RESULTS_DIR": ("output", "results_dir"),
    "QE_WORK_DIR": ("output", "work_dir"),
}


def find_repo_root(start: str | Path | None = None) -> Path:
    """Return the repository root (the directory containing ``_shared``)."""
    cur = Path(start or __file__).resolve()
    if cur.is_file():
        cur = cur.parent
    for parent in [cur, *cur.parents]:
        if (parent / "_shared").is_dir():
            return parent
    return Path.cwd()


def _deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    out = {k: (dict(v) if isinstance(v, dict) else v) for k, v in base.items()}
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(out.get(key), Mapping):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def _read_config_file(path: Path) -> dict[str, Any]:
    text = path.read_text()
    if path.suffix in (".json",):
        import json

        return json.loads(text)
    if yaml is None:  # pragma: no cover
        raise RuntimeError("PyYAML is required to read YAML config files")
    return yaml.safe_load(text) or {}


def _apply_env(cfg: dict[str, Any]) -> dict[str, Any]:
    for env_name, (section, key) in ENV_MAP.items():
        value = os.environ.get(env_name)
        if value:
            cfg.setdefault(section, {})[key] = value
    return cfg


def load_config(config_path: str | Path | None = None) -> dict[str, Any]:
    """Load and merge configuration from file, environment and defaults."""
    cfg = _deep_merge(DEFAULTS, {})
    base = Path.cwd()

    candidate = config_path or os.environ.get("QE_BENCH_CONFIG")
    if candidate is None and (base / "config.yaml").is_file():
        candidate = base / "config.yaml"
    if candidate is not None:
        path = Path(candidate).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"config file not found: {path}")
        cfg = _deep_merge(cfg, _read_config_file(path))
        base = path.resolve().parent

    cfg = _apply_env(cfg)
    cfg["_base_dir"] = str(base)
    return cfg


def get_config(config_path: str | Path | None = None) -> dict[str, Any]:
    """Alias kept for readability in scripts."""
    return load_config(config_path)


def resolve_path(value: str | Path | None, cfg: Mapping[str, Any]) -> Path | None:
    """Resolve a config path relative to ``root``/config directory."""
    if value in (None, ""):
        return None
    p = Path(str(value)).expanduser()
    if p.is_absolute():
        return p
    root = Path(str(cfg.get("root") or ".")).expanduser()
    if not root.is_absolute():
        root = Path(str(cfg.get("_base_dir", "."))) / root
    return (root / p).resolve()


def qe_bin_dir(cfg: Mapping[str, Any], kind: str = "cpu") -> Path:
    key = {"cpu": "bin_dir", "gpu": "gpu_bin_dir", "libxc": "libxc_bin_dir"}[kind]
    path = resolve_path(cfg.get("qe", {}).get(key), cfg)
    if path is None:
        raise ValueError(
            f"qe.{key} is not configured; set it in config.yaml or QE_*_BIN_DIR"
        )
    return path


def qe_executable(cfg: Mapping[str, Any], program: str = "pw.x", kind: str = "cpu") -> Path:
    return qe_bin_dir(cfg, kind) / program


def mpirun(cfg: Mapping[str, Any]) -> str:
    return str(cfg.get("qe", {}).get("mpirun") or "mpirun")


def pseudo_dir(cfg: Mapping[str, Any], kind: str = "sssp_efficiency") -> Path:
    path = resolve_path(cfg.get("pseudo", {}).get(kind), cfg)
    if path is None:
        raise ValueError(f"pseudo.{kind} is not configured")
    return path


def pseudo_json(cfg: Mapping[str, Any], kind: str = "efficiency") -> Path:
    key = f"sssp_{kind}_json"
    path = resolve_path(cfg.get("pseudo", {}).get(key), cfg)
    if path is None:
        raise ValueError(f"pseudo.{key} is not configured")
    return path


def dataset_path(cfg: Mapping[str, Any], kind: str) -> Path:
    path = resolve_path(cfg.get("datasets", {}).get(kind), cfg)
    if path is None:
        raise ValueError(
            f"datasets.{kind} is not configured; set it in config.yaml or "
            f"QE_DATASETS_{kind.upper()}"
        )
    return path


def work_dir(cfg: Mapping[str, Any], subdir: str | None = None) -> Path:
    base = resolve_path(cfg.get("output", {}).get("work_dir") or "runs", cfg)
    assert base is not None
    return base if subdir is None else base / subdir


def results_dir(cfg: Mapping[str, Any]) -> Path:
    path = resolve_path(cfg.get("output", {}).get("results_dir"), cfg)
    assert path is not None
    path.mkdir(parents=True, exist_ok=True)
    return path


def build_env(cfg: Mapping[str, Any], gpu: int | str | None = None) -> dict[str, str]:
    """Build a subprocess environment for pw.x (CPU or GPU build)."""
    env = os.environ.copy()
    qe = cfg.get("qe", {})
    nvhpc = resolve_path(qe.get("nvhpc_root"), cfg)
    conda_lib = resolve_path(qe.get("conda_lib"), cfg)
    if nvhpc is not None:
        nv = nvhpc
        env["PATH"] = (
            f"{nv}/compilers/bin:{nv}/comm_libs/openmpi4/bin:"
            f"{nv}/cuda/12.9/bin:" + env.get("PATH", "")
        )
        libs = [str(nv / "compilers/lib"), str(nv / "comm_libs/openmpi4/lib"),
                str(nv / "cuda/12.9/lib64"),
                str(nv / "math_libs/12.9/targets/x86_64-linux/lib")]
        if conda_lib is not None:
            libs.insert(0, str(conda_lib))
        env["LD_LIBRARY_PATH"] = ":".join(libs)
        env["NVCOMPILERS"] = str(nv / "compilers")
        env["CUDA_HOME"] = str(nv / "cuda/12.9")
    env["OMPI_ALLOW_RUN_AS_ROOT"] = "1"
    env["OMPI_ALLOW_RUN_AS_ROOT_CONFIRM"] = "1"
    env["OMP_NUM_THREADS"] = str(qe.get("omp_num_threads", 8))
    if gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    return env
