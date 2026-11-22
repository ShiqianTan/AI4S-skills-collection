"""Parity metrics used throughout the benchmark skills."""

from __future__ import annotations

from typing import Iterable

import numpy as np


def mae(values: Iterable[float]) -> float:
    return float(np.mean(np.abs(np.asarray(list(values), dtype=float))))


def rmse(values: Iterable[float]) -> float:
    return float(np.sqrt(np.mean(np.asarray(list(values), dtype=float) ** 2)))


def max_abs(values: Iterable[float]) -> float:
    return float(np.max(np.abs(np.asarray(list(values), dtype=float))))


def parity_metrics(reference, predicted) -> dict:
    ref = np.asarray(reference, dtype=float)
    pred = np.asarray(predicted, dtype=float)
    diff = pred - ref
    return {
        "n": int(ref.size),
        "MAE": mae(diff),
        "RMSE": rmse(diff),
        "bias": float(np.mean(diff)),
        "max_abs": max_abs(diff),
    }


def force_metrics(qe_forces, ref_forces) -> dict:
    qe = np.asarray(qe_forces, dtype=float)
    ref = np.asarray(ref_forces, dtype=float)
    diff = qe - ref
    return {
        "force_MAE_eV_A": float(np.mean(np.abs(diff))),
        "force_RMSE_eV_A": float(np.sqrt(np.mean(diff**2))),
        "force_max_err_eV_A": float(np.max(np.abs(diff))),
        "force_mae_norm_eV_A": float(
            np.mean(np.abs(np.linalg.norm(qe, axis=1) - np.linalg.norm(ref, axis=1)))
        ),
    }


def relative_energy_metrics(qe_energy_per_atom, ref_energy_per_atom) -> dict:
    qe = np.asarray(qe_energy_per_atom, dtype=float)
    ref = np.asarray(ref_energy_per_atom, dtype=float)
    rel_err = (qe - qe.min()) - (ref - ref.min())
    return {
        "energy_offset_eV_atom": float(np.mean(qe - ref)),
        "rel_energy_MAE_meV_atom": float(np.mean(np.abs(rel_err)) * 1000.0),
        "rel_energy_max_meV_atom": float(np.max(np.abs(rel_err)) * 1000.0),
    }


def group_parity_metrics(rows: list[dict], group_key: str, value_key: str) -> dict:
    """Aggregate a value column per material group."""
    out: dict[str, dict] = {}
    for key in sorted({r[group_key] for r in rows}):
        vals = [r[value_key] for r in rows if r[group_key] == key and r.get(value_key) is not None]
        if vals:
            out[key] = {"n": len(vals), "MAE": mae(vals), "RMSE": rmse(vals), "max_abs": max_abs(vals)}
    return out
