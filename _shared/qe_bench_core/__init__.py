"""Shared utilities for the Quantum ESPRESSO VASP-parity skills.

The modules in this package are intentionally dependency-light and are imported
by every ``qe-*`` skill in the collection.  Import the submodules you need::

    from qe_bench_core import paths, qeio, runner, parse, metrics

No module performs I/O at import time, so importing the package is always safe.
"""

from __future__ import annotations

__all__ = [
    "paths",
    "structure",
    "qeio",
    "runner",
    "parse",
    "metrics",
    "mp2020",
    "mprelax",
    "plotting",
    "pseudos",
]

__version__ = "0.1.0"
