#!/usr/bin/env python3
"""
List available QE input templates (relax / scf / ph / dynmat).

Step 0: placeholder. Templates are populated in Step 2.
"""

from __future__ import annotations

import sys
from pathlib import Path

from qe_common import TEMPLATE_DIR


def main(argv: list[str] | None = None) -> int:
    templates = sorted(TEMPLATE_DIR.glob("*.tmpl")) if TEMPLATE_DIR.is_dir() else []
    if not templates:
        print("no templates yet (populated in Step 2) under %s" % TEMPLATE_DIR)
        return 0
    for t in templates:
        print(t.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
