#!/usr/bin/env python3
"""Run the qe-powder-raman test suite (unit + integration, no QE/MPI needed)."""
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def main() -> int:
    loader = unittest.TestLoader()
    suite = loader.discover(str(HERE), pattern="test_*.py")
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
