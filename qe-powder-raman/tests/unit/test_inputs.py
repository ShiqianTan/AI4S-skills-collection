"""Unit tests for qe_inputs (template rendering + helpers)."""
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
PSEUDO = Path(__file__).resolve().parents[2] / "pseudo"
sys.path.insert(0, str(SCRIPTS))

from ase import Atoms  # noqa: E402

from qe_inputs import (  # noqa: E402
    fmt_num,
    generate_inputs,
    kpoint_grid,
    render,
    species_order,
)


class TestFmtNum(unittest.TestCase):
    def test_bool(self):
        self.assertEqual(fmt_num(True), ".true.")
        self.assertEqual(fmt_num(False), ".false.")

    def test_int(self):
        self.assertEqual(fmt_num(3), "3")

    def test_float_integer(self):
        self.assertEqual(fmt_num(60.0), "60")
        self.assertEqual(fmt_num(1.0), "1")
        self.assertEqual(fmt_num(0.0), "0")

    def test_sci_notation(self):
        self.assertEqual(fmt_num(1e-10), "1.0e-10")
        self.assertEqual(fmt_num(1e-14), "1.0e-14")

    def test_regular_float(self):
        self.assertEqual(fmt_num(0.5), "0.5")


class TestRender(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(render("prefix = {prefix}", {"prefix": "abc"}), "prefix = abc")

    def test_missing(self):
        self.assertEqual(render("a {x} b", {}), "a  b")

    def test_untouched_braces(self):
        # bash ${...} must not be touched
        self.assertEqual(render("x ${VAR} y", {"VAR": "ignored"}), "x ${VAR} y")


class TestSpeciesOrder(unittest.TestCase):
    def test_first_appearance(self):
        atoms = Atoms("HOH", positions=[[0, 0, 0], [0, 0, 1], [1, 1, 1]])
        self.assertEqual(species_order(atoms), ["H", "O"])


class TestKpointGrid(unittest.TestCase):
    def test_cubic(self):
        cell = [[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]]
        k = kpoint_grid(cell, density=30.0)
        self.assertEqual(k, (6, 6, 6))

    def test_anisotropic(self):
        cell = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 15.0]]
        k = kpoint_grid(cell, density=30.0)
        self.assertEqual(k, (10, 10, 2))

    def test_min_one(self):
        cell = [[20.0, 0, 0], [0, 20.0, 0], [0, 0, 20.0]]
        k = kpoint_grid(cell, density=30.0)
        self.assertEqual(k, (2, 2, 2))


class TestLoValidation(unittest.TestCase):
    def _cfg(self):
        from qe_common import load_defaults
        return load_defaults()

    def _atoms(self):
        return Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]],
                     cell=[5, 5, 5], pbc=True)

    def test_zero_direction_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                generate_inputs(self._atoms(), self._cfg(), d, PSEUDO, "t",
                                lo_dirs=[(0.0, 0.0, 0.0)])

    def test_non_finite_direction_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                generate_inputs(self._atoms(), self._cfg(), d, PSEUDO, "t",
                                lo_dirs=[(float("nan"), 0.0, 1.0)])

    def test_valid_direction_written(self):
        with tempfile.TemporaryDirectory() as d:
            generate_inputs(self._atoms(), self._cfg(), d, PSEUDO, "t",
                            lo_dirs=[(0.0, 0.0, 1.0)])
            text = (Path(d) / "dynmat_lo.in").read_text()
            self.assertIn("= 0 0 1", text)


if __name__ == "__main__":
    unittest.main()
