"""Unit tests for qe_pp (UPF header parsing + validation + mapping)."""
import sys
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
PSEUDO = Path(__file__).resolve().parents[2] / "pseudo"
sys.path.insert(0, str(SCRIPTS))

from qe_pp import (  # noqa: E402
    check_pseudo,
    element_pp_map,
    parse_upf,
)


class TestParseUpf(unittest.TestCase):
    def test_li_hgh(self):
        upf = parse_upf(PSEUDO / "Li.pz-hgh.UPF")
        self.assertEqual(upf.element, "Li")
        self.assertEqual(upf.pseudo_type, "NC")
        self.assertEqual(upf.functional, "SLA-PZ-NOGX-NOGC")
        self.assertEqual(upf.z_valence, 1.0)
        self.assertFalse(upf.is_ultrasoft)
        self.assertFalse(upf.is_paw)

    def test_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            parse_upf(PSEUDO / "NoSuchElement.pz-hgh.UPF")


class TestCheckPseudo(unittest.TestCase):
    def test_li_ok(self):
        msgs = check_pseudo(PSEUDO / "Li.pz-hgh.UPF")
        self.assertFalse(any(lvl == "FAIL" for lvl, _, _ in msgs))
        self.assertTrue(any(lvl == "PASS" and key == "type" for lvl, key, _ in msgs))
        self.assertTrue(any(lvl == "PASS" and key == "xc" for lvl, key, _ in msgs))


class TestElementPpMap(unittest.TestCase):
    def test_maps(self):
        m = element_pp_map(PSEUDO, ["Li", "F"])
        self.assertIn("Li", m)
        self.assertIn("F", m)
        self.assertIn("Li.pz-hgh.UPF", m["Li"])
        self.assertIn("F.pz-hgh.UPF", m["F"])

    def test_missing(self):
        m = element_pp_map(PSEUDO, ["Uuo"])
        self.assertEqual(m, {})


if __name__ == "__main__":
    unittest.main()
