"""Unit tests for qe_pp (UPF header parsing + validation + mapping)."""
import shutil
import sys
import tempfile
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
        self.assertEqual(len(m["Li"]), 1)
        self.assertEqual(len(m["F"]), 1)
        self.assertEqual(Path(m["Li"][0]).name, "Li.pz-hgh.UPF")
        self.assertEqual(Path(m["F"][0]).name, "F.pz-hgh.UPF")

    def test_missing(self):
        m = element_pp_map(PSEUDO, ["Uuo"])
        self.assertEqual(m, {"Uuo": []})

    def test_ambiguous(self):
        tmp = Path(tempfile.mkdtemp(prefix="qer_pp_"))
        try:
            shutil.copy(PSEUDO / "Li.pz-hgh.UPF", tmp / "Li_a.UPF")
            shutil.copy(PSEUDO / "Li.pz-hgh.UPF", tmp / "Li_b.UPF")
            m = element_pp_map(tmp, ["Li"])
            self.assertEqual(len(m["Li"]), 2)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class TestCheckPseudoElementMatch(unittest.TestCase):
    def test_element_mismatch(self):
        tmp = Path(tempfile.mkdtemp(prefix="qer_pp_"))
        try:
            dst = tmp / "Li.fake.UPF"
            shutil.copy(PSEUDO / "F.pz-hgh.UPF", dst)
            msgs = check_pseudo(dst, expected_element="Li")
            self.assertTrue(
                any(lvl == "FAIL" and key == "element_mismatch" for lvl, key, _ in msgs)
            )
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_element_match_ok(self):
        msgs = check_pseudo(PSEUDO / "Li.pz-hgh.UPF", expected_element="Li")
        self.assertFalse(
            any(lvl == "FAIL" and key == "element_mismatch" for lvl, key, _ in msgs)
        )


if __name__ == "__main__":
    unittest.main()
