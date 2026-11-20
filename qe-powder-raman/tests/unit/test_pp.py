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
    build_pseudopotential_manifest,
    check_pseudo,
    element_pp_map,
    parse_upf,
    raman_xc_issues,
)

PBE = PSEUDO / "pbe"


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


class TestRamanXc(unittest.TestCase):
    def test_pbe_nc_classified_gga(self):
        msgs = check_pseudo(PBE / "Li.pbe-fr-oncvpsp.UPF", expected_element="Li")
        self.assertFalse(any(lvl == "FAIL" for lvl, _k, _t in msgs))
        self.assertTrue(any(lvl == "PASS" and key == "xc" and "GGA" in text
                            for lvl, key, text in msgs))

    def test_pbe_manifest_metadata(self):
        m = build_pseudopotential_manifest(["Li", "F"], PBE, None)
        self.assertFalse(m["errors"])
        self.assertEqual(m["xc_families"], ["GGA"])
        self.assertFalse(m["raman_supported"])
        self.assertTrue(any("PseudoDojo" in s for s in m["pseudo_source"]))

    def test_pbe_blocked_for_raman(self):
        m = build_pseudopotential_manifest(["Li", "F"], PBE, None)
        issues = raman_xc_issues(m, raman_enabled=True)
        self.assertTrue(any(lvl == "FAIL" and key == "raman_gga"
                            for lvl, key, _t in issues))

    def test_pbe_allowed_without_raman(self):
        m = build_pseudopotential_manifest(["Li", "F"], PBE, None)
        issues = raman_xc_issues(m, raman_enabled=False)
        self.assertFalse(any(lvl == "FAIL" for lvl, _k, _t in issues))
        self.assertTrue(any(lvl == "INFO" for lvl, _k, _t in issues))

    def test_lda_supported_for_raman(self):
        m = build_pseudopotential_manifest(["Li", "F"], PSEUDO, None)
        issues = raman_xc_issues(m, raman_enabled=True)
        self.assertFalse(any(lvl == "FAIL" for lvl, _k, _t in issues))


class TestManifest(unittest.TestCase):
    def test_manifest_records_sha256(self):
        m = build_pseudopotential_manifest(["Li", "F"], PSEUDO, None)
        self.assertFalse(m["errors"])
        self.assertEqual(m["n_atoms"], 2)
        for el in ("Li", "F"):
            entry = m["pseudopotentials"][el]
            self.assertEqual(len(entry["sha256"]), 64)
        self.assertEqual(m["n_valence_electrons"], 8.0)


if __name__ == "__main__":
    unittest.main()
