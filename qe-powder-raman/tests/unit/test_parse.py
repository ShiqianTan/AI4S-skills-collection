"""Unit tests for qe_parse (low-level parsers + statistics)."""
import sys
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
sys.path.insert(0, str(SCRIPTS))

from qe_parse import (  # noqa: E402
    _stats,
    parse_asr,
    parse_dielectric,
    parse_dynmat,
    parse_scf,
)


def _read(system, name):
    return (FIXTURES / system / name).read_text(errors="replace")


class TestParseScf(unittest.TestCase):
    def test_lif_scf(self):
        d = parse_scf(_read("LiF_mp-1138", "scf.out"))
        self.assertEqual(d["nelec"], 8.0)
        self.assertEqual(d["nbnd"], 4)
        self.assertEqual(d["nk"], 32)
        self.assertEqual(d["ecutwfc"], 60.0)
        self.assertEqual(d["sym_ops"], 12)
        self.assertEqual(d["scf_iter"], 20)
        self.assertAlmostEqual(d["volume_a3"], 17.022, places=3)
        self.assertAlmostEqual(d["etot_ry"], -47.47899621, places=6)

    def test_alas_scf(self):
        d = parse_scf(_read("AlAs_mp-2172", "scf.out"))
        self.assertEqual(d["nelec"], 8.0)
        self.assertEqual(d["sym_ops"], 6)  # rhombohedral primitive cell
        self.assertAlmostEqual(d["etot_ry"], -17.02810451, places=6)


class TestParseDynmat(unittest.TestCase):
    def test_alas(self):
        d = parse_dynmat(_read("AlAs_mp-2172", "dynmat.out"))
        self.assertAlmostEqual(d["cm_factor"], 0.210866, places=5)
        self.assertAlmostEqual(d["raman_cm_factor"], 0.044465, places=5)
        # polarizability diagonal must NOT leak the cm_factor (regression)
        pol = d["polarizability"]
        self.assertEqual(len(pol), 3)
        self.assertAlmostEqual(pol[0][0], 42.102858, places=5)
        self.assertAlmostEqual(pol[1][1], 42.102858, places=5)
        self.assertAlmostEqual(pol[2][2], 42.103279, places=5)
        self.assertTrue(d["has_raman"])
        self.assertEqual(len(d["modes"]), 6)
        m4 = d["modes"][3]
        self.assertEqual(m4["mode"], 4)
        self.assertAlmostEqual(m4["freq"], 332.33, places=2)
        self.assertAlmostEqual(m4["raman"], 354.9075, places=3)
        self.assertAlmostEqual(m4["depol"], 0.75, places=3)

    def test_lif_no_raman(self):
        d = parse_dynmat(_read("LiF_mp-1138", "dynmat.out"))
        self.assertFalse(d["has_raman"])
        self.assertEqual(len(d["modes"]), 6)
        for m in d["modes"]:
            self.assertEqual(m["raman"], 0.0)
            self.assertIsNone(m["depol"])

    def test_asr(self):
        self.assertAlmostEqual(parse_asr(_read("LiF_mp-1138", "dynmat.out")),
                               0.01301689, places=5)


class TestParseDielectric(unittest.TestCase):
    def test_lif(self):
        eps = parse_dielectric(_read("LiF_mp-1138", "ph.out"))
        self.assertIsNotNone(eps)
        self.assertAlmostEqual(eps[0][0], 2.2247, places=4)
        self.assertAlmostEqual(eps[2][2], 2.2247, places=4)

    def test_alas(self):
        eps = parse_dielectric(_read("AlAs_mp-2172", "ph.out"))
        self.assertIsNotNone(eps)
        self.assertAlmostEqual(eps[0][0], 12.2518, places=4)


class TestStats(unittest.TestCase):
    def test_alas(self):
        modes = parse_dynmat(_read("AlAs_mp-2172", "dynmat.out"))["modes"]
        s = _stats(modes)
        self.assertEqual(s["total"], 6)
        self.assertEqual(s["acoustic"], [1, 2, 3])
        self.assertEqual(s["optical"], [4, 5, 6])
        self.assertEqual(s["raman_active"], [4, 5, 6])
        self.assertAlmostEqual(s["sum_raman"], 1064.4515, places=4)
        self.assertEqual(s["soft_modes"], [])

    def test_mos2_soft_mode(self):
        modes = parse_dynmat(_read("MoS2_mp-2815", "dynmat.out"))["modes"]
        s = _stats(modes)
        self.assertEqual(s["total"], 18)
        self.assertTrue(s["soft_modes"])  # -17.29 cm-1 soft optical mode
        self.assertAlmostEqual(s["sum_raman"], 2004.0624, places=4)

    def test_lif_no_raman(self):
        modes = parse_dynmat(_read("LiF_mp-1138", "dynmat.out"))["modes"]
        s = _stats(modes)
        self.assertEqual(s["raman_active"], [])
        self.assertEqual(s["sum_raman"], 0.0)


if __name__ == "__main__":
    unittest.main()
