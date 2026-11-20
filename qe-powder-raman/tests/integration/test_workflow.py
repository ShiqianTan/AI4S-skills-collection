"""Dry-run integration tests: parse validated logs end-to-end (no MPI/QE)."""
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
sys.path.insert(0, str(SCRIPTS))

from qe_parse import parse_outdir, write_outputs  # noqa: E402


class TestParseOutdir(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_int_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _copy(self, system):
        for name in ("scf.out", "ph.out", "dynmat.out"):
            shutil.copy(FIXTURES / system / name, self.tmp / name)
        src = FIXTURES / system
        for lo in src.glob("dynmat_lo*.out"):
            shutil.copy(lo, self.tmp / lo.name)

    def test_alas_success_with_lo(self):
        self._copy("AlAs_mp-2172")
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "success")
        self.assertTrue(rec["ph_job_done"])
        self.assertEqual(rec["stats"]["raman_active"], [4, 5, 6])
        self.assertIn("lo", rec["lo"])
        self.assertAlmostEqual(rec["lo"]["lo"]["modes"][5]["freq"], 365.0, places=2)
        write_outputs(rec)
        self.assertTrue((self.tmp / "raman_modes.csv").exists())
        self.assertTrue((self.tmp / "summary.json").exists())

    def test_lif_success_no_raman(self):
        self._copy("LiF_mp-1138")
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "success")
        self.assertEqual(rec["stats"]["raman_active"], [])
        self.assertEqual(rec["stats"]["sum_raman"], 0.0)

    def test_mos2_warning_soft_mode(self):
        self._copy("MoS2_mp-2815")
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "warning")
        self.assertTrue(rec["stats"]["soft_modes"])
        self.assertAlmostEqual(rec["stats"]["sum_raman"], 2004.0624, places=4)

    def test_empty_dir_failed(self):
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")


if __name__ == "__main__":
    unittest.main()
