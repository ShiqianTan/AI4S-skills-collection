"""Unit tests for qe_common.validate_config + config defaults (P0-03)."""
import sys
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from qe_common import load_defaults, validate_config  # noqa: E402


class TestValidateConfig(unittest.TestCase):
    def test_defaults_ok(self):
        msgs = validate_config(load_defaults())
        self.assertFalse(any(lvl == "FAIL" for lvl, _, _ in msgs))

    def test_occupations_smearing_fail(self):
        cfg = load_defaults()
        cfg["occupations"] = "smearing"
        msgs = validate_config(cfg)
        self.assertTrue(
            any(lvl == "FAIL" and key == "occupations" for lvl, key, _ in msgs)
        )

    def test_nan_fwhm_fail(self):
        cfg = load_defaults()
        cfg["fwhm"] = float("nan")
        msgs = validate_config(cfg)
        self.assertTrue(any(lvl == "FAIL" and key == "fwhm" for lvl, key, _ in msgs))

    def test_ph_symmetry_switch(self):
        from qe_inputs import generate_ph

        cfg = load_defaults()
        self.assertIn("search_sym = .true.", generate_ph(cfg, ["Li"], "test"))
        cfg["ph_search_sym"] = False
        self.assertIn("search_sym = .false.", generate_ph(cfg, ["Li"], "test"))
        cfg["ph_search_sym"] = "false"
        msgs = validate_config(cfg)
        self.assertTrue(any(lvl == "FAIL" and key == "ph_search_sym" for lvl, key, _ in msgs))

    def test_zero_ecutwfc_fail(self):
        cfg = load_defaults()
        cfg["ecutwfc"] = 0.0
        msgs = validate_config(cfg)
        self.assertTrue(any(lvl == "FAIL" and key == "ecutwfc" for lvl, key, _ in msgs))

    def test_zero_np_fail(self):
        cfg = load_defaults()
        cfg["np"] = 0
        msgs = validate_config(cfg)
        self.assertTrue(any(lvl == "FAIL" and key == "np" for lvl, key, _ in msgs))


class TestDefaults(unittest.TestCase):
    def test_outdir_qe_outdir_split(self):
        cfg = load_defaults()
        self.assertEqual(cfg["outdir"], "runs")
        self.assertEqual(cfg["qe_outdir"], "./tmp")


class TestPhaseCValidation(unittest.TestCase):
    def test_low_ecutrho_warns(self):
        cfg = load_defaults()
        cfg["ecutwfc"] = 60.0
        cfg["ecutrho"] = 120.0
        msgs = validate_config(cfg)
        self.assertTrue(
            any(lvl == "WARN" and key == "ecutrho" for lvl, key, _ in msgs)
        )
        self.assertFalse(any(lvl == "FAIL" and key == "ecutrho" for lvl, key, _ in msgs))

    def test_ecutrho_ok(self):
        cfg = load_defaults()
        cfg["ecutwfc"] = 60.0
        cfg["ecutrho"] = 240.0
        msgs = validate_config(cfg)
        self.assertFalse(any(key == "ecutrho" for _, key, _ in msgs))

    def test_kpt_grid_bad_fail(self):
        for bad in ([0, 2, 2], [1, 2], [-1, 2, 2]):
            cfg = load_defaults()
            cfg["kpt_grid"] = bad
            msgs = validate_config(cfg)
            self.assertTrue(
                any(lvl == "FAIL" and key == "kpt_grid" for lvl, key, _ in msgs),
                "expected kpt_grid FAIL for %r" % (bad,),
            )

    def test_kpt_grid_ok(self):
        cfg = load_defaults()
        cfg["kpt_grid"] = [2, 3, 4]
        msgs = validate_config(cfg)
        self.assertFalse(any(key == "kpt_grid" for _, key, _ in msgs))


if __name__ == "__main__":
    unittest.main()
