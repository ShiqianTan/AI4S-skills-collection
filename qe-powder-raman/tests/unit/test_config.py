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


if __name__ == "__main__":
    unittest.main()
