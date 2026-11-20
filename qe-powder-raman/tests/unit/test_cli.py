"""Unit tests for CLI helpers (version parsing, k-point overrides)."""
import argparse
import os
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from qe_common import QE_MIN_VERSION  # noqa: E402
from qe_raman import _apply_cli_overrides, _check_binaries, _parse_version  # noqa: E402


class TestVersionParsing(unittest.TestCase):
    def test_develop_version(self):
        self.assertEqual(_parse_version("Program PWSCF v.8.0dev starts"), (8, 0))

    def test_release_version(self):
        self.assertEqual(_parse_version("Program PWSCF v.7.2"), (7, 2))

    def test_unparseable(self):
        self.assertIsNone(_parse_version("some other output"))

    def test_min_version_constant(self):
        self.assertEqual(QE_MIN_VERSION, (7, 0))


class TestCliOverrides(unittest.TestCase):
    def _args(self, **kwargs):
        defaults = {"kpt_density": None, "kpt_grid": None}
        defaults.update(kwargs)
        return argparse.Namespace(**defaults)

    def test_density_override(self):
        cfg = _apply_cli_overrides({"kpt_density": 30.0}, self._args(kpt_density=6.0))
        self.assertEqual(cfg["kpt_density"], 6.0)

    def test_grid_override(self):
        cfg = _apply_cli_overrides(
            {"kpt_density": 30.0}, self._args(kpt_grid=(2, 3, 4))
        )
        self.assertEqual(cfg["kpt_grid"], [2, 3, 4])

    def test_no_override(self):
        cfg = _apply_cli_overrides({"kpt_density": 30.0}, self._args())
        self.assertEqual(cfg["kpt_density"], 30.0)
        self.assertNotIn("kpt_grid", cfg)


class TestPreflightScratchDir(unittest.TestCase):
    def test_version_probe_does_not_pollute_cwd(self):
        tmp = Path(tempfile.mkdtemp(prefix="qer_cli_"))
        bin_dir = tmp / "bin"
        bin_dir.mkdir()
        for name in ("ph.x", "dynmat.x"):
            exe = bin_dir / name
            exe.write_text("#!/bin/sh\nexit 0\n")
            exe.chmod(0o755)
        pw = bin_dir / "pw.x"
        pw.write_text(
            "#!/bin/sh\n"
            "echo 'Program PWSCF v.8.0dev starts'\n"
            "printf 'crash' > CRASH\n"
            "printf '' > input_tmp.in\n"
        )
        pw.chmod(0o755)
        work = tmp / "work"
        work.mkdir()
        old = Path.cwd()
        os.chdir(work)
        try:
            msgs = _check_binaries(bin_dir)
        finally:
            os.chdir(old)
        self.assertFalse((work / "CRASH").exists())
        self.assertFalse((work / "input_tmp.in").exists())
        self.assertTrue(any(k == "qe_version" and lvl == "PASS" for lvl, k, _ in msgs))


if __name__ == "__main__":
    unittest.main()
