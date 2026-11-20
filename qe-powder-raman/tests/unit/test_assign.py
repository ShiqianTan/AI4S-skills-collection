"""Unit tests for qe_assign (E-04 mode assignment + experimental comparison)."""
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from qe_assign import (  # noqa: E402
    assign_modes,
    detect_peaks,
    fit_single_peak,
    match_experimental,
    read_experimental,
    write_assignment,
    _lorentz,
)


def _mode(mode, freq, *, ir=0.0, raman=0.0, depol=0.75):
    return {"mode": mode, "freq": freq, "thz": freq / 33.35641,
            "ir": ir, "raman": raman, "depol": depol}


MODES = [
    _mode(1, -0.01),
    _mode(2, 0.0),
    _mode(3, 0.01),
    _mode(4, -20.0),
    _mode(5, 100.0, ir=5.0),
    _mode(6, 200.0, raman=3.0),
    _mode(7, 300.0, ir=2.0, raman=4.0),
    _mode(8, 300.4, ir=2.0, raman=4.0),
    _mode(9, 400.0),
]


class TestAssignModes(unittest.TestCase):
    def test_activity_labels_and_degeneracy(self):
        result = assign_modes(MODES, threshold=1e-4, degeneracy_tol=1.0)
        by_mode = {e["mode"]: e for e in result["entries"]}
        self.assertEqual(len(result["acoustic"]), 3)
        self.assertEqual(by_mode[4]["activity"], "soft")
        self.assertEqual(by_mode[5]["activity"], "ir")
        self.assertEqual(by_mode[6]["activity"], "raman")
        self.assertEqual(by_mode[7]["activity"], "both")
        self.assertEqual(by_mode[9]["activity"], "silent")
        self.assertEqual(by_mode[7]["degeneracy"], 2)
        self.assertEqual(by_mode[7]["partners"], [7, 8])
        self.assertEqual(by_mode[7]["group"], by_mode[8]["group"])

    def test_threshold_changes_labels(self):
        result = assign_modes(MODES, threshold=4.5)
        by_mode = {e["mode"]: e for e in result["entries"]}
        self.assertEqual(by_mode[6]["activity"], "silent")
        self.assertEqual(by_mode[5]["activity"], "ir")


class TestExperimental(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_assign_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_read_peak_table(self):
        path = self.tmp / "exp.csv"
        path.write_text("freq_cm1,intensity,error\n100.0,2.0,0.1\n200,3,0.2\n")
        rows = read_experimental(path)
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(rows[0]["freq"], 100.0)
        self.assertAlmostEqual(rows[1]["error"], 0.2)

    def test_match_with_tolerance(self):
        result = assign_modes(MODES)
        matches = match_experimental(
            result["entries"],
            [{"freq": 101.0, "intensity": 1.0, "error": None},
             {"freq": 500.0, "intensity": 1.0, "error": None}],
            tol_cm1=5.0,
        )
        self.assertEqual(matches["n_matched"], 1)
        self.assertEqual(matches["matches"][0]["predicted_mode"], 5)
        self.assertAlmostEqual(matches["matches"][0]["delta"], 1.0)
        self.assertFalse(matches["matches"][1]["matched"])
        self.assertAlmostEqual(matches["rms_delta"], 1.0)


class TestFitting(unittest.TestCase):
    def setUp(self):
        self.x = np.linspace(250.0, 350.0, 1001)
        self.y = _lorentz(self.x, 300.0, 5.0, 10.0)

    def test_detect_peaks(self):
        peaks = detect_peaks(self.x, self.y)
        self.assertEqual(len(peaks), 1)
        self.assertAlmostEqual(peaks[0], 300.0, delta=0.5)

    def test_fit_single_peak(self):
        fit = fit_single_peak(self.x, self.y, 300.0)
        self.assertIsNotNone(fit)
        self.assertAlmostEqual(fit["center"], 300.0, delta=1.0)
        self.assertAlmostEqual(fit["fwhm"], 10.0, delta=2.0)
        self.assertGreater(fit["amplitude"], 4.0)

    def test_write_outputs(self):
        tmp = Path(tempfile.mkdtemp(prefix="qer_assign_w_"))
        self.addCleanup(shutil.rmtree, tmp, True)
        result = assign_modes(MODES)
        matching = match_experimental(
            result["entries"], [{"freq": 300.2, "intensity": 1.0, "error": 0.1}])
        written = write_assignment(tmp, result, matching, None, None,
                                   formats=("png", "csv"))
        names = {p.name for p in written}
        self.assertIn("mode_assignment.csv", names)
        self.assertIn("assignment.md", names)
        self.assertIn("assignment_overlay.png", names)
        text = (tmp / "assignment.md").read_text()
        self.assertIn("模式指认", text)
        self.assertIn("实验峰匹配", text)


if __name__ == "__main__":
    unittest.main()
