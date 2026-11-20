"""Unit tests for qe_spectrum (physics helpers + spectrum construction)."""
import sys
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import numpy as np  # noqa: E402

from qe_spectrum import (  # noqa: E402
    bose_plus_one,
    build_spectrum,
    gauss,
    lorentz,
    select_peaks,
    stokes_prefactor,
)


class TestBose(unittest.TestCase):
    def test_thermal(self):
        n = bose_plus_one(300.0, 300.0)
        self.assertGreater(n, 1.0)
        self.assertLess(n, 4.0)

    def test_limits(self):
        self.assertEqual(bose_plus_one(0.0, 300.0), 1.0)
        self.assertEqual(bose_plus_one(-5.0, 300.0), 1.0)
        self.assertEqual(bose_plus_one(300.0, 0.0), 1.0)
        # x > 50 -> return 1.0
        self.assertEqual(bose_plus_one(20000.0, 300.0), 1.0)

    def test_monotonic(self):
        self.assertGreater(bose_plus_one(100.0, 300.0),
                           bose_plus_one(500.0, 300.0))


class TestStokesPrefactor(unittest.TestCase):
    def test_positive(self):
        self.assertGreater(stokes_prefactor(1000.0, 532.0), 0.0)

    def test_finite_zero_freq(self):
        v = stokes_prefactor(0.0, 532.0)
        self.assertTrue(np.isfinite(v))
        self.assertGreater(v, 0.0)


class TestBroadening(unittest.TestCase):
    def test_lorentz_peak_equals_amp(self):
        x = np.array([100.0])
        self.assertAlmostEqual(float(lorentz(x, 100.0, 5.0, 8.0)[0]), 5.0, places=6)

    def test_gauss_peak_equals_amp(self):
        x = np.array([100.0])
        self.assertAlmostEqual(float(gauss(x, 100.0, 5.0, 8.0)[0]), 5.0, places=6)


class TestSelectPeaks(unittest.TestCase):
    def test_filter(self):
        modes = [
            {"freq": 0.0, "raman": 0.0},       # acoustic
            {"freq": -17.0, "raman": 2.0},     # soft (imaginary)
            {"freq": 100.0, "raman": 0.0},     # silent
            {"freq": 200.0, "raman": 3.0},     # active
            {"freq": 300.0, "raman": 1e-6},    # below threshold
        ]
        peaks = select_peaks(modes, threshold=1e-4, min_freq=1.0)
        self.assertEqual(peaks, [(200.0, 3.0)])


class TestBuildSpectrum(unittest.TestCase):
    def test_activity(self):
        x, y, scaled = build_spectrum([(100.0, 5.0), (200.0, 3.0)],
                                      spectrum_type="activity")
        self.assertEqual(len(x), 4000)
        self.assertEqual(len(y), 4000)
        self.assertEqual(scaled, [(100.0, 5.0), (200.0, 3.0)])
        self.assertGreater(float(y.max()), 0.0)

    def test_stokes_larger_than_activity(self):
        _, y_act, _ = build_spectrum([(300.0, 5.0)], spectrum_type="activity")
        _, y_sto, _ = build_spectrum([(300.0, 5.0)], spectrum_type="stokes",
                                     temperature=300.0, laser_wavelength=532.0)
        self.assertGreater(float(y_sto.max()), float(y_act.max()))

    def test_no_peaks(self):
        x, y, scaled = build_spectrum([], spectrum_type="activity")
        self.assertIsNone(x)
        self.assertIsNone(y)
        self.assertEqual(scaled, [])


if __name__ == "__main__":
    unittest.main()
