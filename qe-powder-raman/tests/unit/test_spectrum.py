"""Unit tests for qe_spectrum (physics helpers + spectrum construction)."""
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import numpy as np  # noqa: E402

from qe_spectrum import (  # noqa: E402
    bose_plus_one,
    build_spectrum,
    gauss,
    generate_spectra,
    lo_direction_labels,
    lorentz,
    select_ir_peaks,
    select_peaks,
    stokes_prefactor,
    validate_spectrum_params,
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


class TestParameterValidation(unittest.TestCase):
    def test_fwhm_zero_raises(self):
        with self.assertRaises(ValueError):
            build_spectrum([(300.0, 1.0)], fwhm=0.0)

    def test_bad_params(self):
        with self.assertRaises(ValueError):
            validate_spectrum_params(fwhm=8.0, temperature=-1.0,
                                     laser_wavelength=532.0, threshold=1e-4,
                                     npoints=4000, formats=("png",))
        with self.assertRaises(ValueError):
            validate_spectrum_params(fwhm=8.0, temperature=300.0,
                                     laser_wavelength=0.0, threshold=1e-4,
                                     npoints=4000, formats=("png",))
        with self.assertRaises(ValueError):
            validate_spectrum_params(fwhm=8.0, temperature=300.0,
                                     laser_wavelength=532.0, threshold=1e-4,
                                     npoints=4000, formats=("pdf",))

    def test_ok(self):
        validate_spectrum_params(fwhm=8.0, temperature=300.0,
                                 laser_wavelength=532.0, threshold=1e-4,
                                 npoints=4000, formats=("png", "csv"))


class TestSelectPeaksClassification(unittest.TestCase):
    def test_low_frequency_optical_kept(self):
        modes = [
            {"mode": 1, "freq": 0.0, "raman": 0.0},
            {"mode": 2, "freq": 0.0, "raman": 0.0},
            {"mode": 3, "freq": 0.0, "raman": 0.0},
            {"mode": 4, "freq": 2.0, "raman": 5.0},   # low-freq optical
            {"mode": 5, "freq": 100.0, "raman": 6.0},
        ]
        peaks = select_peaks(modes, threshold=1e-4, min_freq=1.0)
        self.assertEqual([f for f, _ in peaks], [2.0, 100.0])


class TestSelectIrPeaks(unittest.TestCase):
    def test_filter(self):
        modes = [
            {"mode": 1, "freq": 0.0, "ir": 0.0},
            {"mode": 2, "freq": 0.0, "ir": 0.0},
            {"mode": 3, "freq": 0.0, "ir": 0.0},
            {"mode": 4, "freq": -20.0, "ir": 9.0},   # soft (imaginary)
            {"mode": 5, "freq": 100.0, "ir": 5.0},   # IR active
            {"mode": 6, "freq": 200.0, "ir": 0.0},   # silent
        ]
        self.assertEqual(select_ir_peaks(modes), [(100.0, 5.0)])


class TestLoDirectionLabels(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_spec_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_multiple_directions(self):
        (self.tmp / "inputs_manifest.json").write_text(
            json.dumps({"lo_dirs": [[0, 0, 1], [1, 0, 0]]}))
        labels = lo_direction_labels(self.tmp, ["lo_1", "lo_2"])
        self.assertEqual(labels["lo_1"], "q=(0.000, 0.000, 1.000)")
        self.assertEqual(labels["lo_2"], "q=(1.000, 0.000, 0.000)")

    def test_single_direction(self):
        (self.tmp / "inputs_manifest.json").write_text(
            json.dumps({"lo_dirs": [[1, 0, 0]]}))
        self.assertEqual(lo_direction_labels(self.tmp, ["lo"])["lo"],
                         "q=(1.000, 0.000, 0.000)")

    def test_missing_manifest(self):
        self.assertEqual(lo_direction_labels(self.tmp, ["lo_1"]), {})


class TestGenerateSpectra(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_spec_"))
        rows = ["table,mode,freq_cm1,freq_thz,ir,raman,depol"]
        for table, shift in (("to", 0.0), ("lo_1", 30.0)):
            for mode in range(1, 4):
                rows.append("%s,%d,0.000,0.000,0.000,0.000,0.750" % (table, mode))
            rows += [
                "%s,4,%.3f,%.3f,5.000,10.000,0.750" % (table, 100.0 + shift, 3.0),
                "%s,5,%.3f,%.3f,0.000,4.000,0.750" % (table, 200.0 + shift, 6.0),
                "%s,6,%.3f,%.3f,2.000,0.000,0.000" % (table, 300.0 + shift, 9.0),
            ]
        (self.tmp / "raman_modes.csv").write_text("\n".join(rows) + "\n")
        (self.tmp / "inputs_manifest.json").write_text(
            json.dumps({"lo_dirs": [[1, 0, 0]]}))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _names(self, written):
        return {p.name for p in written}

    def test_ir_and_lo_outputs(self):
        names = self._names(generate_spectra(self.tmp, formats=("png", "csv")))
        for expected in ("spectrum_activity.csv", "spectrum_activity.png",
                         "spectrum_ir.csv", "spectrum_ir.png",
                         "spectrum_raman_ir.png",
                         "spectrum_lo_1.csv", "spectrum_lo_1.png",
                         "spectrum_to_vs_lo_1.png"):
            self.assertIn(expected, names)

    def test_ir_spectrum_has_peaks(self):
        generate_spectra(self.tmp, formats=("csv",))
        lines = (self.tmp / "spectrum_ir.csv").read_text().splitlines()
        self.assertGreater(len(lines), 2)

    def test_disable_ir_and_lo(self):
        names = self._names(generate_spectra(
            self.tmp, formats=("csv",), include_ir=False, include_lo=False))
        self.assertIn("spectrum_activity.csv", names)
        self.assertNotIn("spectrum_ir.csv", names)
        self.assertNotIn("spectrum_lo_1.csv", names)

    def test_raman_blank_but_ir_present(self):
        rows = ["table,mode,freq_cm1,freq_thz,ir,raman,depol"]
        for mode in range(1, 4):
            rows.append("to,%d,0.000,0.000,0.000,0.000,0.750" % mode)
        rows += ["to,4,280.000,8.400,5.000,0.000,0.750",
                 "to,5,280.500,8.415,5.000,0.000,0.750",
                 "to,6,280.500,8.415,5.000,0.000,0.750"]
        (self.tmp / "raman_modes.csv").write_text("\n".join(rows) + "\n")
        generate_spectra(self.tmp, formats=("csv",), include_lo=False)
        self.assertEqual(len((self.tmp / "spectrum_activity.csv")
                             .read_text().splitlines()), 1)
        self.assertGreater(len((self.tmp / "spectrum_ir.csv")
                               .read_text().splitlines()), 2)


if __name__ == "__main__":
    unittest.main()
