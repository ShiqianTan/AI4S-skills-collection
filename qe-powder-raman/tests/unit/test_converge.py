"""Unit tests for parameter convergence scans (E-01)."""
from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from qe_converge import (  # noqa: E402
    CRITERIA,
    DEFAULT_SWEEPS,
    _write_outputs,
    compare_metrics,
    evaluate_fwhm,
    evaluate_group,
    metrics_from_rec,
    render_markdown,
    sweep_values,
    variant_config,
    variant_tag,
)
from qe_report import render_report  # noqa: E402

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def _rec(freqs, ramans, etot=-10.0, eps=5.0):
    modes = [{"mode": i + 1, "freq": f, "raman": a, "ir": 0.0, "thz": f / 33.356}
             for i, (f, a) in enumerate(zip(freqs, ramans))]
    return {
        "status": "success",
        "scf": {"etot_ry": etot},
        "electronic": {"gap_ev": 2.0},
        "dielectric": [[eps, 0, 0], [0, eps, 0], [0, 0, eps]],
        "stats": {
            "optical": [m["mode"] for m in modes[3:]],
            "max_raman": max(ramans[3:]) if len(ramans) > 3 else None,
            "max_raman_freq": freqs[3] if len(freqs) > 3 else None,
            "sum_raman": sum(ramans[3:]),
            "acoustic_max_abs": 0.0,
            "soft_modes": [],
        },
        "modes": modes,
        "warnings": [],
    }


class TestSweeps(unittest.TestCase):
    def test_defaults(self):
        self.assertEqual(sweep_values("ecutwfc"), DEFAULT_SWEEPS["ecutwfc"])
        self.assertEqual(sweep_values("fwhm"), DEFAULT_SWEEPS["fwhm"])

    def test_explicit_sorted_unique(self):
        self.assertEqual(sweep_values("ecutwfc", [70, 50, 50]), [50.0, 70.0])

    def test_invalid(self):
        with self.assertRaises(ValueError):
            sweep_values("tr2_ph", [0.0])
        with self.assertRaises(ValueError):
            sweep_values("ecutwfc", [-1])
        with self.assertRaises(ValueError):
            sweep_values("nope")


class TestVariantConfig(unittest.TestCase):
    def test_ecutwfc_keeps_ratio(self):
        cfg = {"ecutwfc": 60.0, "ecutrho": 240.0}
        out = variant_config(cfg, "ecutwfc", 40)
        self.assertEqual(out["ecutwfc"], 40.0)
        self.assertEqual(out["ecutrho"], 160.0)
        self.assertEqual(cfg["ecutwfc"], 60.0)

    def test_kpt_density_drops_explicit_grid(self):
        out = variant_config({"kpt_density": 30, "kpt_grid": [2, 2, 2]},
                             "kpt_density", 15)
        self.assertEqual(out["kpt_density"], 15.0)
        self.assertNotIn("kpt_grid", out)

    def test_tags(self):
        self.assertEqual(variant_tag("ecutwfc", 60), "ecutwfc_60")
        self.assertEqual(variant_tag("kpt_density", 20), "kpt_20")
        self.assertEqual(variant_tag("tr2_ph", 1e-10), "tr2_1e-10")


class TestCompare(unittest.TestCase):
    def setUp(self):
        self.ref = metrics_from_rec(_rec(
            [0, 0, 0, 300, 400, 500], [0, 0, 0, 10, 20, 30]), n_atoms=2)

    def test_converged(self):
        other = metrics_from_rec(_rec(
            [0, 0, 0, 301.5, 400.5, 500.5], [0, 0, 0, 10.3, 20.6, 30.9]),
            n_atoms=2)
        cmp = compare_metrics(other, self.ref)
        self.assertTrue(cmp["available"])
        self.assertAlmostEqual(cmp["max_freq_dev"], 1.5)
        self.assertLess(cmp["max_activity_rel"], 0.05)
        self.assertTrue(cmp["converged"])

    def test_frequency_fail(self):
        other = metrics_from_rec(_rec(
            [0, 0, 0, 305, 400, 500], [0, 0, 0, 10, 20, 30]), n_atoms=2)
        cmp = compare_metrics(other, self.ref)
        self.assertFalse(cmp["checks"]["freq"])
        self.assertFalse(cmp["converged"])

    def test_activity_fail(self):
        other = metrics_from_rec(_rec(
            [0, 0, 0, 300, 400, 500], [0, 0, 0, 12, 20, 30]), n_atoms=2)
        cmp = compare_metrics(other, self.ref)
        self.assertFalse(cmp["checks"]["activity"])

    def test_missing_variant(self):
        cmp = compare_metrics({"n_modes": 0}, self.ref)
        self.assertFalse(cmp["available"])
        self.assertFalse(cmp["converged"])


class TestEvaluateGroup(unittest.TestCase):
    def test_reference_and_recommendation(self):
        values = [40.0, 50.0, 60.0]
        good = metrics_from_rec(_rec([0, 0, 0, 300, 400, 500],
                                     [0, 0, 0, 10, 20, 30]), n_atoms=2)
        group = evaluate_group("ecutwfc", values, [good, good, good])
        self.assertTrue(group["converged"])
        self.assertEqual(group["reference_value"], 60.0)
        self.assertEqual(group["recommended"], 40.0)
        self.assertTrue(group["variants"][2]["is_reference"])

    def test_tr2_reference_is_lowest(self):
        good = metrics_from_rec(_rec([0, 0, 0, 300, 400, 500],
                                     [0, 0, 0, 10, 20, 30]), n_atoms=2)
        group = evaluate_group("tr2_ph", [1e-8, 1e-12], [good, good])
        self.assertEqual(group["reference_value"], 1e-12)
        self.assertEqual(group["recommended"], 1e-8)

    def test_unconverged(self):
        values = [40.0, 60.0]
        coarse = metrics_from_rec(_rec([0, 0, 0, 320, 400, 500],
                                       [0, 0, 0, 10, 20, 30]), n_atoms=2)
        fine = metrics_from_rec(_rec([0, 0, 0, 300, 400, 500],
                                     [0, 0, 0, 10, 20, 30]), n_atoms=2)
        group = evaluate_group("ecutwfc", values, [coarse, fine])
        self.assertFalse(group["converged"])
        self.assertIsNone(group["recommended"])


class TestFwhm(unittest.TestCase):
    def test_peak_height_and_area(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            rows = [("to", 1, 0.0, 0.0, 0.0, 0.0, ""),
                    ("to", 2, 0.0, 0.0, 0.0, 0.0, ""),
                    ("to", 3, 0.0, 0.0, 0.0, 0.0, ""),
                    ("to", 4, 300.0, 9.0, 0.0, 10.0, "0.75")]
            with (root / "raman_modes.csv").open("w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["table", "mode", "freq_cm1", "freq_thz",
                            "ir", "raman", "depol"])
                w.writerows(rows)
            group = evaluate_fwhm([4.0, 8.0], root)
            a, b = group["variants"]
            self.assertTrue(a["available"])
            self.assertAlmostEqual(a["max_intensity"], b["max_intensity"],
                                   places=4)
            self.assertGreater(b["integral"], a["integral"])


class TestOutputs(unittest.TestCase):
    def test_write_and_render(self):
        good = metrics_from_rec(_rec([0, 0, 0, 300, 400, 500],
                                     [0, 0, 0, 10, 20, 30]), n_atoms=2)
        group = evaluate_group("ecutwfc", [60.0, 70.0], [good, good])
        result = {
            "structure": "x.cif", "outdir": "/tmp/x", "execute": True,
            "criteria": dict(CRITERIA), "groups": [group],
        }
        text = render_markdown(result)
        self.assertIn("收敛扫描报告", text)
        self.assertIn("推荐最小可用值", text)
        with tempfile.TemporaryDirectory() as td:
            paths = _write_outputs(Path(td), result)
            names = {p.name for p in paths}
            self.assertEqual(names, {"convergence.json", "convergence.csv",
                                     "convergence.md"})
            data = json.loads((Path(td) / "convergence.json").read_text())
            self.assertEqual(data["groups"][0]["group"], "ecutwfc")

    def test_report_section(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            (out / "convergence").mkdir()
            (out / "convergence" / "convergence.json").write_text(json.dumps({
                "groups": [{"group": "ecutwfc", "label": "截断能",
                            "reference_value": 70, "recommended": 60,
                            "converged": True,
                            "variants": [{"max_freq_dev": 0.4,
                                          "max_activity_rel": 0.01}]}]
            }))
            text = render_report(out)
            self.assertIn("## 8. 收敛扫描", text)
            self.assertIn("截断能", text)
            self.assertIn("1.00%", text)


class TestConvergeCli(unittest.TestCase):
    def test_dry_run_prepares_variants(self):
        sys.path.insert(0, str(SCRIPTS))
        import qe_raman

        with tempfile.TemporaryDirectory() as td:
            root = Path(td) / "scan"
            rc = qe_raman.main([
                "converge", str(EXAMPLES / "LiF_mp-1138.cif"),
                "--outdir", str(root), "--groups", "ecutwfc",
                "--ecutwfc", "40", "60", "--np", "1",
            ])
            self.assertEqual(rc, 0)
            for tag, ecut, ecutrho in (("ecutwfc_40", 40, 160),
                                       ("ecutwfc_60", 60, 240)):
                scf = (root / "ecutwfc" / tag / "scf.in").read_text()
                self.assertIn("ecutwfc      = %d" % ecut, scf)
                self.assertIn("ecutrho      = %d" % ecutrho, scf)
                self.assertTrue((root / "ecutwfc" / tag
                                 / "pseudopotential_manifest.json").exists())
            self.assertTrue((root / "convergence.md").exists())
            self.assertTrue((root / "convergence.json").exists())


if __name__ == "__main__":
    unittest.main()
