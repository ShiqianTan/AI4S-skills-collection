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
from qe_report import render_report  # noqa: E402
from qe_spectrum import generate_spectra  # noqa: E402


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

    def test_missing_requested_lo_fails(self):
        self._copy("AlAs_mp-2172")
        (self.tmp / "inputs_manifest.json").write_text('{"lo_dirs": [[0, 0, 1]]}')
        (self.tmp / "dynmat_lo.out").unlink()
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")
        self.assertIn("dynmat_lo.out missing", rec["status_reason"])

    def test_truncated_lo_fails_even_with_completion_marker(self):
        self._copy("AlAs_mp-2172")
        log = self.tmp / "dynmat_lo.out"
        text = log.read_text().splitlines()
        header = next(i for i, line in enumerate(text) if line.startswith("# mode"))
        log.write_text("\n".join(text[:header + 2]) + "\nJOB DONE.\n")
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")
        self.assertIn("expected 6 consecutive LO modes; found 1", rec["status_reason"])

    def test_unrequested_lo_log_is_ignored(self):
        self._copy("AlAs_mp-2172")
        (self.tmp / "inputs_manifest.json").write_text('{"lo_dirs": []}')
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "success")
        self.assertEqual(rec["lo"], {})

    def test_lif_success_no_raman(self):
        self._copy("LiF_mp-1138")
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "success")
        self.assertEqual(rec["stats"]["raman_active"], [])
        self.assertEqual(rec["stats"]["sum_raman"], 0.0)

    def test_alas_ir_and_per_lo_spectra(self):
        self._copy("AlAs_mp-2172")
        (self.tmp / "inputs_manifest.json").write_text('{"lo_dirs": [[1, 0, 0]]}')
        write_outputs(parse_outdir(self.tmp))
        generate_spectra(self.tmp, formats=("png", "csv"))
        for name in ("spectrum_ir.png", "spectrum_raman_ir.png",
                     "spectrum_lo.png", "spectrum_to_vs_lo.png"):
            self.assertTrue((self.tmp / name).is_file(), name)
        self.assertGreater(len((self.tmp / "spectrum_lo.csv")
                               .read_text().splitlines()), 2)

    def test_lif_ir_spectrum_present_without_raman(self):
        self._copy("LiF_mp-1138")
        write_outputs(parse_outdir(self.tmp))
        written = generate_spectra(self.tmp, formats=("csv",), include_lo=False)
        names = {p.name for p in written}
        self.assertIn("spectrum_ir.csv", names)
        self.assertEqual(len((self.tmp / "spectrum_activity.csv")
                             .read_text().splitlines()), 1)
        self.assertGreater(len((self.tmp / "spectrum_ir.csv")
                               .read_text().splitlines()), 2)

    def test_mos2_warning_soft_mode(self):
        self._copy("MoS2_mp-2815")
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "warning")
        self.assertTrue(rec["stats"]["soft_modes"])
        self.assertAlmostEqual(rec["stats"]["sum_raman"], 2004.0624, places=4)

    def test_empty_dir_failed(self):
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")

    def test_delete_scf_failed(self):
        self._copy("LiF_mp-1138")
        (self.tmp / "scf.out").unlink()
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")

    def test_stale_dynmat_state_rejects_complete_old_log(self):
        self._copy("LiF_mp-1138")
        (self.tmp / "state.json").write_text('{"stages": {"dynmat": {"status": "stale"}}}')
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")
        self.assertIn("dynmat state is stale", rec["status_reason"])

    def test_missing_job_done_dynmat_failed(self):
        self._copy("AlAs_mp-2172")
        txt = (self.tmp / "dynmat.out").read_text()
        (self.tmp / "dynmat.out").write_text(txt.replace("JOB DONE.", ""))
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")

    def test_truncated_modes_manifest_mismatch_failed(self):
        self._copy("AlAs_mp-2172")
        header = "# mode   [cm-1]    [THz]      IR          Raman   depol.fact\n"
        rows = (
            "    1    332.33   9.9631   6.1453  354.9075   0.7500\n"
            "    2    332.38   9.9645   6.1453  354.7720   0.7500\n"
        )
        (self.tmp / "dynmat.out").write_text(header + rows + "\nJOB DONE.\n")
        (self.tmp / "pseudopotential_manifest.json").write_text('{"n_atoms": 2}')
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")

    def test_manifest_n_atoms_mismatch_failed(self):
        self._copy("LiF_mp-1138")
        (self.tmp / "pseudopotential_manifest.json").write_text('{"n_atoms": 3}')
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")

    def test_metallic_failed(self):
        self._copy("LiF_mp-1138")
        txt = (self.tmp / "scf.out").read_text()
        (self.tmp / "scf.out").write_text(
            txt + "\n     lowest unoccupied level (ev):    -1.0000\n")
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["electronic"]["classification"], "metallic")
        self.assertEqual(rec["status"], "failed")

    def test_undeterminable_is_warning_not_failure(self):
        self._copy("AlAs_mp-2172")
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["electronic"]["classification"], "undeterminable")
        self.assertEqual(rec["status"], "success")
        self.assertTrue(any("undeterminable" in w for w in rec["warnings"]))

    def test_lraman_tensor_missing_failed(self):
        self._copy("AlAs_mp-2172")
        ph = self.tmp / "ph.out"
        text = "\n".join(
            line for line in ph.read_text().splitlines()
            if "raman tensor" not in line.lower()
        )
        ph.write_text(text + "\n")
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["status"], "failed")
        self.assertIn("Raman tensor", rec["status_reason"])

    def test_relax_manifest_in_summary_and_report(self):
        self._copy("AlAs_mp-2172")
        (self.tmp / "relax_manifest.json").write_text(
            '{"final_force": 0.0123, "final_stress": 0.45,'
            ' "n_atoms": 2, "spacegroup": "F-43m"}')
        rec = parse_outdir(self.tmp)
        self.assertEqual(rec["relax"]["n_atoms"], 2)
        self.assertAlmostEqual(rec["relax"]["final_force"], 0.0123)
        write_outputs(rec)
        summary = (self.tmp / "summary.json").read_text()
        self.assertIn('"relax"', summary)
        text = render_report(self.tmp)
        self.assertIn("结构弛豫", text)
        self.assertIn("0.0123", text)

    def test_lraman_off_parse(self):
        self._copy("LiF_mp-1138")
        (self.tmp / "ph.in").write_text("&inputph\n  lraman    = .false.\n/\n")
        rec = parse_outdir(self.tmp)
        self.assertFalse(rec["lraman"])
        self.assertFalse(rec["lraman_tensor"])
        self.assertEqual(rec["status"], "success")
        self.assertTrue(any("lraman disabled" in w for w in rec["warnings"]))

    def test_report_records_pseudo_xc_and_source(self):
        self._copy("AlAs_mp-2172")
        (self.tmp / "ph.in").write_text("&inputph\n  lraman    = .false.\n/\n")
        (self.tmp / "effective_config.json").write_text('{"raman": false}')
        (self.tmp / "pseudopotential_manifest.json").write_text(
            '{"n_atoms": 2, "xc_families": ["GGA"], "raman_supported": false,'
            ' "pseudo_source": ["PseudoDojo ONCVPSP-PBE-FR-PDv0.4"],'
            ' "pseudopotentials": {"Al": {"file": "Al.pbe-fr-oncvpsp.UPF",'
            ' "functional": "PBE", "sha256": "%s"}}}' % ("b" * 64))
        rec = parse_outdir(self.tmp)
        write_outputs(rec)
        text = render_report(self.tmp)
        self.assertIn("赝势 XC 家族", text)
        self.assertIn("GGA", text)
        self.assertIn("PseudoDojo", text)
        self.assertIn("仅声子/IR", text)

    def test_report_sections(self):
        self._copy("AlAs_mp-2172")
        (self.tmp / "inputs_manifest.json").write_text(
            '{"lo_dirs": [[1.0, 0.0, 0.0]]}')
        (self.tmp / "effective_config.json").write_text('{"fwhm": 8.0}')
        (self.tmp / "structure_manifest.json").write_text(
            '{"original": {"spacegroup_number": 216, "spacegroup_symbol": "F-43m"},'
            ' "refined": {"spacegroup_number": 216, "spacegroup_symbol": "F-43m"},'
            ' "standardized": false, "cellpar": [4.05, 4.05, 4.05, 60, 60, 60]}')
        (self.tmp / "pseudopotential_manifest.json").write_text(
            '{"n_atoms": 2, "pseudopotentials": {"Al": {"file": "Al.pz-hgh.UPF",'
            ' "sha256": "%s"}}}' % ("a" * 64))
        rec = parse_outdir(self.tmp)
        write_outputs(rec)
        text = render_report(self.tmp)
        self.assertIn("电子结构适用性", text)
        self.assertIn("Raman 张量生成", text)
        self.assertIn("## 7. 可复现信息", text)
        self.assertIn("Raman CM 因子", text)
        self.assertIn("q=(1.000", text)


class TestAssignmentAndLayers(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_int_e04_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_assign_and_report_section(self):
        from qe_assign import (assign_modes, match_experimental,
                               write_assignment)
        from qe_spectrum import load_modes

        for name in ("scf.out", "ph.out", "dynmat.out", "dynmat_lo.out"):
            shutil.copy(FIXTURES / "AlAs_mp-2172" / name, self.tmp / name)
        rec = parse_outdir(self.tmp)
        write_outputs(rec)
        assignment = assign_modes(load_modes(self.tmp)["to"])
        matching = match_experimental(
            assignment["entries"],
            [{"freq": 331.0, "intensity": 1.0, "error": 0.1},
             {"freq": 364.0, "intensity": 0.5, "error": 0.1}])
        write_assignment(self.tmp, assignment, matching)
        self.assertTrue((self.tmp / "mode_assignment.csv").exists())
        self.assertTrue((self.tmp / "assignment_overlay.png").exists())
        text = render_report(self.tmp)
        self.assertIn("## 9. 模式指认", text)
        self.assertIn("Raman+IR", text)

    def test_report_layered_warning(self):
        for name in ("scf.out", "ph.out", "dynmat.out"):
            shutil.copy(FIXTURES / "MoS2_mp-2815" / name, self.tmp / name)
        write_outputs(parse_outdir(self.tmp))
        (self.tmp / "structure_manifest.json").write_text(
            '{"layers": {"is_layered": true, "is_slab": false, "axis": "c",'
            ' "vacuum_thickness": 4.31, "layer_count": 2,'
            ' "gap_elements": ["S", "S"]}}')
        text = render_report(self.tmp)
        self.assertIn("层状/低维提示", text)
        self.assertIn("层状/范德华间隙", text)


if __name__ == "__main__":
    unittest.main()
