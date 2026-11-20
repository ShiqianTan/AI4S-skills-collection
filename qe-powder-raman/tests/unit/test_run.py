"""Tests for qe_run.run_workflow (env pass-through, skip/resume) + LO cleanup."""
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
PSEUDO = Path(__file__).resolve().parents[2] / "pseudo"
sys.path.insert(0, str(SCRIPTS))

from qe_run import _artifacts, run_workflow  # noqa: E402


class TestRunWorkflowEnvAndSkip(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_run_"))
        self.outdir = self.tmp / "run"
        self.outdir.mkdir()
        (self.outdir / "scf.in").write_text("&CONTROL\n/\n")
        self.bin_dir = self.tmp / "bin"
        self.bin_dir.mkdir()
        self.counter = self.tmp / "counter.txt"
        pw = self.bin_dir / "pw.x"
        pw.write_text(
            "#!/bin/sh\n"
            "echo run >> {counter}\n"
            "echo OMP_NUM_THREADS=$OMP_NUM_THREADS\n"
            "echo JOB DONE.\n"
            "exit 0\n".format(counter=self.counter)
        )
        pw.chmod(0o755)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_env_passthrough(self):
        state, ok = run_workflow(
            self.outdir, self.bin_dir, PSEUDO, np=1, omp_num_threads=1, stage="scf"
        )
        self.assertTrue(ok)
        self.assertTrue(state)
        log = (self.outdir / "scf.out").read_text()
        self.assertIn("OMP_NUM_THREADS=1", log)

    def test_skip_existing_stable(self):
        for _ in range(3):
            run_workflow(
                self.outdir, self.bin_dir, PSEUDO, np=1,
                stage="scf", skip_existing=True,
            )
        lines = self.counter.read_text().strip().splitlines()
        self.assertEqual(len(lines), 1)
        st = json.loads((self.outdir / "state.json").read_text())
        self.assertEqual(st["stages"]["scf"]["status"], "skipped")

    def test_changed_scf_input_reruns_and_invalidates_downstream(self):
        def run(stage):
            return run_workflow(self.outdir, self.bin_dir, PSEUDO, np=1,
                                stage=stage, skip_existing=True)

        (self.outdir / "ph.in").touch()
        ph = self.bin_dir / "ph.x"
        ph.write_text("#!/bin/sh\necho JOB DONE.\n")
        ph.chmod(0o755)
        self.assertTrue(run("scf")[1])
        self.assertTrue(run("ph")[1])

        (self.outdir / "scf.in").write_text("&CONTROL\n prefix='changed'\n/\n")
        self.assertFalse(run("ph")[1])
        state, ok = run("scf")
        self.assertTrue(ok)
        self.assertEqual(state["stages"]["ph"]["status"], "stale")
        self.assertEqual(len(self.counter.read_text().splitlines()), 2)
        self.assertTrue(run("ph")[1])

    def test_interrupt_persists_completed_scf_and_running_ph(self):
        (self.outdir / "ph.in").touch()
        (self.outdir / "dynmat.in").touch()
        from subprocess import run as real_run

        def interrupt_ph(cmd, **kwargs):
            if cmd[-1] == "ph.in":
                raise KeyboardInterrupt
            return real_run(cmd, **kwargs)

        with patch("qe_run.subprocess.run", side_effect=interrupt_ph):
            with self.assertRaises(KeyboardInterrupt):
                run_workflow(self.outdir, self.bin_dir, PSEUDO, np=1)
        state = json.loads((self.outdir / "state.json").read_text())
        self.assertEqual(state["stages"]["scf"]["status"], "success")
        self.assertEqual(state["stages"]["ph"]["status"], "running")
        for name in ("ph.x", "dynmat.x"):
            exe = self.bin_dir / name
            exe.write_text("#!/bin/sh\necho JOB DONE.\n")
            exe.chmod(0o755)
        state, ok = run_workflow(self.outdir, self.bin_dir, PSEUDO, np=1, skip_existing=True)
        self.assertTrue(ok)
        self.assertEqual(state["stages"]["scf"]["status"], "skipped")
        self.assertEqual(state["stages"]["ph"]["status"], "success")

    def test_zero_exit_with_incomplete_log_fails(self):
        (self.bin_dir / "pw.x").write_text("#!/bin/sh\necho interrupted\n")
        state, ok = run_workflow(self.outdir, self.bin_dir, PSEUDO, np=1, stage="scf")
        self.assertFalse(ok)
        self.assertEqual(state["stages"]["scf"]["status"], "failed")

    def test_ph_is_blocked_for_overlapping_bands(self):
        (self.outdir / "ph.in").touch()
        (self.outdir / "scf.out").write_text(
            "convergence has been achieved in   4 iterations\n"
            "highest occupied, lowest unoccupied level (ev): 2.5 2.4\n"
            "JOB DONE.\n"
        )
        ph = self.bin_dir / "ph.x"
        ph.write_text("#!/bin/sh\necho JOB DONE.\n")
        ph.chmod(0o755)

        state, ok = run_workflow(self.outdir, self.bin_dir, PSEUDO, np=1, stage="ph")

        self.assertFalse(ok)
        self.assertEqual(state["stages"]["ph"]["status"], "blocked")
        self.assertEqual(state["electronic_structure"]["gap_ev"], -0.1)
        self.assertFalse((self.outdir / "ph.out").exists())


class TestLoCleanup(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_lo_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_stale_lo_removed(self):
        from ase import Atoms

        from qe_inputs import generate_inputs

        atoms = Atoms("Li", positions=[(0, 0, 0)], cell=[3.0, 3.0, 3.0], pbc=True)
        cfg = {
            "qe_outdir": "./tmp",
            "ecutwfc": 60.0,
            "ecutrho": 240.0,
            "occupations": "fixed",
            "conv_thr": 1e-10,
            "mixing_beta": 0.5,
            "verbosity": "high",
            "tr2_ph": 1e-14,
            "kpt_density": 30.0,
        }
        generate_inputs(atoms, cfg, self.tmp, PSEUDO, "raman", lo_dirs=[(1.0, 0.0, 0.0)])
        self.assertTrue(list(self.tmp.glob("dynmat_lo*.in")))
        generate_inputs(atoms, cfg, self.tmp, PSEUDO, "raman", lo_dirs=[])
        self.assertFalse(list(self.tmp.glob("dynmat_lo*.in")))
        self.assertTrue((self.tmp / "inputs_manifest.json").exists())


_RELAX_OUT = (
    "bfgs converged in 2 scf cycles and 1 bfgs steps\n"
    "Begin final coordinates\n"
    "CELL_PARAMETERS (angstrom)\n"
    " 4.0 0.0 0.0\n 0.0 4.0 0.0\n 0.0 0.0 4.0\n"
    "ATOMIC_POSITIONS (crystal)\n"
    "Li 0.0 0.0 0.0\nF 0.5 0.5 0.5\n"
    "End final coordinates\n"
)


class TestRelaxStage(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_relax_"))
        self.outdir = self.tmp / "run"
        self.outdir.mkdir()
        (self.outdir / "relax.in").write_text("&CONTROL\n/\n")
        (self.outdir / "effective_config.json").write_text(
            json.dumps({"prefix": "t", "ecutwfc": 60.0, "ecutrho": 240.0}))
        (self.outdir / "inputs_manifest.json").write_text(
            json.dumps({"prefix": "t", "lo_dirs": []}))
        self.bin_dir = self.tmp / "bin"
        self.bin_dir.mkdir()
        (self.bin_dir / "dynmat.x").write_text("#!/bin/sh\nexit 0\n")
        (self.bin_dir / "dynmat.x").chmod(0o755)
        pw = self.bin_dir / "pw.x"
        pw.write_text(
            "#!/bin/sh\n"
            "cat <<'EOF'\n" + _RELAX_OUT + "EOF\n"
            "echo 'JOB DONE.'\nexit 0\n"
        )
        pw.chmod(0o755)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_relax_regenerates_scf(self):
        state, ok = run_workflow(
            self.outdir, self.bin_dir, PSEUDO, np=1, stage="relax")
        self.assertTrue(ok)
        self.assertTrue((self.outdir / "relax_manifest.json").exists())
        scf = (self.outdir / "scf.in").read_text()
        self.assertIn("ATOMIC_POSITIONS", scf)
        self.assertIn("F", scf)

    def test_unconverged_relax_does_not_regenerate_scf(self):
        scf = self.outdir / "scf.in"
        scf.write_text("original input\n")
        pw = self.bin_dir / "pw.x"
        pw.write_text(pw.read_text().replace(
            _RELAX_OUT.splitlines()[0], "The maximum number of steps has been reached."))
        state, ok = run_workflow(self.outdir, self.bin_dir, PSEUDO, np=1, with_relax=True)
        self.assertFalse(ok)
        self.assertEqual(state["stages"]["relax"]["status"], "failed")
        self.assertIn("maximum number of steps", state["stages"]["relax"]["error"])
        self.assertEqual(scf.read_text(), "original input\n")
        self.assertFalse((self.outdir / "scf.out").exists())


if __name__ == "__main__":
    unittest.main()


class TestScfArtifacts(unittest.TestCase):
    """QE may store the save directory as .dat (classic) or .hdf5."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_art_"))
        self.text = "&CONTROL\n prefix = 'LiF', outdir = './tmp'\n/\n"

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _save(self, names):
        save = self.tmp / "tmp" / "LiF.save"
        save.mkdir(parents=True)
        for name in names:
            (save / name).write_bytes(b"x")
        return save

    def test_classic_dat_build(self):
        self._save(["data-file-schema.xml", "charge-density.dat",
                    "wfc1.dat", "wfc2.dat"])
        art = _artifacts("scf", self.text, self.tmp)
        self.assertIsNotNone(art)
        self.assertTrue(any(k.endswith("charge-density.dat") for k in art))
        self.assertTrue(any(k.endswith("wfc2.dat") for k in art))

    def test_hdf5_build(self):
        self._save(["data-file-schema.xml", "charge-density.hdf5",
                    "wfc1.hdf5", "wfc2.hdf5"])
        art = _artifacts("scf", self.text, self.tmp)
        self.assertIsNotNone(art)
        self.assertEqual(len(art), 4)
        self.assertTrue(any(k.endswith("charge-density.hdf5") for k in art))
        self.assertTrue(any(k.endswith("wfc2.hdf5") for k in art))

    def test_missing_charge_density_returns_none(self):
        self._save(["data-file-schema.xml", "wfc1.dat"])
        self.assertIsNone(_artifacts("scf", self.text, self.tmp))


class TestHdf5UpstreamArtifacts(unittest.TestCase):
    """Regression: HDF5 save files must satisfy the scf->ph artifact gate."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qer_h5_"))
        self.outdir = self.tmp / "run"
        self.outdir.mkdir()
        (self.outdir / "scf.in").write_text(
            "&CONTROL\n prefix = 'LiF', outdir = './tmp'\n/\n")
        (self.outdir / "ph.in").write_text(
            "&inputph\n prefix = 'LiF', outdir = './tmp', fildyn = 'LiF.dyn'\n/\n")
        self.bin_dir = self.tmp / "bin"
        self.bin_dir.mkdir()
        (self.bin_dir / "pw.x").write_text(
            "#!/bin/sh\n"
            "mkdir -p tmp/LiF.save\n"
            "printf x > tmp/LiF.save/data-file-schema.xml\n"
            "printf x > tmp/LiF.save/charge-density.hdf5\n"
            "printf x > tmp/LiF.save/wfc1.hdf5\n"
            "echo JOB DONE.\n")
        (self.bin_dir / "ph.x").write_text(
            "#!/bin/sh\nprintf x > LiF.dyn\necho JOB DONE.\n")
        for name in ("pw.x", "ph.x"):
            (self.bin_dir / name).chmod(0o755)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_ph_runs_after_hdf5_scf(self):
        _, ok = run_workflow(self.outdir, self.bin_dir, PSEUDO, np=1, stage="scf")
        self.assertTrue(ok)
        _, ok = run_workflow(self.outdir, self.bin_dir, PSEUDO, np=1, stage="ph")
        self.assertTrue(ok)
        state = json.loads((self.outdir / "state.json").read_text())
        self.assertEqual(state["stages"]["scf"]["status"], "success")
        self.assertEqual(state["stages"]["ph"]["status"], "success")
