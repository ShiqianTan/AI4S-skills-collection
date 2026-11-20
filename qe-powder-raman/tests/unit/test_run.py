"""Tests for qe_run.run_workflow (env pass-through, skip/resume) + LO cleanup."""
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
PSEUDO = Path(__file__).resolve().parents[2] / "pseudo"
sys.path.insert(0, str(SCRIPTS))

from qe_run import run_workflow  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
