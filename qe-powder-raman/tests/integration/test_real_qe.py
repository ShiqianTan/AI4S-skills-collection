"""Optional real-QE integration smoke test (opt-in, default skipped).

Run a tiny LiF job through the full CLI (`all --execute`) and check the negative
control (rocksalt LiF has no first-order Raman). This test needs a working
Quantum ESPRESSO build, so it only runs when both are true:

  * environment variable ``QE_REAL_TESTS=1`` is set, and
  * ``pw.x`` is found in ``QE_BIN`` (or the skill's default bin dir).

Typical use on the validated local build:

    QE_REAL_TESTS=1 QE_BIN=/path/to/install_cpu/bin \
        python tests/run_tests.py

It uses ``kpt_density: 6`` (a coarse mesh) so it stays reasonably fast; it is a
workflow smoke test, not a converged physics reference.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SKILL = Path(__file__).resolve().parents[2]
SCRIPTS = SKILL / "scripts"
EXAMPLE = SKILL / "examples" / "LiF_mp-1138.cif"
sys.path.insert(0, str(SCRIPTS))

try:
    from qe_common import DEFAULT_BIN_DIR
except Exception:  # noqa: BLE001
    DEFAULT_BIN_DIR = None


def _bin_dir():
    env = os.environ.get("QE_BIN")
    if env:
        return Path(env).expanduser()
    if DEFAULT_BIN_DIR is not None:
        return Path(DEFAULT_BIN_DIR).expanduser()
    return None


@unittest.skipUnless(
    os.environ.get("QE_REAL_TESTS") == "1",
    "set QE_REAL_TESTS=1 (and QE_BIN) to run the real-QE smoke test",
)
class TestRealQeLiF(unittest.TestCase):
    def test_lif_negative_control(self):
        bin_dir = _bin_dir()
        if bin_dir is None or not (bin_dir / "pw.x").exists():
            self.skipTest("pw.x not found (set QE_BIN to a QE install)")
        if not EXAMPLE.exists():
            self.skipTest("example CIF missing: %s" % EXAMPLE)

        tmp = Path(tempfile.mkdtemp(prefix="qer_realqe_"))
        try:
            cfg = tmp / "cfg.yaml"
            cfg.write_text("kpt_density: 6.0\n")
            run = tmp / "run"
            cmd = [
                sys.executable, str(SCRIPTS / "qe_raman.py"), "all",
                str(EXAMPLE),
                "--config", str(cfg),
                "--outdir", str(run),
                "--prefix", "LiF",
                "--bin-dir", str(bin_dir),
                "--np", "1",
                "--execute",
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            self.assertEqual(
                proc.returncode, 0,
                "all --execute failed:\nstdout:\n%s\nstderr:\n%s"
                % (proc.stdout[-4000:], proc.stderr[-4000:]),
            )

            summary = json.loads((run / "summary.json").read_text())
            self.assertNotEqual(summary["status"], "failed")
            self.assertEqual(summary["stats"]["total"], 6)  # 3N, N=2
            # LiF is a negative control: no first-order Raman-active mode.
            self.assertEqual(summary["stats"]["raman_active"], [])
            self.assertEqual(summary["stats"]["sum_raman"], 0.0)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
