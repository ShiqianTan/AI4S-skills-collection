"""Unit tests for qe_struct (structure checks + refinement)."""
import sys
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
sys.path.insert(0, str(SCRIPTS))

from qe_struct import check_structure, refine_structure  # noqa: E402


class TestCheckStructure(unittest.TestCase):
    def test_lif_ok(self):
        res = check_structure(EXAMPLES / "LiF_mp-1138.cif")
        self.assertFalse(res.has_fail)
        self.assertEqual(res.n_atoms, 2)
        self.assertEqual(sorted(set(res.elements)), ["F", "Li"])
        self.assertEqual(res.spacegroup_number, 225)
        self.assertTrue(res.n_sym_ops > 0)

    def test_cds_hex_noise_warn(self):
        # original (un-refined) CIF has gamma=119.99998586 deg
        res = check_structure(EXAMPLES / "CdS_mp-672.cif")
        self.assertFalse(res.has_fail)
        self.assertTrue(res.has_warn)
        keys = [k for _, k, _ in res.messages]
        self.assertIn("hex_angle_noise", keys)


class TestRefineStructure(unittest.TestCase):
    def test_cds_cleanup(self):
        import ase  # noqa: F401
        from qe_struct import read_structure

        atoms = read_structure(EXAMPLES / "CdS_mp-672.cif")
        refined, info = refine_structure(atoms, symprec=1e-3)
        self.assertEqual(len(refined), len(atoms))
        self.assertTrue(info["standardized"])
        # gamma must be cleaned to exactly 120.0
        g = info["cellpar"][5]
        self.assertAlmostEqual(g, 120.0, places=5)
        self.assertEqual(info["refined"]["spacegroup_number"], 186)

    def test_lif_no_change(self):
        from qe_struct import read_structure

        atoms = read_structure(EXAMPLES / "LiF_mp-1138.cif")
        _, info = refine_structure(atoms, symprec=1e-3)
        self.assertFalse(info["standardized"])


if __name__ == "__main__":
    unittest.main()
