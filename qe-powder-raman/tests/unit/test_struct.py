"""Unit tests for qe_struct (structure checks + refinement)."""
import sys
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
sys.path.insert(0, str(SCRIPTS))

from qe_struct import (  # noqa: E402
    StructureCheckResult,
    _check_layers,
    _check_partial_occupancy,
    analyze_layers,
    check_structure,
    parse_final_structure,
    refine_structure,
    relax_metrics,
    validate_structure,
)


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
        self.assertIn("fallback_reason", info)

    def test_lif_no_change(self):
        from qe_struct import read_structure

        atoms = read_structure(EXAMPLES / "LiF_mp-1138.cif")
        _, info = refine_structure(atoms, symprec=1e-3)
        self.assertFalse(info["standardized"])


class TestValidateStructure(unittest.TestCase):
    def _atoms(self, **info):
        from ase import Atoms

        a = Atoms(
            "LiF",
            positions=[(0, 0, 0), (2.0, 2.0, 2.0)],
            cell=[4.0, 4.0, 4.0],
            pbc=True,
        )
        a.info.update(info)
        return a

    def test_partial_occupancy_fail(self):
        a = self._atoms(occupancy={"0": {"Li": 0.5}, "1": {"F": 1.0}})
        res = StructureCheckResult()
        _check_partial_occupancy(a, res)
        self.assertTrue(
            any(lvl == "FAIL" and key == "occupancy" for lvl, key, _ in res.messages)
        )

    def test_mixed_site_fail(self):
        a = self._atoms(occupancy={"0": {"Li": 0.5, "F": 0.5}, "1": {"F": 1.0}})
        res = StructureCheckResult()
        _check_partial_occupancy(a, res)
        self.assertTrue(
            any(lvl == "FAIL" and key == "mixed_occupancy" for lvl, key, _ in res.messages)
        )

    def test_pbc_fail(self):
        from ase import Atoms

        a = Atoms("Li", positions=[(0, 0, 0)], cell=[4.0, 4.0, 4.0], pbc=False)
        res = StructureCheckResult()
        validate_structure(a, res)
        self.assertTrue(
            any(lvl == "FAIL" and key == "non_3d" for lvl, key, _ in res.messages)
        )


_RELAX_BLOCK = """Begin final coordinates
     new unit-cell volume = 100.0
CELL_PARAMETERS (angstrom)
    4.000000   0.000000   0.000000
    0.000000   4.000000   0.000000
    0.000000   0.000000   4.000000
ATOMIC_POSITIONS (crystal)
Li  0.000000  0.000000  0.000000
F   0.500000  0.500000  0.500000
End final coordinates
"""


class TestParseFinalStructure(unittest.TestCase):
    def test_parse(self):
        atoms = parse_final_structure(_RELAX_BLOCK)
        self.assertIsNotNone(atoms)
        self.assertEqual(len(atoms), 2)
        self.assertEqual(list(atoms.get_chemical_symbols()), ["Li", "F"])
        self.assertAlmostEqual(float(atoms.get_cell()[0][0]), 4.0, places=6)

    def test_missing_block(self):
        self.assertIsNone(parse_final_structure("no coordinates here"))

    def test_latest_alat_block_scales_cell_and_positions(self):
        latest = _RELAX_BLOCK.replace("(angstrom)", "(alat= 2.0)").replace("(crystal)", "(alat)")
        atoms = parse_final_structure(_RELAX_BLOCK + latest)
        self.assertAlmostEqual(atoms.cell[0, 0], 8 * 0.529177210903)
        self.assertAlmostEqual(atoms.positions[1, 0], 0.529177210903)

    def test_alat_without_scale_is_rejected(self):
        self.assertIsNone(parse_final_structure(_RELAX_BLOCK.replace("(angstrom)", "(alat)")))

    def test_relax_metrics(self):
        m = relax_metrics("step 1\nTotal force = 0.0500\nTotal force = 0.0123\n")
        self.assertAlmostEqual(m["final_force"], 0.0123)


class TestLayers(unittest.TestCase):
    def test_mos2_layered(self):
        res = check_structure(EXAMPLES / "MoS2_mp-2815.cif")
        self.assertFalse(res.has_fail)
        self.assertIsNotNone(res.layers)
        self.assertTrue(res.layers["is_layered"])
        self.assertEqual(res.layers["layer_count"], 2)
        self.assertEqual(res.layers["gap_elements"], ["S", "S"])
        keys = [k for _, k, _ in res.messages]
        self.assertIn("layered_gap", keys)

    def test_lif_not_layered(self):
        res = check_structure(EXAMPLES / "LiF_mp-1138.cif")
        self.assertFalse(res.layers["is_layered"])
        self.assertFalse(res.layers["is_slab"])
        keys = [k for _, k, _ in res.messages]
        self.assertNotIn("layered_gap", keys)
        self.assertNotIn("slab_vacuum", keys)

    def test_synthetic_slab_vacuum(self):
        from ase import Atoms

        atoms = Atoms("Si2", positions=[(0, 0, 0), (1.35, 1.35, 0)],
                      cell=[3.0, 3.0, 20.0], pbc=True)
        info = analyze_layers(atoms)
        self.assertTrue(info["is_slab"])
        self.assertGreater(info["vacuum_thickness"], 8.0)
        res = StructureCheckResult()
        _check_layers(atoms, res)
        self.assertTrue(any(k == "slab_vacuum" for _, k, _ in res.messages))


if __name__ == "__main__":
    unittest.main()
