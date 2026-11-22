"""MPRelaxSet (pymatgen c4998d92) -> QE parameter mapping.

Shared by ``qe-mprelax-mapping`` and the benchmark skills so every run uses the
same conversion.  See ``qe-mprelax-mapping/references/vasp-qe-mapping.md``.
"""
from __future__ import annotations

RY_TO_EV = 13.605693122994
ENCUT_EV = 520.0
EDIFF_PER_ATOM_EV = 5.0e-5
SIGMA_EV = 0.05
NELM = 100
NSTEP = 99
FORC_CONV_THR_RY_BOHR = 5.0e-4
ETOT_CONV_THR_RY = 1.0e-5

MP_SETTINGS = {
    "smearing": "gaussian",
    "mixing_beta": 0.7,
    "electron_maxstep": NELM,
    "nstep": NSTEP,
    "forc_conv_thr": FORC_CONV_THR_RY_BOHR,
    "etot_conv_thr": ETOT_CONV_THR_RY,
    "nspin": 1,
    "ion_dynamics": "bfgs",
    "cell_dynamics": "bfgs",
}


def mp_kpoints(structure) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    from pymatgen.io.vasp.sets import MPRelaxSet

    kp = MPRelaxSet(structure).kpoints
    if kp is None:
        raise RuntimeError("MPRelaxSet returned no k-points")
    grid = tuple(int(x) for x in kp.kpts[0])
    style = str(kp.style).lower()
    shift = (0, 0, 0) if "gamma" in style else (1, 1, 1)
    return grid, shift


def mp_cutoffs(structure, sssp: dict) -> tuple[float, float]:
    elements = [str(e) for e in structure.composition.elements]
    encut_ry = ENCUT_EV / RY_TO_EV
    ecutwfc = max(encut_ry, max(sssp[e]["cutoff_wfc"] for e in elements))
    ecutrho = max(8.0 * ecutwfc, max(sssp[e]["cutoff_rho"] for e in elements))
    return round(ecutwfc, 1), round(ecutrho, 1)


def mp_electrons(natoms: int) -> tuple[float, float]:
    conv_thr = max(1.0e-8, EDIFF_PER_ATOM_EV * natoms / RY_TO_EV)
    degauss = SIGMA_EV / RY_TO_EV
    return conv_thr, degauss


def mapping_for(structure, sssp: dict) -> dict:
    grid, shift = mp_kpoints(structure)
    ecutwfc, ecutrho = mp_cutoffs(structure, sssp)
    conv_thr, degauss = mp_electrons(len(structure))
    return {
        "formula": structure.composition.reduced_formula,
        "natoms": len(structure),
        "ecutwfc_Ry": ecutwfc,
        "ecutrho_Ry": ecutrho,
        "conv_thr_Ry": conv_thr,
        "smearing": MP_SETTINGS["smearing"],
        "degauss_Ry": degauss,
        "k_grid": grid,
        "k_shift": shift,
        "electron_maxstep": MP_SETTINGS["electron_maxstep"],
        "nstep": MP_SETTINGS["nstep"],
        "forc_conv_thr_Ry_bohr": MP_SETTINGS["forc_conv_thr"],
        "etot_conv_thr_Ry": MP_SETTINGS["etot_conv_thr"],
        "nspin": MP_SETTINGS["nspin"],
        "calculation": "vc-relax",
        "ion_dynamics": MP_SETTINGS["ion_dynamics"],
        "cell_dynamics": MP_SETTINGS["cell_dynamics"],
    }
