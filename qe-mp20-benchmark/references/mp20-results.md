# MP20 results

QE 7.6 + SSSP 1.3.0 PBE efficiency, GPU build, vc-relax from the MP CIF, dense
nscf band gap, QE elemental references, MP2020 corrections. MP/VASP references
from the Materials Project.

## MP20-5 (initial systems)

Systems: AlAs, ZnTe, AgI, AlSb, MgTe.

| Metric | Value |
|---|---:|
| Band-gap MAE | **0.117 eV** |
| Band-gap RMSE | 0.177 eV |
| Max gap error | 0.341 eV (AgI) |
| Volume MAE | **0.166 %** |
| Raw formation-energy MAE | 0.135 eV/atom |

ZnTe/AlSb/MgTe gap errors <0.03 eV; AgI 0.34 eV, AlAs 0.20 eV. CdTe `mp-406`
was excluded because its MP gap (1.2285 eV) contradicts the other CdTe
polymorphs (~0.62 eV).

## MP20-20 (expanded systems)

Systems span III-V, II-VI, oxides, nitrides, halides and metal/covalent
compounds.

| Metric | SSSP | MP-like |
|---|---:|---:|
| Band-gap MAE | 0.1686 eV | **0.1680 eV** |
| Band-gap RMSE | 0.2410 eV | 0.2415 eV |
| Max gap error | 0.5514 eV (PbS) | 0.5528 eV |
| Volume MAE | 0.4347 % | **0.4053 %** |
| Max volume error | 1.1822 % | 1.1239 % |
| Raw formation-energy MAE | 0.1724 eV/atom | 0.1730 eV/atom |
| **MP2020 formation-energy MAE** | 0.0145 eV/atom | **0.0131 eV/atom** |
| MP2020 max formation-energy error | 0.0505 eV/atom | 0.0449 eV/atom |

Main gap outliers: PbS/PbSe ~ -0.5 eV, InN -0.48 eV, AgI -0.34 eV, GaN
-0.32 eV.

Convergence and precision checks: AlAs/AgI gaps change <0.001 eV over k = 6-14
and ecut = 35-65 Ry; SSSP precision moves the gaps by only ~6 meV. The residual
gap disagreement is therefore not a QE numerical error.

## MPRelaxSet-like vs SSSP

| Metric | Old SSSP | New MP-like |
|---|---:|---:|
| MP20 gap MAE | 0.1686 eV | 0.1680 eV |
| MP20 volume MAE | 0.4347 % | 0.4053 % |
| MP20 MP2020 Ef MAE | 0.0145 eV/atom | 0.0131 eV/atom |
| MPTraj force MAE | 0.1428 eV/Å | 0.1454 eV/Å |
| MPTraj stress MAE | 2.467e-3 eV/Å³ | 2.481e-3 eV/Å³ |

The two configurations give nearly identical MP20 accuracy; the MP-like one is
slightly better for volume and formation energy. MPTraj forces/stress are
insensitive to the configuration.

## Interpretation

- Same-basin geometry: volume MAE 0.17-0.43 % (max ~1.2 %) — good.
- MP2020-corrected formation energy: ~13 meV/atom — good for thermodynamic
  screening.
- Band gap: 0.168 eV MAE — trend-level only until a controlled same-structure,
  same-functional VASP-PBE reference is available.
- Raw total energies are not comparable (pseudopotential zeros differ).
