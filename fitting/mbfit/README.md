# fitting/mbfit — PIP fitting notebooks

Jupyter notebooks for fitting permutationally invariant polynomials (PIPs)
to the two- and three-body PMF training data.

## Notebooks

- `NP_Dimer_fitting.ipynb` — fits W2(d12) with a 7th-order PIP (Eq. 3 in paper)
- `NP_Trimer_fitting.ipynb` — fits ΔW3(d12,d13,d23) with a 5th-order PIP (Eq. 6)

## Framework credit

Both notebooks use the MB-Fit framework and MBX force evaluator
developed by the Paesani group (UCSD):

- Two-body basis: Babin, Leforestier & Paesani, *JCTC* 9, 5395 (2013)
- Three-body basis: Babin, Medders & Paesani, *JCTC* 10, 1599 (2014)
- MB-Fit software: https://github.com/paesanilab/MB-Fit
- MBX software: https://github.com/paesanilab/MBX

MB-Fit was originally developed for the MB-pol water potential at
quantum chemistry accuracy (CCSD(T)/CBS). These notebooks adapt it
to polymer-grafted nanoparticle systems.

## Key adaptations from MB-pol water

| Parameter | Water (MB-pol) | NPs (this work) |
|---|---|---|
| Monomer | H₂O (3 atoms, 3 sites) | Single-site particle (1 atom) |
| SMILES | `O` | `P` (phosphorus as proxy) |
| Distance range | ~1.5–7.5 Å | 6–14 σ (mesoscale) |
| Energy scale | up to 25 kcal/mol (ab initio) | ~20 kBT (PMF) |
| Electrostatics | charges + polarizability | all zero |
| Dispersion | C6, C8 terms | zero |
| Training data | CCSD(T) energies | PMF from CG MD |

## NP-specific parameter choices (in notebooks)

```python
min_d_2b = 6.0, max_d_2b = 14.0   # NP separation range (σ), not Å
polynomial_order = 7                # matches paper Eq. 3 / Supp. Table 1
num_training_configs = 121          # matches PMF grid spacing
mon_ids = ["p", "p"]               # must match p.json monomer name
```

The `bind_emax` value is set large: NPs have no bond dissociation limit
(PMF → 0 at large separation), unlike molecules where energy → +∞.

## Prerequisites

- MB-Fit installed: https://github.com/paesanilab/MB-Fit
- MBX installed: https://github.com/paesanilab/MBX
- Update `MBX_HOME` path in the last cell of each notebook

## Input

PMF data from `training/extract_pmf.py`:
- `W2.csv` → dimer notebook
- `DeltaW3_grid.csv` → trimer notebook

## Output

Fitted PIP parameter files consumed by MBX at runtime via
`simulation/3_assembly/mbx.json` and `p.json`.
