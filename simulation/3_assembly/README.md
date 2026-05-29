# 3_assembly

LAMMPS input and MBX configuration files for running assembly simulations
using the fitted many-body potential.

## Files

- `assembly_mbx.lammps` — LAMMPS run script using MBX force field
- `mbx.json` — MBX runtime configuration (cutoffs, solver settings)
- `p.json` — NP monomer definition for MBX
- `job.sh` — SLURM submission template (update account and MBX path)

## How this differs from 2_pmf_calculation

| | PMF calculation | Assembly simulation |
|---|---|---|
| System | ~22,000 atoms (explicit polymer) | 125 NP centers only |
| Force field | `lj/expand` (CG explicit) | `pair_style mbx` (fitted PIP) |
| Timestep | 0.002τ | 0.02τ (10× larger) |
| Thermostat | NVT Nosé-Hoover | Langevin |
| Runtime | ~2.75M steps | 200M steps |

The 10× timestep increase is enabled by the implicit polymer treatment —
no polymer degrees of freedom means no fast bond vibrations to constrain.

## MBX configuration notes

**`mbx.json`**
- `twobody_cutoff: 12.0` matches Ro = 12σ in the two-body PIP switching function
- `threebody_cutoff: 10.0` matches Ro = 10σ in the three-body PIP switching function
- `dipole_method: "cg"` — conjugate gradient solver for induced dipoles
  (dipoles are zero here since pol = 0.0 in p.json, so this is a no-op)

**`p.json`**
Defines NPs as single-site monomers with all electrostatic and dispersion
terms set to zero. This is the key adaptation from the MB-pol water
potential: water has charges, polarizability, and dispersion; NPs have none.
Only the fitted short-range PIP terms (W2 + ΔW3) contribute to the energy.

```json
"charges" : [0.0],    // NPs are electrically neutral
"pol"     : [0.0],    // no polarizability — no induction energy
"c6lr"    : [0.0]     // no long-range dispersion
```

## Usage

```bash
# 1. Update job.sh: set your account and path to lmp_mpi_mbx
# 2. Ensure mbx.json and p.json are in the run directory
# 3. Submit
sbatch job.sh
```

## Langevin damping

Damping parameter 0.0864 (mσ²/ε)^(1/2) was calibrated so that NP
self-diffusivity from Langevin dynamics matches that from explicit CG MD:
D_xyz = 0.0005 σ²/τ in both cases. See Supplementary Fig. 10 of the paper.
