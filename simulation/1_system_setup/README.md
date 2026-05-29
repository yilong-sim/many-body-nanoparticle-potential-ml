# 1_system_setup

Generates the LAMMPS data file (atom positions, bonds, masses) for the
polymer-grafted NP system used in PMF calculations.

## Files

- `param.py` — system parameters (edit this before running)
- `param_template.py` — fully documented template with all options
- `generate_lammps_data.py` — generates `lammps.input` from `param.py`

## Usage

```bash
# 1. Edit param.py to set your system
# 2. Run the generator
python generate_lammps_data.py

# Output: lammps.input (gitignored — regenerate as needed)
```

## What it generates

A LAMMPS data file containing:
- NP core atom (type 3, mass = diameter³)
- NP shell beads placed on Fibonacci sphere (type 1) — graft anchor points
- Graft chain beads (type 2), extending radially outward
- Matrix chain beads (type 4), placed on a grid

## Atom type legend

| Type | Description | Mass |
|---|---|---|
| 1 | NP shell bead (graft anchor, Fibonacci placement) | 1m |
| 2 | Graft chain bead | 1m |
| 3 | NP core (rigid) | diameter³ × m |
| 4 | Matrix chain bead | 1m |
| 5 | Matrix end bead (attraction variants) | 1m |
