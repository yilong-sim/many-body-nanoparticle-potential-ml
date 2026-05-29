# param_template.py — fully documented parameter template
# Copy to param.py and edit before running generate_lammps_data.py.
#
# SYSTEM VARIANTS USED IN PAPER (npj Comp. Mat. 2023):
#   Bare NP (globular aggregate):  diameter=6, grftdens=0.0
#   Grafted 0.15 (2D sheets):      diameter=6, grftdens=0.15
#   Grafted 0.3  (1D strings):     diameter=6, grftdens=0.3
#   Grafted 0.4  (dispersed):      diameter=6, grftdens=0.4

# ── Job identification ────────────────────────────────────────────
pbsname = '6SigTemplate'   # label for PBS/SLURM job name

# ── Reproducibility ───────────────────────────────────────────────
rseed = 1                  # integer random seed

# ── Nanoparticle geometry ─────────────────────────────────────────
diameter = 6               # NP diameter in units of σ
                           # Core mass auto-set to diameter^3
                           # Set to 0 for pure polymer system (no NP)

# ── Polymer grafting ──────────────────────────────────────────────
grftdens = 0.3             # Grafting density (chains/σ²)
                           # Number of grafts = round(π × diameter² × grftdens)
                           # Set to 0 for bare NP

ltether = 10               # Graft chain length (number of beads)
                           # Paper uses Lg = 10 throughout

# ── Matrix polymer ────────────────────────────────────────────────
lmatrix = 20               # Matrix chain length (number of beads)
                           # Paper uses Lm = 20 throughout
                           # lmatrix > ltether → matrix penetrates grafts (wetting)

# ── Attractive interactions ───────────────────────────────────────
attraction = 0             # Number of sticky beads at graft chain end
                           # 0 = standard LJ grafts (used in paper)

matrixattraction = 0       # Number of sticky beads at matrix chain end
                           # 0 = standard LJ matrix (used in paper)

# ── Thermodynamic conditions ──────────────────────────────────────
systemdens = 0.85          # Bead number density (beads/σ³)
                           # 0.85 = polymer melt conditions
                           # Temperature fixed at T = ε/kB in LAMMPS input

wtprcnt = 9.8              # Target NP weight percent (%)
                           # Determines number of matrix chains added
                           # Actual wt% printed during generation

# ── Simulation length ─────────────────────────────────────────────
runtime = 27               # Simulation runtime (millions of timesteps)
                           # With timestep 0.002τ: 27M steps = 54,000τ

# ── Box replication ───────────────────────────────────────────────
lammpsreplicate = '3 2 2'  # Replication in x y z directions
                           # Total NPs = product = 12 NPs for PMF calc
                           # Use '1 1 1' for single-NP test runs
