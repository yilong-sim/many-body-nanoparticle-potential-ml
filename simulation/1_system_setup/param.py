# param.py — System parameters for LAMMPS data file generation
# Pass to generate_lammps_data.py by placing in the working directory.

pbsname        = '6SigTemplate'  # job name label
rseed          = 1               # random seed for simulation

# Nanoparticle geometry
diameter       = 6               # NP diameter (σ); radius = diameter/2
                                 # NP core mass = diameter^3 = 216m

# Polymer grafting
grftdens       = 0.3             # grafting density (chains/σ²)
                                 # 0.0 = bare NP, 0.15 = weak brush,
                                 # 0.3 = moderate brush, 0.4 = dense brush
ltether        = 10              # graft chain length (beads)

# Matrix polymer
lmatrix        = 20              # matrix chain length (beads)
                                 # lmatrix > ltether → wetting regime

# Attractive interactions (set 0 for standard grafted NP system)
attraction       = 0             # sticky beads per graft end
matrixattraction = 0             # sticky beads per matrix chain end

# System thermodynamic conditions
systemdens    = 0.85             # bead number density (beads/σ³)
                                 # 0.85 = melt-like conditions
wtprcnt       = 9.8              # target NP weight percent

# Simulation runtime
runtime       = 27               # simulation length (millions of timesteps)
                                 # timestep = 0.002τ → total = 54×10⁶ τ units

# Box replication — controls number of NPs in system
lammpsreplicate = '3 2 2'        # (x y z) replication of unit cell
                                 # total NPs = product of three integers = 12
