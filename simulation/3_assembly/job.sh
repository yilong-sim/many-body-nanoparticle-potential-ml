#!/bin/bash
# SLURM submission script for MBX assembly simulation
# Update the three lines marked UPDATE before submitting.

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 8:00:00
#SBATCH -J np_assembly
#SBATCH -A YOUR_ACCOUNT_HERE       # UPDATE: your allocation account
#SBATCH --export=ALL

module load openmpi/4.0.4 fftw gcc

export OMP_NUM_THREADS=9

# UPDATE: set path to your lmp_mpi_mbx executable
LMP_EXEC=/path/to/lmp_mpi_mbx

srun --mpi=pmi2 -n 4 $LMP_EXEC -in assembly_mbx.lammps

exit
