#!/bin/bash
#SBATCH --partition=standart
#SBATCH --job-name=vasp-graphene
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=2
#SBATCH --share
#SBATCH --output=vasp_%j.out
#SBATCH --error=vasp_%j.err

module load vasp/5.4
module load intel/17.0.5.239


ulimit -s unlimited

mpirun -np $SLURM_NTASKS vasp_std 
