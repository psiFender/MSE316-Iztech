#!/bin/bash

#SBATCH -p heavy
#SBATCH --nodes=1 
#SBATCH -J DynMatLAMM #çalışacak işin adı 
#SBATCH --cpus-per-task=40
##SBATCH --time=0 #uygulamanin en fazla calisma suresi 
#SBATCH --output=slurm-%j.out 
#SBATCH --error=slurm-%j.err 


export FI_PROVIDER=sockets
#export I_MPI_OFI_PROVIDER=mlx
export OMP_NUM_THREADS=4
ulimit -s unlimited
ulimit -l unlimited
ulimit -m unlimited

#export LD_LIBRARY_PATH=/usr/lib64/:$LD_LIBRARY_PATH
#module load lammps/29Oct20
#export /opt/source/lammps-29Oct20/potentials/
source /opt/intel/oneapi/setvars.sh

mpirun -np 10  /opt/lammps-29Oct20/src/lmp_intel_cpu_intelmpi -in input.lammps > lammps_out.txt


