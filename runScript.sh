#!/bin/bash
#SBATCH --ntasks=15
#SBATCH --nodes=1

module load gcc/7.2.0
module load openmpi-1.8/gcc
mpirun -np 15 ./MPIGillespie.exe