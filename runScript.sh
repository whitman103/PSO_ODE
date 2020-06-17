#!bin/bash
#SBATCH --nodes=15

module load gcc/7.2.0
module load openmpi-1.8/gcc
mpic++ -O2 -g MPIGillespie.cpp fuzzyDef.cpp GillespieFunctions.cpp -o MPIGillespie.exe
mpirun -np 15 ./MPIGillespie.exe