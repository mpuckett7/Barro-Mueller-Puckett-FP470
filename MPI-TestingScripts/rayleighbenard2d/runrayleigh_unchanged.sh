#!/bin/bash
#SBATCH --job-name=rayleigh_unchanged-MPI_NUM_TASKS
#SBATCH --output=rayleigh_unchanged-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../..

# make clean; make

cd openLB_original/examples/thermal/rayleighBenard2d

make

salloc -Qn MPI_NUM_TASKS mpirun ./rayleighBenard2d