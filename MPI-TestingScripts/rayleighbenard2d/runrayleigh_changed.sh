#!/bin/bash
#SBATCH --job-name=rayleigh_changed-MPI_NUM_TASKS
#SBATCH --output=rayleigh_changed-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS
#SBATCH --ntasks-per-node=4

cd ../../
cd OpenLB/examples/thermal/rayleighBenard2d

make

salloc -Qn MPI_NUM_TASKS mpirun ./rayleighBenard2d