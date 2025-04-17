#!/bin/bash
#SBATCH --job-name=rayleigh_changed-MPI_NUM_TASKS
#SBATCH --output=rayleigh_changed-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS


cd ../../
cd OpenLB/examples/thermal/rayleighBenard2d

make

salloc -Qn MPI_NUM_TASKS --ntasks-per-node=4 mpirun ./rayleighBenard2d