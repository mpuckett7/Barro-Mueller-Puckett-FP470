#!/bin/bash
#SBATCH --job-name=rayleigh1-MPI_NUM_TASKS
#SBATCH --output=rayleigh1-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS


cd ../../
cd OpenLB/examples/thermal/rayleighBenard2d

make

salloc -Qn MPI_NUM_TASKS mpirun ./rayleighBenard2d