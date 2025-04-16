#!/bin/bash
#SBATCH --job-name=cavity3d_changed-MPI_NUM_TASKS
#SBATCH --output=cavity3d_changed-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../../
cd OpenLB/examples/laminar/cavity3d

make

salloc -Qn MPI_NUM_TASKS mpirun ./cavity3d