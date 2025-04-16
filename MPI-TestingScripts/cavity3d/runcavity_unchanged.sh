#!/bin/bash
#SBATCH --job-name=cavity3d_unchanged-MPI_NUM_TASKS
#SBATCH --output=cavity3d_unchanged-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../..

# make clean; make

cd openLB_original/examples/laminar/cavity3d

make

salloc -Qn MPI_NUM_TASKS mpirun ./cavity3d