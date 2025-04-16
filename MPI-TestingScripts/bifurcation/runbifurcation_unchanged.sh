#!/bin/bash
#SBATCH --job-name=bifurcation3d_unchanged-MPI_NUM_TASKS
#SBATCH --output=bifurcation3d_unchanged-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../../

# make clean; make

cd openLB_original/examples/particles/bifurcation3d/eulerLagrange

make

salloc -Qn MPI_NUM_TASKS mpirun ./bifurcation3d