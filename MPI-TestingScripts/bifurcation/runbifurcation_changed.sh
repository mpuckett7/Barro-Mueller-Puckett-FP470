#!/bin/bash
#SBATCH --job-name=bifurcation3d_changed-MPI_NUM_TASKS
#SBATCH --output=bifurcation3d_changed-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../../
cd OpenLB/examples/particles/bifurcation3d/eulerLagrange

make

salloc -Qn MPI_NUM_TASKS mpirun ./bifurcation3d