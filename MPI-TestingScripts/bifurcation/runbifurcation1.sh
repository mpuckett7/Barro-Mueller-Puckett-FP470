#!/bin/bash
#SBATCH --job-name=bifurcation3d1-MPI_NUM_TASKS
#SBATCH --output=bifurcation3d1-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../../
cd OpenLB/

make clean; make

cd examples/particles/bifurcation3d/eulerLagrange

make

salloc -Qn MPI_NUM_TASKS mpirun ./bifurcation3d