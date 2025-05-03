#!/bin/bash
#SBATCH --job-name=bifurcation3d-hybrid-NUM_TASKS
#SBATCH --output=bifurcation3d-hybrid-NUM_TASKS.txt
#SBATCH --ntasks=NUM_TASKS

cd ../../
cd OpenLB/

make clean; make

cd examples/particles/bifurcation3d/eulerLagrange

make

OMP_NUM_THREADS=NUM_TASKS salloc -Qn 2 mpirun ./bifurcation3d