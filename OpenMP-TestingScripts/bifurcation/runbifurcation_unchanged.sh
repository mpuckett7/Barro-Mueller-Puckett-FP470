#!/bin/bash
#SBATCH --job-name=bifurcation3d_unchanged-OMP_NUM_THREADS
#SBATCH --output=bifurcation3d_unchanged-OMP_NUM_THREADS.txt
#SBATCH --ntasks=OMP_NUM_THREADS

cd ../..
cd openLB_original/examples/particles/bifurcation3d/eulerLagrange
make clean
make

OMP_NUM_THREADS=OMP_NUM_THREADS ./bifurcation3d