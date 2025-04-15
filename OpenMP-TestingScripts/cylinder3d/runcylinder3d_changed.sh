#!/bin/bash
#SBATCH --job-name=cylinder3d_changed-OMP_NUM_THREADS
#SBATCH --output=cylinder3d_changed-OMP_NUM_THREADS.txt
#SBATCH --ntasks=OMP_NUM_THREADS

cd ../../
cd OpenLB/examples/laminar/cylinder3d
make clean
make

OMP_NUM_THREADS=OMP_NUM_THREADS ./cylinder3d