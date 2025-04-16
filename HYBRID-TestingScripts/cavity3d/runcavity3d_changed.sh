#!/bin/bash
#SBATCH --job-name=cavity3d_changed-OMP_NUM_THREADS
#SBATCH --output=cavity3d_changed-OMP_NUM_THREADS.txt
#SBATCH --ntasks=OMP_NUM_THREADS

cd ../../
cd OpenLB/examples/laminar/cavity3d
make clean
make

OMP_NUM_THREADS=OMP_NUM_THREADS ./cavity3d