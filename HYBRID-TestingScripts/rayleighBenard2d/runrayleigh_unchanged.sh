#!/bin/bash
#SBATCH --job-name=rayleigh_unchanged-OMP_NUM_THREADS
#SBATCH --output=rayleigh_unchanged-OMP_NUM_THREADS.txt
#SBATCH --ntasks=OMP_NUM_THREADS

cd ../..
cd openLB_original/examples/thermal/rayleighBenard2d
make clean
make

OMP_NUM_THREADS=OMP_NUM_THREADS ./rayleighBenard2d