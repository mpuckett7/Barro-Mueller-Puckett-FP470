#!/bin/bash
#SBATCH --job-name=rayleigh_changed-OMP_NUM_THREADS
#SBATCH --output=rayleigh_changed-OMP_NUM_THREADS.txt
#SBATCH --ntasks=OMP_NUM_THREADS

cd ../../
cd OpenLB/examples/thermal/rayleighBenard2d
make clean
make

OMP_NUM_THREADS=OMP_NUM_THREADS ./rayleighBenard2d