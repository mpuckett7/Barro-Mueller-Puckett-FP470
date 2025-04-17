#!/bin/bash
#SBATCH --job-name=rayleigh_1
#SBATCH --output=rayleigh_1.txt

export OMP_NUM_THREADS=1
./rayleighBenard2d
