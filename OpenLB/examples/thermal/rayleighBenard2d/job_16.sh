#!/bin/bash
#SBATCH --job-name=rayleigh_16
#SBATCH --output=rayleigh_16.txt

export OMP_NUM_THREADS=16
./rayleighBenard2d
