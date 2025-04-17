#!/bin/bash
#SBATCH --job-name=rayleigh_2
#SBATCH --output=rayleigh_2.txt

export OMP_NUM_THREADS=2
./rayleighBenard2d
