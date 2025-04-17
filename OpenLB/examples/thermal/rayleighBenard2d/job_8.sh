#!/bin/bash
#SBATCH --job-name=rayleigh_8
#SBATCH --output=rayleigh_8.txt

export OMP_NUM_THREADS=8
./rayleighBenard2d
