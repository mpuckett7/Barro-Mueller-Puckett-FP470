#!/bin/bash
#SBATCH --job-name=rayleigh_4
#SBATCH --output=rayleigh_4.txt

export OMP_NUM_THREADS=4
./rayleighBenard2d
