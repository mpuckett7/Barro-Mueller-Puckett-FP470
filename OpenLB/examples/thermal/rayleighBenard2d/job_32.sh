#!/bin/bash
#SBATCH --job-name=rayleigh_32
#SBATCH --output=rayleigh_32.txt

export OMP_NUM_THREADS=32
./rayleighBenard2d
