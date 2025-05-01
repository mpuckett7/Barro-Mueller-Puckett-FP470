#!/bin/bash
#SBATCH --job-name=rayleigh-hybrid-NUM_TASKS
#SBATCH --output=rayleigh-hybrid-NUM_TASKS.txt
#SBATCH --ntasks=NUM_TASKS

cd ../../
cd OpenLB/

make clean; make 

cd examples/thermal/rayleighBenard2d

make

OMP_NUM_THREADS=NUM_TASKS salloc -Qn 4 mpirun ./rayleighBenard2d