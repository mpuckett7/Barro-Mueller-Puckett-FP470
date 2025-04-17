#!/bin/bash
# #SBATCH --job-name=RAY_OMP_TESTING
# #SBATCH --nodes=1

cd ../../
cd OpenLB/
make clean; make
cd examples/thermal/rayleighBenard2d
make

for t in 1 2 4 8 16 32; do
    OMP_NUM_THREADS=$t salloc ./rayleighBenard2d
done