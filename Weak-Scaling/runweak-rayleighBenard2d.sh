#!/bin/bash
# #SBATCH --job-name=RAY_OMP_TESTING
# #SBATCH --nodes=1

cd ../../
cd OpenLB/
make clean; make
cd examples/thermal/rayleighBenard2d
make

for t in 2 5 10 20; do
    OMP_NUM_THREADS=16 salloc ./rayleighBenard2d -r $t &> weak_ray_r=$t.txt
done