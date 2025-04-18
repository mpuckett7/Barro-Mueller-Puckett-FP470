#!/bin/bash
# #SBATCH --job-name=RAY_OMP_TESTING
# #SBATCH --nodes=1

cd ../
cd OpenLB/
module load mpi
make clean; make
cd examples/thermal/rayleighBenard2d
make

thread_count=(1 2 4 8 16 32)
sim_res=(2 3 4 5 6 7)

for t in 0 1 2 3 4 5; do
    OMP_NUM_THREADS=${thread_count[$t]} salloc ./rayleighBenard2d -r ${sim_res[$t]} &> weak_ray_r=$t.txt
done