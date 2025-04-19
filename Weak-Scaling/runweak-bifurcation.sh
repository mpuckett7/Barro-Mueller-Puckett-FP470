#!/bin/bash
# #SBATCH --job-name=BIFUR_OMP_TESTING
# #SBATCH --nodes=1

cd ../
cd OpenLB/
module load mpi
make clean; make
cd examples/particles/bifurcation3d/eulerLagrange
make


thread_count=(1 2 4 8 16 32)
sim_res=(10 11 12 14 16 19)

for t in 0 1 2 3 4 5; do
    OMP_NUM_THREADS=${thread_count[$t]} salloc ./bifurcation3d ${sim_res[$t]} &> weak_bifur_r=$t.txt
done