#!/bin/bash
# #SBATCH --job-name=BIFUR_OMP_TESTING
# #SBATCH --nodes=1

cd ../
cd OpenLB/
module load mpi
make clean; make
cd examples/particles/bifurcation3d/eulerLagrange
make

for t in 2 5 10 20; do
    OMP_NUM_THREADS=8 salloc ./bifurcation3d $t &> weak_bifur_r=$t.txt
done