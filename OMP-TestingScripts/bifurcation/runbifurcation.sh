#!/bin/bash
# #SBATCH --job-name=BIFURCATION_OMP_TESTS
# #SBATCH --nodes=1

cd ../../
cd OpenLB/
module load mpi
make clean; make
cd examples/particles/bifurcation3d/eulerLagrange
make

for t in 1 2 4 8 16 32; do
    OMP_NUM_THREADS=$t salloc ./bifurcation3d &> bifurcation$t.txt
done