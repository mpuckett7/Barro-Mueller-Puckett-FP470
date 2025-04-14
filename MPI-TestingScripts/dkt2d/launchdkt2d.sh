#!/bin/bash
# Run to submit all jobs for dkt2d both changed and unchanged tests using MPI.
# Make sure the MPI config file provided is the one that is named config and in the outer directories
# That would be openLB_original and OpenLB

cd ../..
cd openLB_original/examples/particles/dkt2d
make clean
make

for n in 1 2 4 8 16 32 64; do
    sed -e "s/MPI_NUM_TASKS/$n/g" rundkt2d_unchanged.sh | sbatch
done

cd -
cd OpenLB/examples/particles/dkt2d
make clean
make

for n in 1 2 4 8 16 32 64; do
    sed -e "s/MPI_NUM_TASKS/$n/g" rundkt2d_changed.sh | sbatch
done