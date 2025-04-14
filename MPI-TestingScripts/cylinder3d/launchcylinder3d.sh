#!/bin/bash
# Run to submit all jobs for cylinder3d both changed and unchanged tests using MPI.
# Make sure the MPI config file provided is the one that is named config and in the outer directories
# That would be openLB_original and OpenLB

cd ../..
cd openLB_original/examples/laminar/cylinder3d
make clean
make

for n in 1 2 4 8 16 32 64; do
    sed -e "s/MPI_NUM_TASKS/$n/g" runcylinder3d_unchanged.sh | sbatch
done

cd -
cd OpenLB/examples/laminar/cylinder3d
make clean
make

for n in 1 2 4 8 16 32 64; do
    sed -e "s/MPI_NUM_TASKS/$n/g" runcylinder3d_changed.sh | sbatch
done