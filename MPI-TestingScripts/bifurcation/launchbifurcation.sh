#!/bin/bash
# Run to submit all jobs for bifurcation3d both changed and unchanged tests using MPI.
# Make sure the MPI config file provided is the one that is named config and in the outer directories
# That would be openLB_original and OpenLB

module load mpi

# for n in 1 2 4 8 16 32 64; do
#     sed -e "s/MPI_NUM_TASKS/$n/g" runbifurcation_unchanged.sh | sbatch
# done

for n in 1 2 4 8 16 32 64; do
    sed -e "s/MPI_NUM_TASKS/$n/g" runbifurcation_changed.sh | sbatch
done