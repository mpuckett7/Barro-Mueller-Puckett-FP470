#!/bin/bash
# Run to submit all jobs for RayleighBrenard2d both changed and unchanged tests using MPI.
# Make sure the MPI config file provided is the one that is named config and in the outer directories
# That would be openLB_original and OpenLB

module load mpi

sed -e "s/MPI_NUM_TASKS/1/g" runrayleigh1.sh | sbatch

sed -e "s/MPI_NUM_TASKS/2/g" runrayleigh2.sh | sbatch

for n in 4 8 16 32; do
    sed -e "s/MPI_NUM_TASKS/$n/g" runrayleigh_changed.sh | sbatch
done