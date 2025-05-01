#!/bin/bash
# Run to submit all jobs for bifurcation3d both changed and unchanged tests using MPI.
# Make sure the MPI config file provided is the one that is named config and in the outer directories

module load mpi

sed -e "s/MPI_NUM_TASKS/1/g" runbifurcation1.sh | sbatch

for n in 2 4 8 16 32; do
    sed -e "s/MPI_NUM_TASKS/$n/g" runbifurcation_changed.sh | sbatch
done