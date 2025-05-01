#!/bin/bash
#Run to submit all jobs for bifurcation using hybrid parallelism.

module load mpi

for n in 2 4 8; do
    sed -e "s/NUM_TASKS/$n/g" runbifurcation.sh | sbatch
done
