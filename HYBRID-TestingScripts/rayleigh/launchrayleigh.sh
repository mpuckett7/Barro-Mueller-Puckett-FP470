#!/bin/bash
#Run to submit all jobs for rayleighbenard2d with hybrid paralleism

module load mpi

for n in 2 4 8; do
    sed -e "s/NUM_TASKS/$n/g" runrayleigh.sh | sbatch
done

sed -e "s/NUM_TASKS/16/g" runrayleigh2-16.sh | sbatch
