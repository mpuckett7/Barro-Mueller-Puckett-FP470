#!/bin/bash
# Run to submit all jobs for bifurcation both changed and unchanged tests using OpenMP.
# Make sure the OpenMP config file provided is the one that is named config and in the outer directories
# That would be openLB_original and OpenLB

for n in 1 2 4 8 16 32 64; do
    sed -e "s/OMP_NUM_THREADS/$n/g" runbifurcation_changed.sh | sbatch
done
