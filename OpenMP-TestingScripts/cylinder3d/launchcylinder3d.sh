#!/bin/bash
# Run to submit all jobs for cylinder3d both changed and unchanged tests using OpenMP.
# Make sure the OpenMP config file provided is the one that is named config and in the outer directories
# That would be openLB_original and OpenLB


for n in 1 8 16 32 64 128; do
    sed -e "s/OMP_NUM_THREADS/$n/g" run.sh | sbatch
done