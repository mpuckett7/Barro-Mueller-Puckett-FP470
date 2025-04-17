#!/bin/bash
# Run to submit all jobs for rayleighBenard2d both changed and unchanged tests using OpenMP.
# Make sure the OpenMP config file provided is the one that is named config and in the outer directories
# That would be openLB_original and OpenLB


for n in 1 do
    sed -e "s/N/$n/g" runrayleigh_changed.sh | sbatch
done
