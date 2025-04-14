#!/bin/bash
# Run to submit all jobs for EXAMPLE-NAME both changed and unchanged tests using OpenMP.
# Make sure the OpenMP config file provided is the one that is named config and in the outer directories
# That would be openLB_original and OpenLB

cd ../..
cd openLB_original/examples/INSERT_PATH
make clean
make

for n in 1 2 4 8 16 32 64; do
    sed -e "s/OMP_NUM_THREADS/$n/g" runEXAMPLE-NAME_unchanged.sh | sbatch
done

cd -
cd OpenLB/examples/INSERT_PATH
make clean
make

for n in 1 2 4 8 16 32 64; do
    sed -e "s/OMP_NUM_THREADS/$n/g" runEXAMPLE-NAME_changed.sh | sbatch
done
