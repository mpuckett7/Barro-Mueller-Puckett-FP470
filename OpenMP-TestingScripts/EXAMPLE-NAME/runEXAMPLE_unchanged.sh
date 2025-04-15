#!/bin/bash
#SBATCH --job-name=EXAMPLE-NAME_unchanged-OMP_NUM_THREADS
#SBATCH --output=EXAMPLE-NAME_unchanged-OMP_NUM_THREADS.txt
#SBATCH --ntasks=OMP_NUM_THREADS

cd ../..
cd openLB_original/examples/INSERT_PATH
make clean
make

OMP_NUM_THREADS=OMP_NUM_THREADS ./EXENAME