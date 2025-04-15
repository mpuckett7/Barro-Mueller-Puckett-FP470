#!/bin/bash
#SBATCH --job-name=dkt2d_unchanged-MPI_NUM_TASKS
#SBATCH --output=dkt2d_unchanged-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../..
cd openLB_original/examples/particles/dkt2d
make clean
make

module load mpi
salloc -Qn MPI_NUM_TASKS mpirun ./dkt2d