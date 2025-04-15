#!/bin/bash
#SBATCH --job-name=cylinder3d_unchanged-MPI_NUM_TASKS
#SBATCH --output=cylinder3d_unchanged-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../..
cd openLB_original/examples/laminar/cylinder3d
make clean
make

module load mpi
salloc -Qn MPI_NUM_TASKS mpirun ./cylinder3d