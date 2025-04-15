#!/bin/bash
#SBATCH --job-name=clyinder3d_changed-MPI_NUM_TASKS
#SBATCH --output=cylinder3d_changed-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../../
cd OpenLB/examples/laminar/cylinder3d
make clean
make

module load mpi
mpirun ./cylinder3d