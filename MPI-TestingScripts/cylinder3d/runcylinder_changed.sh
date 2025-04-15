#!/bin/bash
#SBATCH --job-name=clyinder3d_changed-MPI_NUM_TASKS
#SBATCH --output=cylinder3d_changed-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

cd ../../
cd OpenLB/examples/laminar/cylinder3d

make

module load mpi
salloc -Qn 4 MPI_NUM_TASKS mpirun ./cylinder3d