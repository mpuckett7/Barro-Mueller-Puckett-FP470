# DO NOT RUN DIRECTLY

#!/bin/bash
#SBATCH –job-name=bifurcation3d_unchanged-MPI_NUM_TASKS
#SBATCH --output=bifurcation3d_unchanged-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

module load mpi
salloc -Qn MPI_NUM_TASKS mpirun ./bifurcation3d