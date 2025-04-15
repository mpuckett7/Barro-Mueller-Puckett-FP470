# DO NOT RUN DIRECTLY

#!/bin/bash
#SBATCH –job-name=bifurcation3d_changed-MPI_NUM_TASKS
#SBATCH --output=bifurcation3d_changed-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

module load mpi
mpirun ./bifurcation3d