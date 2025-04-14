# DO NOT RUN DIRECTLY

#!/bin/bash
#SBATCH –job-name=clyinder3d_changed-MPI_NUM_TASKS
#SBATCH --output=cylinder3d_changed-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

module load mpi
salloc -Qn MPI_NUM_TASKS mpirun ./cylinder3d