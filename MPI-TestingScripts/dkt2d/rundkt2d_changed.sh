# DO NOT RUN DIRECTLY

#!/bin/bash
#SBATCH –job-name=dkt2d_changed-MPI_NUM_TASKS
#SBATCH --output=dkt2d_changed-MPI_NUM_TASKS.txt
#SBATCH --ntasks=MPI_NUM_TASKS

module load mpi
salloc -Qn MPI_NUM_TASKS mpirun ./dkt2d