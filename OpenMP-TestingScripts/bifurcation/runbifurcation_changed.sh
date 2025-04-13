# DO NOT RUN DIRECTLY

#!/bin/bash
#SBATCH –job-name=bifurcation3d_changed-OMP_NUM_THREADS
#SBATCH --output=bifurcation3d_changed-OMP_NUM_THREADS.txt
#SBATCH --ntasks=OMP_NUM_THREADS

OMP_NUM_THREADS=OMP_NUM_THREADS ./bifurcation3d