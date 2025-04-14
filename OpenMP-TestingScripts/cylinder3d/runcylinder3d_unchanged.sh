# DO NOT RUN DIRECTLY

#!/bin/bash
#SBATCH –job-name=cylinder3d_unchanged-OMP_NUM_THREADS
#SBATCH --output=cylinder3d_unchanged-OMP_NUM_THREADS.txt
#SBATCH --ntasks=OMP_NUM_THREADS

OMP_NUM_THREADS=OMP_NUM_THREADS ./cylinder3d