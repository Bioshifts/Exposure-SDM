#!/bin/bash 
#SBATCH --time=01:00:00	 # max 1hr each job
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10 # for mclapply functions
#SBATCH --mem=100G

#SBATCH -J velocity
#SBATCH -e /storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/7_velocity/%J-%j-err.log
#SBATCH -o /storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/7_velocity/%J-%j-out.log
#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr
#SBATCH --mail-type=FAIL

# set directories
JOB_DIR="/storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/7_velocity"
IMG_DIR="/storage/simple/projects/t_cesab/brunno"

# load module
module purge
module load singularity/3.5

singularity exec $IMG_DIR/brunnospatial.sif Rscript $JOB_DIR/1_get_velocity.R

