#!/bin/bash 
#SBATCH --time=00:40:00	 # max 40 min each job
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12 
#SBATCH --mem=64G

#SBATCH -J ORAS_proj
#SBATCH -e /storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/ORAS/slurm_log/%J-%j-err.log
#SBATCH -o /storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/ORAS/slurm_log/%J-%j-out.log
#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr
#SBATCH --mail-type=FAIL

# set directories
JOB_DIR="/storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/ORAS"
IMG_DIR="/storage/simple/projects/t_cesab/brunno"

# load module
module purge
module load singularity/3.5

singularity exec $IMG_DIR/brunnospatial.sif Rscript $JOB_DIR/2_unzip_and_save.R
