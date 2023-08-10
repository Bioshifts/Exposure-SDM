#!/bin/bash 
#SBATCH --time=01:00:00	 # max 1hr each job
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20 # for mclapply functions
#SBATCH --mem=64G

#SBATCH -J _proj
#SBATCH -e /storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/CHELSA/slurm-log-proj/%J-%j-err.log
#SBATCH -o /storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/CHELSA/slurm-log-proj/%J-%j-out.log
#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr
#SBATCH --mail-type=FAIL

# get i
i=$1

# set directories
JOB_DIR="/storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/CHELSA"
IMG_DIR="/storage/simple/projects/t_cesab/brunno"

# load module
module purge
module load singularity/3.5

singularity exec $IMG_DIR/brunnospatial.sif Rscript $JOB_DIR/2_CHELSA_bioclimatics_proj.R $i

