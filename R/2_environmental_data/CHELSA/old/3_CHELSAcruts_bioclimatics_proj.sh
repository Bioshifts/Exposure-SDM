#!/bin/bash 
#SBATCH --time=00:40:00	 # max 40 min each job
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 28
#SBATCH --mem=64G

#SBATCH -J cruts_proj
#SBATCH -e /storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/CHELSAcruts/slurm_log/%J-%j-err.log
#SBATCH -o /storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/CHELSAcruts/slurm_log/%J-%j-out.log
#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr
#SBATCH --mail-type=FAIL

# get i
i=$1

# set directories
JOB_DIR="/storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/CHELSAcruts"
IMG_DIR="/storage/simple/projects/t_cesab/brunno"

# load module
module purge
module load singularity/3.5

singularity exec $IMG_DIR/brunnospatial.sif Rscript $JOB_DIR/3_CHELSAcruts_bioclimatics_proj.R $i

