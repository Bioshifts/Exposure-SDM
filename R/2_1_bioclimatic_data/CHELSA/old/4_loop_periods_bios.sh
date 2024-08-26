#! /bin/bash
#SBATCH --time=2-1:00:00
#SBATCH --nodes=1					# Use one node
#SBATCH --ntasks=1					# N jobs running in a single node
#SBATCH --cpus-per-task=1				# N cores (only for parallelization)

set -x # tipo executavel

N=1380 # N periods to go 1380
at_a_time=50 # N jobs that run at once 

# Adjust to 1 less than the starting point
typeset -i count=0

while [[ $count -le $N ]]; do # -le = less or equal
    if [[ $(squeue -u $USER --noheader | wc -l) -ge $at_a_time ]]; then # se a quatidade de jobs for maior que at_a_time (ge = greater than or equal)
        sleep 60
        continue # reinicia o loop
    fi
    
    count+=1
    
    # submit a job here
    JOB_NAME="CHELbios_${count}"
    JOB_NAME="J${count}"
    JOB_DIR="/storage/simple/projects/t_cesab/brunno/Exposure-SDM/R/2_environmental_data/CHELSAcruts"

    sbatch --job-name=$JOB_NAME -N 1 -n 1 -c 1 $JOB_DIR/4_CHELSAcruts_bioclimatics_bios.sh $count

done

