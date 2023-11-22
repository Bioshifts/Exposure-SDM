# load the list of species created with 2_check_if_sdms_ran_test2.R


missing_sps <- read.csv("missing_sps.csv")

dim(missing_sps)

########################
# submit jobs

computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    
    work_dir <- getwd()
    
    script.dir <- here::here(work_dir,"R/5_run_sdms")
}

Rscript_file = here::here(script.dir,"1_run_sdm.R")

# create dir for log files
logdir <- here::here(script.dir,"slurm-log")
if(!dir.exists(logdir)){
    dir.create(logdir,recursive = TRUE)
}

# create dir for job files
jobdir <- here::here(script.dir,"job")
if(!dir.exists(jobdir)){
    dir.create(jobdir,recursive = TRUE)
}

N_jobs_at_a_time = 100

N_Nodes = 1
tasks_per_core = 1
cores = 1
time = "40:00"
memory = "64G"

# wait all jobs to finish before running this one
tmp <- system("squeue -u $USER",intern = T)

while(length(tmp)>0){
    
    Sys.sleep(60)
    tmp <- system("squeue -u $USER",intern = T)
    
}


for(i in 1:nrow(missing_sps)){
    
    sptogo <- missing_sps$sps[i]
    sptogo <- gsub(" ","_",sptogo)
    
    realm = missing_sps$realm[i]
    args = paste(sptogo, realm)
    
    # Start writing to this file
    sink(here::here(jobdir,paste0(sptogo,'.sh')))
    
    # the basic job submission script is a bash script
    cat("#!/bin/bash\n")
    
    cat("#SBATCH -N",N_Nodes,"\n")
    cat("#SBATCH -n",tasks_per_core,"\n")
    cat("#SBATCH -c",cores,"\n")
    
    cat("#SBATCH --job-name=",sptogo,"\n", sep="")
    cat("#SBATCH --output=",here::here(logdir,paste0(sptogo,".out")),"\n", sep="")
    cat("#SBATCH --error=",here::here(logdir,paste0(sptogo,".err")),"\n", sep="")
    cat("#SBATCH --time=",time,"\n", sep="")
    cat("#SBATCH --mem=",memory,"\n", sep="")
    # cat("#SBATCH --mail-type=ALL\n")
    # cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
    
    cat("IMG_DIR='/storage/simple/projects/t_cesab/brunno'\n")
    
    cat("module purge\n")
    cat("module load singularity/3.5\n")
    
    cat("singularity exec --disable-cache $IMG_DIR/brunnospatial.sif Rscript",Rscript_file, args,"\n", sep=" ")
    
    # Close the sink!
    sink()
    
    # Submit to run on cluster
    system(paste("sbatch", here::here(jobdir, paste0(sptogo,'.sh'))))
    
    
}