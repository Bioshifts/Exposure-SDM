
### P.S.: Cannot run this inside singularity container
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)

########################
# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    work_dir <- getwd()
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
    work_dir <- getwd()
}

# source settings
source("R/settings.R")

script.dir <- here::here(work_dir,"R/7_velocity")

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

Rscript_file = here::here(script.dir,"3_global_velocity_maps.R")


# jobs dataframe
jobs_data <- data.frame(ECO = c("Ter","Ter","Mar"),
                        velocity_variables = c("mat","map","mean"))

########################
# submit jobs

N_Nodes = 1
tasks_per_core = 1
cores = 1 # this is not a parallel job


for(i in 1:nrow(jobs_data)){
    
    
    # velocity variables
    velocity_variable <- jobs_data$velocity_variables[i]
    
    # ECO
    ECO <- jobs_data$ECO[i]
    
    # args
    args = c(ECO, velocity_variable)
    
    job_name <- paste("gVel",paste(args,collapse ="_"), sep = "_")
    
    # Start writing to this file
    sink(here::here(jobdir,paste0(job_name,'.sh')))
    
    # the basic job submission script is a bash script
    cat("#!/bin/bash\n")
    
    cat("#SBATCH -N",N_Nodes,"\n")
    cat("#SBATCH -n",tasks_per_core,"\n")
    cat("#SBATCH -c",cores,"\n")
    cat("#SBATCH --mem=",ifelse(ECO == "Ter", "400G", "120G"),"\n", sep="")
    
    ifelse(ECO=="Ter",
        cat("#SBATCH --partition=bigmem\n"),
        cat("#SBATCH --partition=normal\n"))
    
    cat("#SBATCH --job-name=",job_name,"\n", sep="")
    cat("#SBATCH --output=",here::here(logdir,paste0(job_name,".out")),"\n", sep="")
    cat("#SBATCH --error=",here::here(logdir,paste0(job_name,".err")),"\n", sep="")
    cat("#SBATCH --time=",ifelse(ECO == "Ter", "2-20:00:00", "24:00:00"),"\n", sep="")
    cat("#SBATCH --mail-type=ALL\n")
    cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
    
    cat(paste0("IMG_DIR='",singularity_image,"'\n"))
    
    cat("module purge\n")
    cat("module load singularity\n")
    
    cat("singularity exec --disable-cache $IMG_DIR Rscript",Rscript_file, args,"\n", sep=" ")
    
    # Close the sink!
    sink()
    
    # Submit to run on cluster
    system(paste("sbatch", here::here(jobdir, paste0(job_name,'.sh'))))
    
}

