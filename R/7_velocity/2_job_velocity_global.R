
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
source("R/my_functions.R")

# create dir for log files
logdir <- here::here(velocity_script_dir,"slurm-log")
if(!dir.exists(logdir)){
    dir.create(logdir,recursive = TRUE)
}

# create dir for job files
jobdir <- here::here(velocity_script_dir,"job")
if(!dir.exists(jobdir)){
    dir.create(jobdir,recursive = TRUE)
}

Rscript_file = here::here(velocity_script_dir,"2_get_velocity_global.R")


# jobs dataframe
jobs_data <- data.frame(ECO = c("Ter","Ter","Mar"),
                        velocity_variables = c("mat","map","mean"))

########################
# submit jobs
N_Nodes = 1
tasks_per_core = 1

for(i in 1:nrow(jobs_data)){
    
    # velocity variables
    velocity_variable <- jobs_data$velocity_variables[i]
    
    # ECO
    ECO <- jobs_data$ECO[i]
    
    # args
    args = c(ECO, velocity_variable)
    
    job_name <- paste("gVel",paste(args,collapse ="_"), sep = "_")
    
    if(ECO == "Mar"){ # for the Marine use this
        cores = 1
        time = "24:00:00"
        memory = "50G"
        partition = "normal-amd"
    }
    if(ECO == "Ter"){ # for the Terrestrial use this (bigger jobs) 
        cores = 5 # reduce N cores because of out-of-memory issue
        time = "1-24:00:00"
        memory = "500G"
        partition = "bigmem-amd"
    }
    
    slurm_job_singularity(jobdir = jobdir,
                          logdir = logdir, 
                          sptogo = job_name, 
                          args = args,
                          N_Nodes = N_Nodes, 
                          tasks_per_core = tasks_per_core, 
                          cores = cores, 
                          time = time, 
                          memory = memory, 
                          partition = partition, 
                          singularity_image = singularity_image, 
                          Rscript_file = Rscript_file)
    
}

