
### P.S.: Cannot run this inside singularity container
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    
    work_dir <- getwd()
}

# source settings
source("R/settings.R")

script.dir <- velocity_SA_scrit_dir

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

Rscript_file = here::here(script.dir,"1_get_velocity.R")

# Load study areas v3
v3_polygons <- gsub(".shp","",list.files(SA_shps_dir,pattern = ".shp"))

########################
# submit jobs

N_jobs_at_a_time = 50

N_Nodes = 1
tasks_per_core = 1
cores = 5 # this is not a parallel job
# time = "5:00:00"
# memory = "64G"

# for the big jobs
time = "2-20:00:00"
memory = "120G"


# Check if file exists
check_if_file_exists <- TRUE

for(i in 1:length(v3_polygons)){
    
    SAtogo <- v3_polygons[i]
    
    args = SAtogo
    
    ########################
    # Check if file exists
    
    file.test <- here::here(velocity_SA_dir, paste0(SAtogo,".csv"))

    if(check_if_file_exists){
        RUN <- !file.exists(file.test)
    } else {
        RUN <- FALSE
    }
    if(RUN){
        # Start writing to this file
        sink(here::here(jobdir,paste0(SAtogo,'.sh')))
        
        # the basic job submission script is a bash script
        cat("#!/bin/bash\n")
        
        cat("#SBATCH -N",N_Nodes,"\n")
        cat("#SBATCH -n",tasks_per_core,"\n")
        cat("#SBATCH -c",cores,"\n")
        cat("#SBATCH --mem=",memory,"\n", sep="")
        
        cat("#SBATCH --job-name=",SAtogo,"\n", sep="")
        cat("#SBATCH --output=",here::here(logdir,paste0(SAtogo,".out")),"\n", sep="")
        cat("#SBATCH --error=",here::here(logdir,paste0(SAtogo,".err")),"\n", sep="")
        cat("#SBATCH --time=",time,"\n", sep="")
        cat("#SBATCH --mail-type=ALL\n")
        cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
        
        cat("IMG_DIR='/storage/simple/projects/t_cesab/brunno'\n")
        
        cat("module purge\n")
        cat("module load singularity/3.5\n")
        
        cat("singularity exec --disable-cache $IMG_DIR/brunnospatial.sif Rscript",Rscript_file, args,"\n", sep=" ")
        
        # Close the sink!
        sink()
        
        # Submit to run on cluster
        system(paste("sbatch", here::here(jobdir, paste0(SAtogo,'.sh'))))
        
        # check how many jobs in progress
        n_sps_going <- system("squeue -u $USER",intern = T)
        
        # N sps going
        n_sps_going <- length(n_sps_going) -1
        while(n_sps_going>N_jobs_at_a_time){
            
            Sys.sleep(10)
            n_sps_going <- system("squeue -u $USER",intern = T)
            n_sps_going <- length(n_sps_going) -1
            
        }
    }
}


# # check if I got velocities for all SA
# # list of SA we got data
# SA_got <- list.files(velocity_SA_dir)
# SA_got <- gsub(".csv","",SA_got)
# missing_SA <- v3_polygons[!v3_polygons %in% SA_got]
# length(missing_SA)
# head(missing_SA)
# 
# v3_polygons <- missing_SA
# 
# # Check error file
# error_f <- lapply(missing_SA, function(x){
#     tmp <- read.csv(here::here(script.dir,"slurm-log",paste0(x,".err")))
#     tmp[nrow(tmp),]
# })
# error_f


