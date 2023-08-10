
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
    sps.dir <- here::here(work_dir,"Data/GBIF_data")
    script.dir <- here::here(work_dir,"R/4_create_PA")
}

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

Rscript_file = here::here(script.dir,"get_sdm_data_sps.R")

########################

# sp list
all_sps <- list.files(here::here(sps.dir), pattern = '.qs')
all_sps <- gsub("_"," ",all_sps)
all_sps <- gsub(".qs","",all_sps)

N_OCC <- read.csv("Data/n_occ.csv")

terrestrials <- N_OCC$scientificName[which(N_OCC$ECO == "T")]
marines <- N_OCC$scientificName[which(N_OCC$ECO == "M")]

mar_sps <- all_sps[which(all_sps %in% marines)]
mar_sps <- data.frame(sps = mar_sps, realm = "Mar")

ter_sps <- all_sps[which(all_sps %in% terrestrials)]
ter_sps <- data.frame(sps = ter_sps, realm = "Ter")


# pipe line of species: first marines (because they run faster due to coarser resolution data) then terrestrials
# all_sps <- rbind(mar_sps, ter_sps)
all_sps <- ter_sps # running only for terrestrials

########################
# submit jobs

N_jobs_at_a_time = 50

N_Nodes = 1
tasks_per_core = 1
cores = 28
time = "1:30:00"
memory = "64G"

# Check if file exists
check_if_file_exists <- TRUE

for(i in 1:nrow(all_sps)){
    
    sptogo <- all_sps$sps[i]
    sptogo <- gsub(" ","_",sptogo)
    realm <- all_sps$realm[i]
    
    args = paste(sptogo, realm)
    
    ########################
    # Check if file exists
    output_dir <- here::here(work_dir,"Data/Env_data",realm)
    
    file.test <- here::here(output_dir,paste0(gsub(" ","_",sptogo),".qs"))
    
    if(check_if_file_exists){
        RUN <- !file.exists(file.test)
    } else {
        RUN <- TRUE
    }
    if(RUN){
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
        cat("#SBATCH --mail-type=ALL\n")
        cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
        
        cat("IMG_DIR='/storage/simple/projects/t_cesab/brunno'\n")
        
        cat("module purge\n")
        cat("module load singularity/3.5\n")
        
        cat("singularity exec --disable-cache $IMG_DIR/brunnospatial.sif Rscript",Rscript_file, args,"\n", sep=" ")
        
        # Close the sink!
        sink()
        
        # Submit to run on cluster
        system(paste("sbatch", here::here(jobdir, paste0(sptogo,'.sh'))))
        
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


# check if I got env data for all species







