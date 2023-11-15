
### Cannot run this inside singularity
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)

########################
# set computer
computer = "matrics"

for(j in c("Mar","Ter")){
    
    
    if(computer == "muse"){
        setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
        work_dir <- getwd()
    }
    if(computer == "matrics"){
        setwd("/users/boliveira/Exposure-SDM")
        work_dir <- getwd()
    }
    
    script.dir <- here::here(work_dir,"R/6_shifts")
    
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
    
    Rscript_file = here::here(script.dir,"1_get_shifts.R")
    
    ########################
    # load functions
    source("R/my_functions.R")
    source("R/range_shift_functions.R")
    # source settings
    source("R/settings.R")
    
    # create dirs
    realm <- j
    if(!dir.exists(shift_dir(realm))){
        dir.create(shift_dir(realm),recursive = TRUE)
    }
    
    ########################
    # sp list
    all_sps <- list.dirs(here::here(sdm_dir(realm)), recursive = FALSE, full.names = FALSE)
    all_sps <- gsub("_"," ",all_sps)
    
    length(all_sps)
    
    # select only v1 species for now
    biov1 <- read.csv("Data/Bioshifts/biov1_fixednames.csv", header = T)
    biov1 <- biov1[biov1$Type=="LAT",]
    
    biov1_sps <- unique(gsub("_"," ",biov1$sp_name_std_v1))
    
    all_sps <- all_sps[all_sps %in% biov1_sps]
    
    ########################
    # submit jobs
    
    N_jobs_at_a_time = 100
    
    N_Nodes = 1
    tasks_per_core = 1
    cores = 1
    time = "40:00"
    memory = "16G"
    
    cat("Running for", length(all_sps), realm, "species\n")
    
    # all_sps = missing_sps
    
    for(i in 1:length(all_sps)){
        
        sptogo <- all_sps[i]
        sptogo <- gsub(" ","_",sptogo)
        
        args = paste(sptogo, realm)
        
        # Start writing to this file
        sink(here::here(jobdir,paste0(sptogo,'.sh')))
        
        # the basic job submission script is a bash script
        cat("#!/bin/bash\n")
        
        cat("#SBATCH -N",N_Nodes,"\n")
        cat("#SBATCH -n",tasks_per_core,"\n")
        cat("#SBATCH -c",cores,"\n")
        if(computer == "matrics"){
            cat("#SBATCH --partition=normal\n")
        }
        
        cat("#SBATCH --job-name=",sptogo,"\n", sep="")
        cat("#SBATCH --output=",here::here(logdir,paste0(sptogo,".out")),"\n", sep="")
        cat("#SBATCH --error=",here::here(logdir,paste0(sptogo,".err")),"\n", sep="")
        cat("#SBATCH --time=",time,"\n", sep="")
        cat("#SBATCH --mem=",memory,"\n", sep="")
        # cat("#SBATCH --mail-type=ALL\n")
        # cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
        
        if(computer == "muse"){
            cat("IMG_DIR='/storage/simple/projects/t_cesab/brunno'\n")
            
            cat("module purge\n")
            cat("module load singularity/3.5\n")
            
            cat("singularity exec --disable-cache $IMG_DIR/brunnospatial.sif Rscript",Rscript_file, args,"\n", sep=" ")
        }
        if(computer == "matrics"){
            cat(paste0("IMG_DIR='",singularity_image,"'\n"))
            
            cat("module purge\n")
            cat("module load singularity\n")
            
            cat("singularity exec --disable-cache $IMG_DIR Rscript",Rscript_file, args,"\n", sep=" ")
        }
        
        # Close the sink!
        sink()
        
        # Submit to run on cluster
        system(paste("sbatch", here::here(jobdir, paste0(sptogo,'.sh'))))
        
        # check how many jobs in progress
        tmp <- system("squeue -u $USER",intern = T)
        
        while(length(tmp)>N_jobs_at_a_time){
            
            Sys.sleep(10)
            tmp <- system("squeue -u $USER",intern = T)
            
        }
        
    }
    
}


# check if everything run
test <- list.files(shiftplot_dir(realm))
test <- sapply(test, function(x) strsplit(x," ")[[1]][1])
test <- gsub("_"," ",test)
# missing
missing_sps <- all_sps[!all_sps %in% test]
missing_sps