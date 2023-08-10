
### Cannot run this inside singularity
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)

########################
# set computer
computer = "muse"

for(j in c("Mar","Ter")){
    
    
    if(computer == "muse"){
        setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
        
        realm = j
        
        work_dir <- getwd()
        
        shift_dir <- here::here(work_dir,"Data/SHIFT",realm)
        scratch_dir <- "/lustre/oliveirab"
        sdm_dir <- here::here(scratch_dir,"SDMs",realm)
        script.dir <- here::here(work_dir,"R/6_shifts")
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
    
    Rscript_file = here::here(script.dir,"1_get_shifts.R")
    
    ########################
    # load functions
    source("R/my_functions.R")
    source("R/range_shift_functions.R")
    # source settings
    source("R/settings.R")
    
    # create dirs
    if(!dir.exists(shift_dir)){
        dir.create(shift_dir,recursive = TRUE)
    }
    
    ########################
    # sp list
    all_sps <- list.dirs(here::here(sdm_dir), recursive = FALSE, full.names = FALSE)
    all_sps <- gsub("_"," ",all_sps)
    
    length(all_sps)
    
    # select species SDMs worked out
    check <- sapply(1:length(all_sps), function(i) {
        
        # check if dir exists
        tmp = dir.exists(here::here(sdm_dir,
                                    gsub(" ","_",all_sps[i]),
                                    gsub(" ",".",all_sps[i])))
        
        if(!tmp){
            FALSE
        } else {
            # check if there a file with ensemble model outputs
            tmp = list.files(here::here(sdm_dir,
                                        gsub(" ","_",all_sps[i]),
                                        gsub(" ",".",all_sps[i])))
            tmp = any(grepl("SA start ens",tmp))
        }
        
        return(tmp)
        
    })
    table(check)
    all_sps <- all_sps[check]
    
    length(all_sps)
    
    # select only v1 species for now
    biov1 <- read.csv("Data/Bioshifts/biov1_fixednames.csv", header = T)
    biov1 <- biov1[biov1$Type=="LAT",]
    
    biov1_sps <- unique(gsub("_"," ",biov1$sp_name_std_v1))
    
    all_sps <- all_sps[all_sps %in% biov1_sps]
    
    length(all_sps)
    
    # Run for species from which shift was not calculated yet or with errors in shift
    # species shift were calculated
    tmp <- list.files(shift_dir, recursive = T, pattern = " ens shift SA.csv")
    tmp <- gsub(" ens shift SA.csv","",tmp)
    tmp <- gsub("_"," ",tmp)
    # missing 1
    miss1 <- all_sps[!all_sps %in% tmp]
    # Find the ones with errors
    tmp <- list.files(shift_dir, recursive = T, pattern = " ens shift SA.csv")
    tmp <- gsub(" ens shift SA.csv","",tmp)
    tmp <- gsub("_"," ",tmp)
    tmp2 <- sapply(tmp, function(x){
        try({
            test <- read.csv(here::here(shift_dir,paste0(gsub(" ","_",x)," ens shift SA.csv")))
        is.na(test[1,1])
        },silent = TRUE)
    })
    # missing 2
    miss2 <- tmp[tmp2]
    # group
    miss <- na.omit(unique(c(miss1,miss2)))
    
    all_sps <- miss
    
    length(all_sps)
    
    
    ########################
    # submit jobs
    
    N_jobs_at_a_time = 100
    
    N_Nodes = 1
    tasks_per_core = 1
    cores = 1
    time = "40:00"
    memory = "64G"
    
    cat("Running for", length(all_sps), realm, "species\n")
    
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
        
        # check how many jobs in progress
        tmp <- system("squeue -u $USER",intern = T)
        
        while(length(tmp)>N_jobs_at_a_time){
            
            Sys.sleep(10)
            tmp <- system("squeue -u $USER",intern = T)
            
        }
        
    }
    
}

