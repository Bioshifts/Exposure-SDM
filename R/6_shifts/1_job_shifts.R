
### Cannot run this inside singularity
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()

library(tictoc)
library(here)

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

########################
# set computer
computer = "matrics"

realms = c("Mar","Ter")

for(j in 1:length(realms)){
    
    realm <- realms[j]
    
    if(computer == "muse"){
        setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    }
    if(computer == "matrics"){
        setwd("/users/boliveira/Exposure-SDM")
    }
    
    work_dir <- getwd()
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
    if(!dir.exists(shift_dir(realm))){
        dir.create(shift_dir(realm),recursive = TRUE)
    }
    
    ########################
    # sp list I have sdms
    all_sps <- list.dirs(here::here(sdm_dir(realm)), recursive = FALSE, full.names = FALSE)
    all_sps <- gsub("_"," ",all_sps)
    length(all_sps)
    
    # get info from v1
    biov1 <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v1), header = T)
    # only LAT
    biov1$Type[which(biov1$Type=="HOR")] <- "LAT"
    biov1 <- biov1[which(biov1$Type=="LAT"),]
    biov1 <- biov1[which(biov1$sp_name_std_v1 %in%  gsub(" ","_",all_sps)),]
    
    # sp list - Species with full sdms
    I_have_sdms <- c()
    for(i in 1:length(all_sps)) { cat(i, "from", length(all_sps),"\r")
        sp_i <- gsub(" ","_",all_sps[i])
        sp_i_realm <- realm
        # get studies for species i
        ID_i <- unique(biov1$ID[which(biov1$sp_name_std_v1 == sp_i)])
        # for each ID_i, look if has projections for all years
        I_have_sdms[i] <- all(sapply(ID_i, function(x){
            bio_i <- biov1[which(biov1$ID == x),]
            bio_i <- biov1[which(biov1$sp_name_std_v1 == sp_i),]
            years_ID_i <- round(bio_i$START[1],0):round(bio_i$END[1],0)
            # check
            # Focus on SA for now
            sdms_sp_i <- list.files(paste(sdm_dir(sp_i_realm),sp_i,gsub("_",".",sp_i),sep = "/"),pattern = "SA")
            # get ensemble models
            sdms_sp_i_ens <- sdms_sp_i[grep(" ens",sdms_sp_i)]
            # get single algorithms
            # sdms_sp_i <- sdms_sp_i[-grep(" ens",sdms_sp_i)]
            # check if all years exist
            sdms_sp_i_ens <- all(sapply(years_ID_i, function(x){any(grepl(x,sdms_sp_i_ens))}))
            # sdms_sp_i <- all(sapply(years_ID_i, function(x){any(grepl(x,sdms_sp_i))}))
            # all(sdms_sp_i,sdms_sp_i_ens)
            sdms_sp_i_ens
        }))
    }
    I_have_sdms <- all_sps[which(I_have_sdms)]
    length(I_have_sdms)
    
    # select from the list of species I have sdms the ones which we still did not have shift calculated
    I_have_shift <- list.files(shift_dir(realm))
    I_have_shift <- I_have_shift[grep(" ens SA",I_have_shift)]
    I_have_shift <- gsub(" shift ens SA.csv","",I_have_shift)
    length(I_have_shift)
    
    # missing
    missing <- I_have_sdms[!I_have_sdms %in% gsub("_"," ",I_have_shift)]
    length(missing)
    
    # species to go
    all_sps <- missing
    
    ########################
    # submit jobs
    
    N_jobs_at_a_time = 100
    
    N_Nodes = 1
    tasks_per_core = 1
    cores = 1
    time = "40:00"
    memory = "16G"
    
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
        if(computer == "matrics"){
            cat("#SBATCH --partition=normal-amd\n")
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

