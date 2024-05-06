
### Cannot run this inside singularity
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()

library(tictoc)
library(here)
library(dplyr)

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
    length(all_sps)
    
    test <- sapply(all_sps, function(x){
        file.exists(here::here(sdm_dir(realm),x,paste(x,"shift_info.csv")))
    })
    all_sps <- all_sps[test]
    length(all_sps)
    
    # Species with sdms projections for all years
    I_have_sdms <- data.frame()
    for(i in 1:length(all_sps)) { cat(i, "from", length(all_sps),"\r")
        
        sp_i <- all_sps[i]
        sp_i_realm <- realm
        
        # get studies for species i
        shift_info <- read.csv(here::here(sdm_dir(realm),sp_i,paste(sp_i,"shift_info.csv")))
        
        ID_i <- shift_info$ID
        # for each ID_i, look if has projections for all years
        tmp <- data.frame(sps=sp_i,
                          realm=sp_i_realm,
                          ID=ID_i)
        
        tmp$I_have_sdms <- sapply(ID_i, function(x){
            shift_info_i <- shift_info[which(shift_info$ID==x),]
            years_ID_i <- round(shift_info_i$Start,0):round(shift_info_i$End,0)
            # check
            # Focus on SA for now
            sdms_sp_i <- list.files(paste(sdm_dir(sp_i_realm),sp_i,gsub("_",".",sp_i),sep = "/"),pattern = "SA")
            # get ensemble models
            sdms_sp_i_ens <- sdms_sp_i[grep(" ens",sdms_sp_i)]
            # get ID_i models
            sdms_sp_i_ens <- sdms_sp_i_ens[grep(x,sdms_sp_i_ens)]
            # check if all years exist
            sdms_sp_i_ens <- all(sapply(years_ID_i, function(x){any(grepl(x,sdms_sp_i_ens))}))
            sdms_sp_i_ens
        })
        I_have_sdms <- rbind(I_have_sdms,tmp)
    }
    nrow(I_have_sdms) # these are all possible sdms (species + ID)
    I_have_sdms <- I_have_sdms[which(I_have_sdms$I_have_sdms),]
    nrow(I_have_sdms) # these are the ones with projections for all years
    
    # sps I have shifts calculated
    I_have_shift <- list.files(shift_dir(realm))
    I_have_shift <- lapply(I_have_shift, function(x){
        tmp <- strsplit(x,"_")[[1]]
        data.frame(sps=paste(tmp[1],tmp[2],sep="_"),realm=tmp[7],ID=paste(tmp[3],tmp[4],sep="_"))
    })
    I_have_shift <- do.call(rbind,I_have_shift)
    
    nrow(I_have_shift)
    
    # select from the list of species I have sdms the ones which we still did not have shift calculated
    # missing
    missing <- I_have_sdms[which(!paste(I_have_sdms$sps,I_have_sdms$ID) %in% paste(I_have_shift$sps,I_have_shift$ID)),]
    nrow(missing)
    
    missing <- missing %>% distinct(sps, ID, .keep_all = TRUE)
    
    if(nrow(missing)>0){
        ########################
        # submit jobs
        
        N_jobs_at_a_time = 50
        N_Nodes = 1
        tasks_per_core = 1
        
        cat("Running for", nrow(missing), realm, "species\n")
        
        for(i in 1:nrow(missing)){ 
            
            sptogo <- missing$sps[i]
            realmtogo <- missing$realm[i]
            IDtogo <- missing$ID[i]
            
            args <- paste(sptogo,realmtogo,IDtogo)
            
            shifttogo = gsub(" ","_",args)
            
            if(realm == "Mar"){ # for the Marine use this
                cores = 20
                time = "5:00:00"
                memory = "16G"
                partition = "normal"
            }
            if(realm == "Ter"){ # for the Terrestrial use this (bigger jobs) 
                cores = 5 # reduce N cores because of out-of-memory issue
                time = "20:00:00"
                memory = "100G"
                partition = "bigmem"
            }
            
            slurm_job_singularity(jobdir = jobdir,
                                  logdir = logdir, 
                                  sptogo = shifttogo, 
                                  args = args,
                                  N_Nodes = N_Nodes, 
                                  tasks_per_core = tasks_per_core, 
                                  cores = cores, 
                                  time = time, 
                                  memory = memory, 
                                  partition = partition, 
                                  singularity_image = singularity_image, 
                                  Rscript_file = Rscript_file)
            
            # check how many jobs in progress
            tmp <- system("squeue -u $USER",intern = T)
            
            while(length(tmp)>N_jobs_at_a_time){
                
                Sys.sleep(10)
                tmp <- system("squeue -u $USER",intern = T)
                
            }
        }
    }
}

