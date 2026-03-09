
### Cannot run this inside singularity
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()

library(tictoc)
library(here)
library(dplyr)
library(pbapply)

# set computer
computer = "matrics"

if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}
work_dir <- getwd()

# source settings
source("R/settings.R")
source("R/my_functions.R")

script_dir <- here::here(work_dir,"R/0_supplementary")

# create dir for log files
logdir <- here::here(script_dir,"slurm-log")
if(!dir.exists(logdir)){
    dir.create(logdir,recursive = TRUE)
}

# create dir for job files
jobdir <- here::here(script_dir,"job")
if(!dir.exists(jobdir)){
    dir.create(jobdir,recursive = TRUE)
}

Rscript_file = here::here(script_dir,"3_run_sdm.R")

########################
# sp I have environmental data
all_sps_mar <- list.files(here(env_data_dir("Mar"),"Supp"), 
                          pattern = '.qs', 
                          recursive = TRUE)
all_sps_mar <- gsub(".qs","",all_sps_mar)
all_sps_mar <- sapply(all_sps_mar, function(x){
    tmp <- strsplit(x,"_")[[1]]
    return(paste(tmp[1:2],collapse = " "))
})

all_sps_ter <- list.files(here(env_data_dir("Ter"),"Supp"), 
                          pattern = '.qs', 
                          recursive = TRUE)
all_sps_ter <- gsub(".qs","",all_sps_ter)
all_sps_ter <- sapply(all_sps_ter, function(x){
    tmp <- strsplit(x,"_")[[1]]
    return(paste(tmp[1:2],collapse = " "))
})

all_sps <- unique(c(all_sps_ter,all_sps_mar))

length(all_sps) 

#####################
# Separate terrestrials from marines
bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_sdms_selection(bioshifts)
bioshifts <- filter(bioshifts, gsub("_"," ",sp_name_std) %in% all_sps)
# head(bioshifts)

terrestrials <- bioshifts$sp_name_std[which(bioshifts$Eco == "Ter")]
marines <- bioshifts$sp_name_std[which(bioshifts$Eco == "Mar")]

ter_sps <- all_sps[which(all_sps %in% gsub("_"," ",terrestrials))]
ter_sps <- data.frame(sps = ter_sps, realm = "Ter")
nrow(ter_sps) 

mar_sps <- all_sps[which(all_sps %in% gsub("_"," ",marines))]
mar_sps <- data.frame(sps = mar_sps, realm = "Mar")
nrow(mar_sps) 

# pipe line of species: first marines (because they run faster due to coarser resolution data)
all_sps <- rbind(mar_sps, ter_sps)
all_sps$sps <- gsub(" ","_",all_sps$sps)
nrow(all_sps) 

# head(all_sps)

#####################
# 1nd) get missing species
# Species with sdms projections for all possible shifts
my_sdms <- check_if_has_sdms_for_all_shifts(
    all_sps, 
    bioshifts, 
    sdm_dir = sdm_dir_supp)
# head(my_sdms)

# these are all possible sdms (species + ID)
nrow(my_sdms)

# these are the ones with projections for all years
I_have_sdms <- my_sdms[which(my_sdms$I_have_sdms),]
nrow(I_have_sdms)  

# N species I have sdms with projections for all years
length(unique(I_have_sdms$sps)) 

# missing species (== dont have sdms projected for all years)
all_sps <- my_sdms[which(!my_sdms$I_have_sdms),]
length(unique(all_sps$sps))
head(all_sps)

#####################
# 3rd) submit jobs
N_jobs_at_a_time = 100
N_Nodes = 1
tasks_per_core = 1

for(i in 1:nrow(all_sps)){
    
    sptogo <- all_sps$sps[i]
    
    # check if job for this species is running
    test_run <- system("squeue --format='%.50j' --me", intern = TRUE)
    test_run <- gsub(" ","",test_run)
    
    # Submit job
    if(!sptogo %in% test_run){
        
        realm <- all_sps$realm[i]
        
        args = paste(sptogo, realm)
        
        if(realm == "Mar"){ # for the Marine use this
            cores = 20
            time = "24:00:00"
            memory = "16G"
            partition = "normal"
        }
        if(realm == "Ter"){ # for the Terrestrial use this (bigger jobs) 
            # cores = 10 
            # reduce N cores due and increase memory to avoid out-of-memory issue
            cores = 5 
            memory = "100G"
            time = "1-24:00:00"
            partition = select_partition(
                request_mem = as.numeric(gsub("[^0-9.-]", "", memory)), 
                request_cpu = cores, 
                limits=limits)
        }
        
        slurm_job_singularity(jobdir = jobdir,
                              logdir = logdir, 
                              sptogo = sptogo, 
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



