
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

Rscript_file = here::here(script.dir,"1_get_shifts_edge_method.R")

########################
# load functions
source("R/my_functions.R")
source("R/range_shift_functions.R")
# source settings
source("R/settings.R")


########################
# job tunning
N_jobs_at_a_time = 50

cores = 1
partition = "normal"
time = "20:00:00"
N_Nodes = 1
tasks_per_core = 1


########################
# sp list I have sdms
all_sps_ter <- list.dirs(here::here(sdm_dir("Ter")), recursive = FALSE, full.names = FALSE)
length(all_sps_ter)

test <- sapply(all_sps_ter, function(x){
    file.exists(here::here(sdm_dir("Ter"),x,paste(x,"shift_info.csv")))
})
all_sps_ter <- all_sps_ter[test]

all_sps_mar <- list.dirs(here::here(sdm_dir("Mar")), recursive = FALSE, full.names = FALSE)
length(all_sps_mar)

test <- sapply(all_sps_mar, function(x){
    file.exists(here::here(sdm_dir("Mar"),x,paste(x,"shift_info.csv")))
})
all_sps_mar <- all_sps_mar[test]

all_sps <- c(all_sps_mar,all_sps_ter)

#####################
# Separate terrestrials from marines
bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_sdms_selection(bioshifts)
bioshifts <- filter(bioshifts, sp_name_std %in% all_sps)
# head(bioshifts)

terrestrials <- bioshifts$sp_name_std[which(bioshifts$Eco == "Ter")]
marines <- bioshifts$sp_name_std[which(bioshifts$Eco == "Mar")]

ter_sps <- all_sps[which(all_sps %in% terrestrials)]
ter_sps <- data.frame(sps = ter_sps, realm = "Ter")
nrow(ter_sps) 

mar_sps <- all_sps[which(all_sps %in% marines)]
mar_sps <- data.frame(sps = mar_sps, realm = "Mar")
nrow(mar_sps) 

# pipe line of species: first marines (because they run faster due to coarser resolution data)
all_sps <- rbind(mar_sps, ter_sps)
nrow(all_sps) 

# 1nd) get missing species
# Species with sdms projections for all possible shifts
my_sdms <- check_if_has_sdms_for_all_shifts(all_sps,bioshifts)
# head(my_sdms)

# these are all possible sdms (species + ID)
nrow(my_sdms)

# these are the ones with projections for all years
I_have_sdms <- my_sdms[which(my_sdms$I_have_sdms),]
nrow(I_have_sdms)  

head(I_have_sdms)  

########################
# sp list I have shift estimated
I_have_shift_ter <- list.files(shift_dir("Ter"), pattern = "SA_edges.csv")
I_have_shift_ter <- lapply(I_have_shift_ter, function(x){
    tmp <- strsplit(x,"_")[[1]]
    data.frame(sps=paste(tmp[1],tmp[2],sep="_"),realm=tmp[6],ID=paste(tmp[3],tmp[4],sep="_"))
})
I_have_shift_ter <- do.call(rbind,I_have_shift_ter)
I_have_shift_ter <- unique(I_have_shift_ter)

I_have_shift_mar <- list.files(shift_dir("Mar"), pattern = "SA_edges.csv")
I_have_shift_mar <- lapply(I_have_shift_mar, function(x){
    tmp <- strsplit(x,"_")[[1]]
    data.frame(sps=paste(tmp[1],tmp[2],sep="_"),realm=tmp[6],ID=paste(tmp[3],tmp[4],sep="_"))
})
I_have_shift_mar <- do.call(rbind,I_have_shift_mar)
I_have_shift_mar <- unique(I_have_shift_mar)

I_have_shift <- rbind(I_have_shift_mar,I_have_shift_ter)

########################
# missing species
# select from the list of species I have sdms the ones which we still did not have shift calculated
missing <- I_have_sdms[which(!paste(I_have_sdms$sps,I_have_sdms$ID) %in% paste(I_have_shift$sps,I_have_shift$ID)),]
nrow(missing)
# 11/21/24 - 1413

missing <- missing %>% distinct(sps, ID, .keep_all = TRUE)

head(missing)

########################
# submit jobs

if(nrow(missing)>0){
    
    for(i in 1:nrow(missing)){ 
        
        sptogo <- missing$sps[i]
        realmtogo <- missing$realm[i]
        IDtogo <- missing$ID[i]
        
        cat("Running", i, "from", nrow(missing), "species\n")
        
        args <- paste(sptogo,realmtogo,IDtogo)
        
        shifttogo = gsub(" ","_",args)
        
        # check if job for this species is running
        test <- system("squeue --format='%.50j' --me", intern = TRUE)
        test <- gsub(" ","",test)
        
        if(!shifttogo %in% test){
            
            if(realmtogo == "Mar"){ # for the Marine use this
                memory = "16G"
                partition = "normal"
            }
            if(realmtogo == "Ter"){ # for the Terrestrial use this (bigger jobs) 
                memory = "500G"
                partition = select_partition(request_mem = as.numeric(gsub("[^0-9.-]", "", memory)), 
                                             request_cpu = cores, 
                                             limits=limits)
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
        }
        
        # check how many jobs in progress
        tmp <- system("squeue -u $USER",intern = T)
        
        while(length(tmp)>N_jobs_at_a_time){
            
            Sys.sleep(10)
            tmp <- system("squeue -u $USER",intern = T)
            
        }
    }
}


