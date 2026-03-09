
### P.S.: Cannot run this inside singularity container
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)
library(dplyr)

########################
# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}
work_dir <- getwd()

script.dir <- here::here(work_dir,"R/4_bioclimatics_sps")

# source settings
source("R/settings.R")
source("R/my_functions.R")

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

Rscript_file = here::here(script.dir,"get_bioclimatics_sps.R")

########################

# sp with GBIF occurrences
all_sps <- list.files(here::here(occ_dir), pattern = '.qs')
all_sps <- gsub(".qs","",all_sps)
all_sps <- gsub("_"," ",all_sps) 
length(all_sps) # 5154

# Load raw bioshifts
# Load v3
bioshifts <- read.csv(here::here(Bioshifts_dir, Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_sdms_selection(bioshifts)
length(unique(bioshifts$sp_name_std)) # 5495

# select species with GBIF occurrences
bioshifts <- bioshifts %>% filter(sp_name_std %in% all_sps)
length(unique(bioshifts$sp_name_std)) # 5141

# Separate terrestrials from marines
terrestrials <- unique(bioshifts$sp_name_std[which(bioshifts$Eco == "Ter")])
length(terrestrials) # 4538
marines <- unique(bioshifts$sp_name_std[which(bioshifts$Eco == "Mar")])
length(marines) # 603

# Get missing species
I_have_ter <- list.files(env_data_dir("Ter"))
I_have_ter <- gsub("_Ter.qs","",I_have_ter) 
I_have_ter <- gsub("_"," ",I_have_ter) 
I_have_ter <- I_have_ter[I_have_ter %in% terrestrials]
length(I_have_ter) 

I_have_mar <- list.files(env_data_dir("Mar"))
I_have_mar <- gsub("_Mar.qs","",I_have_mar) 
I_have_mar <- gsub("_"," ",I_have_mar) 
I_have_mar <- I_have_mar[I_have_mar %in% marines]
length(I_have_mar)

missing_ter <- terrestrials[!terrestrials %in% I_have_ter]
length(missing_ter)

missing_mar <- marines[!marines %in% I_have_mar]
length(missing_mar)

# get missing species
all_sps <- rbind(
    data.frame(realm = "Mar", sps = missing_mar),
    data.frame(realm = "Ter", sps = missing_ter)
)
# all_sps <- data.frame(realm = "Mar", sps = missing_mar)
# all_sps <- data.frame(realm = "Ter", sps = missing_ter)
nrow(all_sps)

########################
# submit jobs

N_jobs_at_a_time = 50
N_Nodes = 1
tasks_per_core = 1

# Check if file exists
check_if_file_exists <- TRUE

for(i in 1:nrow(all_sps)){
    # for(i in 1:nrow(missing_ter)){
    
    sptogo <- all_sps$sps[i]
    sptogo <- gsub(" ","_",sptogo)
    realm <- all_sps$realm[i]
    
    args = paste(sptogo, realm)
    
    ########################
    # Check if file exists
    
    file.test <- here::here(env_data_dir(realm),paste0(gsub(" ","_",sptogo),"_",realm,".qs"))
    
    if(check_if_file_exists){
        RUN <- !file.exists(file.test)
    } else {
        RUN <- TRUE
    }
    if(RUN){
        
        if(realm == "Mar"){ # for the Marine use this
            cores = 10
            time = "24:00:00"
            memory = "20G"
            partition= "normal"
        }
        if(realm == "Ter"){ # for the Terrestrial use this (bigger jobs) 
            cores = 5 # reduce N cores because of out-of-memory issue
            time = "1-24:00:00"
            memory = "100G"
            partition = "bigmem"
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


# # check if I got env data for all species
# # list of species we got data
# sps_got_ter <- list.files(env_data_dir("Mar"))
# sps_got_ter <- gsub("_Mar.qs","",sps_got_ter)
# sps_got_ter <- gsub("_"," ",sps_got_ter)
# missing_ter <- all_sps[!all_sps$sps %in% sps_got_ter,]
# nrow(missing_ter)
# head(missing_ter)
# 
# 
# "Centrostephanus rodgersii" %in% N_OCC$scientificName
# N_OCC[N_OCC$scientificName=="Centrostephanus rodgersii",]
