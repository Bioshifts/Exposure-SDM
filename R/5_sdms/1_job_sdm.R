
### Cannot run this inside singularity
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)
library(dplyr)

# set computer
computer = "matrics"

if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}
work_dir <- getwd()

# source settings
source("R/settings.R")
source("R/my_functions.R")

sps.dir <- here::here(work_dir,"Data/Env_data")

# create dir for log files
logdir <- here::here(sdm_script_dir,"slurm-log")
if(!dir.exists(logdir)){
    dir.create(logdir,recursive = TRUE)
}

# create dir for job files
jobdir <- here::here(sdm_script_dir,"job")
if(!dir.exists(jobdir)){
    dir.create(jobdir,recursive = TRUE)
}

Rscript_file = here::here(sdm_script_dir,"1_run_sdm.R")

########################
# sp I have environmental data
all_sps <- list.files(here::here(sps.dir), pattern = '.qs', recursive = TRUE)
all_sps <- gsub(".qs","",all_sps)
all_sps <- sapply(all_sps, function(x){
    tmp <- strsplit(x,"/")[[1]][2]
    tmp <- strsplit(tmp,"_")[[1]][1:2]
    return(paste(tmp,collapse = " "))
})

length(all_sps) # 2656

#####################
# Separate terrestrials from marines
bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_sdms_selection(bioshifts)
bioshifts <- filter(bioshifts, gsub("_"," ",sp_name_std) %in% all_sps)
head(bioshifts)

terrestrials <- bioshifts$sp_name_std[which(bioshifts$Eco == "Ter")]
marines <- bioshifts$sp_name_std[which(bioshifts$Eco == "Mar")]

bioshifts_ter <- bioshifts[bioshifts$sp_name_std %in% terrestrials,]
bioshifts_mar <- bioshifts[bioshifts$sp_name_std %in% marines,]

ter_sps <- all_sps[which(all_sps %in% gsub("_"," ",terrestrials))]
ter_sps <- data.frame(sps = ter_sps, realm = "Ter")
nrow(ter_sps) # 2055

mar_sps <- all_sps[which(all_sps %in% gsub("_"," ",marines))]
mar_sps <- data.frame(sps = mar_sps, realm = "Mar")
nrow(mar_sps) # 602

# pipe line of species: first marines (because they run faster due to coarser resolution data) then terrestrials
all_sps <- rbind(mar_sps, ter_sps)
nrow(all_sps)
head(all_sps)

#####################
# 1nd) get missing species
# Species with sdms projections for all possible shifts
my_sdms <- data.frame()
for(i in 1:nrow(all_sps)) { cat(i, "from", nrow(all_sps),"\r")
    
    sp_i <- gsub(" ","_",all_sps$sps[i])
    sp_i_realm <- all_sps$realm[i]
    
    # test if have sdm for sp_i
    test <- list.files(
        here::here(sdm_dir(sp_i_realm),sp_i,gsub("_",".",sp_i)), 
        pattern = "ensemble.models.out")
    
    test <- length(test) > 0
    
    if(test){ # if yes, check if have projections for all possible shifts
        
        # all possible shifts for sp_i?
        shift_info <- filter(bioshifts, sp_name_std == sp_i)
        shift_info <- select(shift_info, c(ID, Start, End))
        shift_info <- unique(shift_info)
        
        ID_i <- unique(shift_info$ID)
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
        
        
    } else {
        tmp <- data.frame(sps = sp_i, realm = sp_i_realm, ID = NA, I_have_sdms = FALSE)
    }
    my_sdms <- rbind(my_sdms,tmp)
}
# these are all possible sdms (species + ID)
nrow(my_sdms) # 3930
# these are the ones with projections for all years
I_have_sdms <- my_sdms[which(my_sdms$I_have_sdms),]
nrow(I_have_sdms) # 1651
# N species I have sdms with projections for all years
length(unique(I_have_sdms$sps)) # 1056


# missing species (== dont have sdms projected for all years)
all_sps <- my_sdms[which(!my_sdms$I_have_sdms),]
length(unique(all_sps$sps)) # 1670

all_sps <- unique(all_sps[,1:2])

# put species with SDMs fitted on the top of the list
I_have_sdms <- all_sps[which(all_sps$sps %in% I_have_sdms$sps),] 
I_dont_have_sdms <- all_sps[which(!all_sps$sps %in% I_have_sdms$sps),] 
all_sps <- rbind(I_have_sdms,I_dont_have_sdms)
nrow(all_sps)

# put marine species on the top of the list
all_sps <- rbind(all_sps[which(all_sps$realm=="Mar"),],
                 all_sps[which(all_sps$realm=="Ter"),])

# select only marine or only terrestrials?
# all_sps <- all_sps[which(all_sps$realm=="Mar"),]
# all_sps <- all_sps[which(all_sps$realm=="Ter"),]
nrow(all_sps)

#####################
# 3rd) submit jobs
N_jobs_at_a_time = 30
N_Nodes = 1
tasks_per_core = 1

for(i in 1:nrow(all_sps)){
    
    sptogo <- all_sps$sps[i]
    
    # check if job for this species is running
    test <- system("squeue --format='%.50j' --me", intern = TRUE)
    test <- gsub(" ","",test)
    
    if(!sptogo %in% test){
        
        realm <- all_sps$realm[i]
        
        args = paste(sptogo, realm)
        
        if(realm == "Mar"){ # for the Marine use this
            cores = 20
            time = "24:00:00"
            memory = "16G"
            partition = "normal"
        }
        if(realm == "Ter"){ # for the Terrestrial use this (bigger jobs) 
            cores = 10 # reduce N cores because of out-of-memory issue
            time = "1-24:00:00"
            memory = "32G"
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
        tmp <- system("squeue -u $USER",intern = T)
        
        while(length(tmp)>N_jobs_at_a_time){
            
            Sys.sleep(10)
            tmp <- system("squeue -u $USER",intern = T)
            
        }
    }
}



