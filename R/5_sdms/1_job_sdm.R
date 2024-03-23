
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
    }
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}
work_dir <- getwd()

# source settings
source("R/settings.R")

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

# sp list
all_sps <- list.files(here::here(sps.dir), pattern = '.qs', recursive = TRUE)
all_sps <- gsub("_"," ",all_sps)
all_sps <- gsub(".qs","",all_sps)
all_sps <- sapply(all_sps, function(x){
    tmp <- strsplit(x,"/")[[1]][2]
    tmp <- strsplit(tmp," ")[[1]][1:2]
    return(paste(tmp,collapse = " "))
})

length(all_sps)

# select only v1 species for now
biov1 <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v1), header = T)
length(unique(biov1$ssp))

# v1 from IDs we have selected for v3
v3_selected <- unique(c(
    list.files(bios_SA_dir("Ter")),
    list.files(bios_SA_dir("Mar"))))
biov1 <- biov1[biov1$ID %in% v3_selected,]

# only LAT
biov1$Type[which(biov1$Type=="HOR")] <- "LAT"
biov1 <- biov1[which(biov1$Type=="LAT"),]
all_sps <- all_sps[all_sps %in% gsub("_"," ",biov1$sp_name_std_v1)]
length(all_sps)

#####################
# Separate terrestrials from marines
terrestrials <- gsub("_"," ",biov1$sp_name_std_v1[which(biov1$ECO == "T")])
marines <- gsub("_"," ",biov1$sp_name_std_v1[which(biov1$ECO == "M")])

# Use temporal period from the environmental data
biov1_ter <- biov1 %>% filter(START >= temporal_range_env_data("Ter")[1] + n_yr_bioclimatic,
                              sp_name_std_v1 %in% terrestrials)
biov1_mar <- biov1 %>% filter(START >= temporal_range_env_data("Mar")[1] + n_yr_bioclimatic,
                              sp_name_std_v1 %in% marines)

biov1 <- rbind(biov1_ter,
               biov1_mar)

all_sps <- all_sps[all_sps %in% gsub("_"," ",biov1$sp_name_std_v1)]
length(all_sps)

ter_sps <- all_sps[which(all_sps %in% terrestrials)]
ter_sps <- data.frame(sps = ter_sps, realm = "Ter")
nrow(ter_sps)

mar_sps <- all_sps[which(all_sps %in% marines)]
mar_sps <- data.frame(sps = mar_sps, realm = "Mar")
nrow(mar_sps)

# pipe line of species: first marines (because they run faster due to coarser resolution data) then terrestrials
all_sps <- rbind(mar_sps, ter_sps)
# all_sps = ter_sps
# all_sps = mar_sps
nrow(all_sps)

#####################
# 1st) delete all duplicated sdms
for(i in 1:nrow(all_sps)){ 
    
    cat("\r",i, "from", nrow(all_sps))
    
    # models
    files_sdms <- list.files(
        here::here(sdm_dir(all_sps$realm[i]),gsub(" ","_",all_sps$sps[i]),gsub(" ",".",all_sps$sps[i]),"models"), 
        full.names = TRUE)  
    
    if(length(files_sdms) > 1){ # if there are > 1 model, keep the most recent one
        del <- order(file.info(files_sdms)$ctime,decreasing = TRUE)[-1]
        
        unlink(files_sdms[del],recursive = TRUE)
    }
    
    # link to models
    files_sdms_model <- list.files(
        here::here(sdm_dir(all_sps$realm[i]),gsub(" ","_",all_sps$sps[i]),gsub(" ",".",all_sps$sps[i])), 
        pattern = "models.out",
        full.names = TRUE)  
    
    if(length(files_sdms)==0){ # if there are mo models, delete all
        unlink(files_sdms_model)
    }
    
    if(length(files_sdms_model) > 2){ # if there are > 2 (simple + ensemble)
        
        which_sdm <- strsplit(files_sdms,split = "/")[[1]]
        which_sdm <- which_sdm[length(which_sdm)]
        
        del <- files_sdms_model[!grepl(which_sdm, files_sdms_model)]
        
        unlink(del)
        
    }
}

#####################
# 2nd) get missing species
# species we have env data but dont have sdms
I_have <- c()
for(i in 1:nrow(all_sps)) { cat(i, "from", nrow(all_sps),"\r")
    sp_i <- gsub(" ","_",all_sps$sps[i])
    sp_i_realm <- gsub(" ","_",all_sps$realm[i])
    # get studies for species i
    ID_i <- unique(biov1$ID[which(biov1$sp_name_std_v1 == sp_i)])
    # for each ID_i, look if has projections for all years
    I_have[i] <- all(sapply(ID_i, function(x){
        bio_i <- biov1[which(biov1$ID == x),]
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

missing_sps <- !I_have

# I have
length(which(I_have))
# 1079

# missing species
all_sps <- all_sps[missing_sps,]
nrow(all_sps)
# 670


#####################
# 3rd) submit jobs
N_jobs_at_a_time = 30
N_Nodes = 1
tasks_per_core = 1

for(i in 1:nrow(all_sps)){
    
    sptogo <- all_sps$sps[i]
    sptogo <- gsub(" ","_",sptogo)
    
    # check the job for this species is running
    test <- system("squeue --format='%.50j' --me", intern = TRUE)
    test <- gsub(" ","",test)
    
    if(!sptogo %in% test){
        
        realm <- all_sps$realm[i]
        
        args = paste(sptogo, realm)
        
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
            partition = "bigmem-amd"
        }
        
        # Start writing to this file
        sink(here::here(jobdir,paste0(sptogo,'.sh')))
        
        # the basic job submission script is a bash script
        cat("#!/bin/bash\n")
        
        cat("#SBATCH -N",N_Nodes,"\n")
        cat("#SBATCH -n",tasks_per_core,"\n")
        cat("#SBATCH -c",cores,"\n")
        cat("#SBATCH --mem=",memory,"\n", sep="")
        cat("#SBATCH --partition=",partition,"\n", sep="")
        
        cat("#SBATCH --job-name=",sptogo,"\n", sep="")
        cat("#SBATCH --output=",here::here(logdir,paste0(sptogo,".out")),"\n", sep="")
        cat("#SBATCH --error=",here::here(logdir,paste0(sptogo,".err")),"\n", sep="")
        cat("#SBATCH --time=",time,"\n", sep="")
        # cat("#SBATCH --mail-type=ALL\n")
        # cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
        
        cat(paste0("IMG_DIR='",singularity_image,"'\n"))
        
        cat("module purge\n")
        cat("module load singularity\n")
        
        cat("singularity exec --disable-cache $IMG_DIR Rscript",Rscript_file, args,"\n", sep=" ")
        
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



