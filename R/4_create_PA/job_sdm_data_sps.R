
### P.S.: Cannot run this inside singularity container
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)

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

script.dir <- here::here(work_dir,"R/4_create_PA")

# source settings
source("R/settings.R")

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
all_sps <- list.files(here::here(occ_dir), pattern = '.qs')
all_sps <- gsub("_"," ",all_sps)
all_sps <- gsub(".qs","",all_sps)
length(all_sps)

# get info from v1
biov1 <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v1), header = T)
# only LAT
biov1$Type[which(biov1$Type=="HOR")] <- "LAT"
biov1 <- biov1[which(biov1$Type=="LAT"),]
nrow(biov1)

biov1 <- biov1[which(biov1$sp_name_std_v1 %in% gsub(" ","_",all_sps)),]

all_sps <- data.frame(sps = biov1$sp_name_std_v1,
                      realm = biov1$ECO)
all_sps <- all_sps[-which(duplicated(all_sps$sps)),]

all_sps$realm[which(all_sps$realm == "T")] = "Ter"
all_sps$realm[which(all_sps$realm == "M")] = "Mar"



# ECO
ter_sps <- all_sps[which(all_sps$realm == "Ter"),]
mar_sps <- all_sps[which(all_sps$realm == "Mar"),]

all_sps <- rbind(mar_sps,
                 ter_sps)

table(all_sps$realm)

# all_sps <- ter_sps # running only for terrestrials

nrow(all_sps)

########################
# submit jobs

N_jobs_at_a_time = 100
N_Nodes = 1
tasks_per_core = 1

# Check if file exists
check_if_file_exists <- TRUE

for(i in 1:nrow(all_sps)){
    # for(i in 1:nrow(missing_ter)){
    
    sptogo <- all_sps$sps[i]
    sptogo <- gsub(" ","_",sptogo)
    realm <- all_sps$realm[i]
    # sptogo <- missing_ter$sps[i]
    # sptogo <- gsub(" ","_",sptogo)
    # realm <- missing_ter$realm[i]
    
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
        # Start writing to this file
        sink(here::here(jobdir,paste0(sptogo,'.sh')))
        
        # the basic job submission script is a bash script
        cat("#!/bin/bash\n")
        
        cat("#SBATCH -N",N_Nodes,"\n")
        cat("#SBATCH -n",tasks_per_core,"\n")
        
        if(realm == "Mar"){ # for the Marine use this
            cores = 20
            time = "24:00:00"
            memory = "8G"
            partition = "normal-amd"
        }
        if(realm == "Ter"){ # for the Terrestrial use this (bigger jobs) 
            cores = 10 # reduce N cores because of out-of-memory issue
            time = "1-24:00:00"
            memory = "32G"
            partition = "bigmem"
        }
        
        cat("#SBATCH -c",cores,"\n")
        cat("#SBATCH --partition",partition,"\n")
        cat("#SBATCH --mem=",memory,"\n", sep="")
        cat("#SBATCH --time=",time,"\n", sep="")
        
        cat("#SBATCH --job-name=",sptogo,"\n", sep="")
        cat("#SBATCH --output=",here::here(logdir,paste0(sptogo,".out")),"\n", sep="")
        cat("#SBATCH --error=",here::here(logdir,paste0(sptogo,".err")),"\n", sep="")
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
