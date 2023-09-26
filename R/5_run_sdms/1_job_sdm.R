
### Cannot run this inside singularity
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    
    work_dir <- getwd()
    sps.dir <- here::here(work_dir,"Data/Env_data")
    script.dir <- here::here(work_dir,"R/5_run_sdms")
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

Rscript_file = here::here(script.dir,"1_run_sdm.R")

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
biov1 <- read.csv("Data/Bioshifts/biov1_fixednames.csv", header = T)

all_sps <- all_sps[all_sps %in% gsub("_"," ",biov1$sp_name_std_v1)]
length(all_sps)

#####################
# N occurrences
N_OCC <- read.csv("Data/n_occ.csv")

terrestrials <- N_OCC$scientificName[which(N_OCC$ECO == "T")]
marines <- N_OCC$scientificName[which(N_OCC$ECO == "M")]

mar_sps <- all_sps[which(all_sps %in% marines)]
mar_sps <- data.frame(sps = mar_sps, realm = "Mar")

ter_sps <- all_sps[which(all_sps %in% terrestrials)]
ter_sps <- data.frame(sps = ter_sps, realm = "Ter")


# pipe line of species: first marines (because they run faster due to coarser resolution data) then terrestrials
# all_sps <- rbind(mar_sps, ter_sps)
all_sps = ter_sps
# all_sps = mar_sps

########################
# 1st) delete all duplicated sdms
for(i in 1:nrow(all_sps)){ cat("\r",i, "from", nrow(all_sps))
    
    # models
    files_sdms <- list.files(
        here::here("/lustre/oliveirab/SDMs",all_sps$realm[i],gsub(" ","_",all_sps$sps[i]),gsub(" ",".",all_sps$sps[i]),"models"), 
        full.names = TRUE)  
    
    if(length(files_sdms) > 1){ # if there are > 1 model, keep the most recent one
        del <- order(file.info(files_sdms)$ctime,decreasing = TRUE)[-1]
        
        unlink(files_sdms[del],recursive = TRUE)
    }
    
    # link to models
    files_sdms_model <- list.files(
        here::here("/lustre/oliveirab/SDMs",all_sps$realm[i],gsub(" ","_",all_sps$sps[i]),gsub(" ",".",all_sps$sps[i])), 
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

########################
# 2nd) submit jobs


N_jobs_at_a_time = 500
N_Nodes = 1
tasks_per_core = 1


for(i in 1:nrow(all_sps)){
    
    sptogo <- all_sps$sps[i]
    sptogo <- gsub(" ","_",sptogo)
    
    realm <- all_sps$realm[i]
    
    args = paste(sptogo, realm)
    
    if(realm == "Mar"){ # for the Marine use this
        cores = 28
        time = "24:00:00"
        memory = "80G"
    }
    if(realm == "Ter"){ # for the Terrestrial use this (bigger jobs) 
        cores = 10 # reduce N cores because of out-of-memory issue
        time = "1-24:00:00"
        memory = "100G"
    }
    
    # Start writing to this file
    sink(here::here(jobdir,paste0(sptogo,'.sh')))
    
    # the basic job submission script is a bash script
    cat("#!/bin/bash\n")
    
    cat("#SBATCH -N",N_Nodes,"\n")
    cat("#SBATCH -n",tasks_per_core,"\n")
    cat("#SBATCH -c",cores,"\n")
    cat("#SBATCH --mem=",memory,"\n", sep="")
    
    cat("#SBATCH --job-name=",sptogo,"\n", sep="")
    cat("#SBATCH --output=",here::here(logdir,paste0(sptogo,".out")),"\n", sep="")
    cat("#SBATCH --error=",here::here(logdir,paste0(sptogo,".err")),"\n", sep="")
    cat("#SBATCH --time=",time,"\n", sep="")
    cat("#SBATCH --mail-type=ALL\n")
    cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
    
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



