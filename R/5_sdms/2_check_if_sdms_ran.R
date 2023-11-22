
rm(list=ls())
gc()


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

source(here::here("R/settings.R"))

work_dir <- getwd()
sps.dir <- here::here(work_dir,"Data/Env_data")
sdm_dir <- here::here(scratch_dir,"SDMs")
script.dir <- here::here(work_dir,"R/5_run_sdms")

Rscript_file = here::here(script.dir,"1_run_sdm.R")

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

# sp list with environmental data
all_sps <- list.files(here::here(sps.dir), pattern = '.qs', recursive = TRUE)
all_sps <- gsub("_"," ",all_sps)
all_sps <- gsub(".qs","",all_sps)
all_sps <- sapply(all_sps, function(x){
    tmp <- strsplit(x,"/")[[1]][2]
    tmp <- strsplit(tmp," ")[[1]][1:2]
    return(paste(tmp,collapse = " "))
})

N_OCC <- read.csv("Data/n_occ.csv")

# select a list of representative terrestrial species
# or run for all terrestrial species (option == "all")

option <- "all"

if(option == "all"){
    terrestrials <- N_OCC$scientificName[which(N_OCC$ECO == "T")]
} else {
    my_terrestrial_list <- readRDS("my_list_sps.RDS")
    N_OCC_T <- N_OCC[N_OCC$scientificName %in% my_terrestrial_list,]
    terrestrials <- N_OCC_T$scientificName
}

marines <- N_OCC$scientificName[which(N_OCC$ECO == "M")]

mar_sps <- all_sps[which(all_sps %in% marines)]
mar_sps <- data.frame(sps = mar_sps, realm = "Mar")

ter_sps <- all_sps[which(all_sps %in% terrestrials)]
ter_sps <- data.frame(sps = ter_sps, realm = "Ter")


## get list of marine species I have SDMs
sdms_mar <- sapply(1:nrow(mar_sps), function(i) {
    
    # check if dir exists
    tmp = dir.exists(here::here(sdm_dir,
                                "Mar",
                                gsub(" ","_",mar_sps$sps[i]),
                                gsub(" ",".",mar_sps$sps[i])))
    
    if(!tmp){
        tmp
    } else {
        # check if there a file with ensemble model outputs
        tmp = list.files(here::here(sdm_dir,
                                    "Mar",
                                    gsub(" ","_",mar_sps$sps[i]),
                                    gsub(" ",".",mar_sps$sps[i])))
        tmp = any(grepl("ensemble.models.out",tmp))
        return(tmp)
    }
    
})

if(any(sdms_mar)){ # any missing dir?
    cat("There are", length(sdms_mar), "marine species:", 
        "\nSDMs fitted for:",length(which(sdms_mar)), "species",
        "\nSDMs missing for:",length(which(!sdms_mar)), "species\n")
} else {
    cat("I have SDMs for all",length(sdms_mar), "marine species from which environmental data was extracted")
}



## get list of terrestrial species I have SDMs
sdms_ter <- sapply(1:nrow(ter_sps), function(i) {
    
    # check if dir is missing
    tmp = dir.exists(here::here(sdm_dir,
                                 "Ter",
                                 gsub(" ","_",ter_sps$sps[i]),
                                 gsub(" ",".",ter_sps$sps[i])))
    
    if(!tmp){
        tmp
    } else {
        # check if there a file with ensemble model outputs
        tmp = list.files(here::here(sdm_dir,
                                    "Ter",
                                    gsub(" ","_",ter_sps$sps[i]),
                                    gsub(" ",".",ter_sps$sps[i])))
        tmp = any(grepl("ensemble.models.out",tmp))
    }
    
})

if(any(sdms_ter)){ # any missing dir?
    cat("There are", length(sdms_ter), "terrestrial species:", 
        "\nSDMs fitted for:",length(which(sdms_ter)), "species",
        "\nSDMs missing for:",length(which(!sdms_ter)), "species\n")
} else {
    cat("I have SDMs for all",length(sdms_ter), "terrestrial species from which environmental data was extracted")
}

# missing species
missing_sps <- rbind(
    data.frame(realm = "Mar",
               sps = mar_sps$sps[which(!sdms_mar)]),
    data.frame(realm = "Ter",
               sps = ter_sps$sps[which(!sdms_ter)])
)

dim(missing_sps)

########################
# submit jobs

N_jobs_at_a_time = 100

N_Nodes = 1
tasks_per_core = 1
cores = 28
time = "40:00"
memory = "64G"

for(i in 1:nrow(missing_sps)){
    
    sptogo <- missing_sps$sps[i]
    sptogo <- gsub(" ","_",sptogo)
    
    realm <- missing_sps$realm[i]
    
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




