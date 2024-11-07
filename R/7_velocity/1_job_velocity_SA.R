
### P.S.: Cannot run this inside singularity container
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()

# res_raster <- "1km"
res_raster <- "25km"

library(tictoc)
library(here)
library(dplyr)

########################
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
source("R/my_functions.R")

# create dir for log files
logdir <- here::here(velocity_script_dir,"slurm-log")
if(!dir.exists(logdir)){
    dir.create(logdir,recursive = TRUE)
}

# create dir for job files
jobdir <- here::here(velocity_script_dir,"job")
if(!dir.exists(jobdir)){
    dir.create(jobdir,recursive = TRUE)
}

Rscript_file = here::here(velocity_script_dir,"1_get_velocity_SA.R")

# Load study areas v3
v3_polygons <- gsub(".shp","",list.files(SA_shps_dir,pattern = ".shp"))
length(v3_polygons)

Bioshifts_DB <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v3))
Bioshifts_DB <- bioshifts_fix_columns(Bioshifts_DB)
Bioshifts_DB <- Bioshifts_DB %>%
    mutate(new_eco = ifelse(
        Eco == "Mar", "Mar", ifelse(Eco != "Mar", "Ter", NA))) %>%
    mutate(Eco = new_eco)
Bioshifts_DB <- Bioshifts_DB %>%
    dplyr::select(ID, Eco) %>%
    unique()

length(unique(Bioshifts_DB$ID))
nrow(Bioshifts_DB)

all(Bioshifts_DB$ID %in% v3_polygons)

##############################
# Rerun just for marines
# Bioshifts_DB <- Bioshifts_DB[which(Bioshifts_DB$Eco=="Mar"),]
##############################
# Rerun just for terrestrials
# Bioshifts_DB <- Bioshifts_DB[which(Bioshifts_DB$Eco=="Ter"),]
##############################

any(!Bioshifts_DB$ID %in% v3_polygons)

# Filter Polygons in Study areas v3
v3_polygons <- v3_polygons[v3_polygons %in% Bioshifts_DB$ID]

I_have <- sapply(1:nrow(Bioshifts_DB), function(x){
    ID <- Bioshifts_DB$ID[x]
    if(res_raster=="25km"){
        ID <- paste(ID,res_raster,sep="_")
    }
    test1 <- file.exists(here::here(velocity_SA_dir, paste0(ID,".csv")))
    if(test1){
        Eco <- Bioshifts_DB$Eco[x]
        if(Eco=="Ter"){
            test2 <- list.files(here::here(velocity_SA_dir),pattern = ID)
            var <- "mat"
            if(res_raster=="25km"){
                var <- paste(var,res_raster,sep="_")
            }
            test2 <- any(grepl(var,test2))
            test3 <- any(grepl("sst",list.files(here::here(velocity_SA_dir),pattern = ID)))
            if(test3){
                tmp <- list.files(here::here(velocity_SA_dir),pattern = ID, full.names = TRUE)
                unlink(tmp)
            }
        } else {
            test2 <- list.files(here::here(velocity_SA_dir),pattern = ID)
            test2 <- any(grepl("sst",test2))
            test3 <- any(grepl("mat",list.files(here::here(velocity_SA_dir),pattern = ID)))
            if(test3){
                tmp <- list.files(here::here(velocity_SA_dir),pattern = ID, full.names = TRUE)
                unlink(tmp)
            }
        }
    } else {
        test2 <- FALSE
    }
    length(which(c(test1,test2)))==2
})
length(which(I_have))

# missing 
missing <- Bioshifts_DB[which(!I_have),]
nrow(missing)
head(missing)

########################
# submit jobs

N_jobs_at_a_time = 50
N_Nodes = 1
tasks_per_core = 1

# Check if file exists
check_if_file_exists <- FALSE

for(i in 1:nrow(Bioshifts_DB)){
    
    SAtogo <- Bioshifts_DB$ID[i]
    ECO <- Bioshifts_DB$Eco[i]
    
    args = c(SAtogo, ECO, res_raster)
    
    ########################
    # Check if file exists
    
    if(res_raster=="1km"){
        file.test <- here::here(velocity_SA_dir, paste0(SAtogo,".csv"))
    } else {
        file.test <- here::here(velocity_SA_dir, paste0(paste(SAtogo,res_raster,sep = "_"),".csv"))
    }
    
    if(check_if_file_exists){
        RUN <- !file.exists(file.test)
    } else {
        RUN <- TRUE
    }
    if(RUN){
        
        # check if job for this SA is running
        test <- system("squeue --format='%.50j' --me", intern = TRUE)
        test <- gsub(" ","",test)
        
        if(!SAtogo %in% test){
            
            if(ECO == "Mar"){ # for the Marine use this
                cores = 1
                time = "24:00:00"
                memory = "8G"
                partition = "bigmem-amd"
            }
            if(ECO == "Ter"){ # for the Terrestrial use this (bigger jobs) 
                cores = 5 
                memory = "500G"
                time = "1-24:00:00"
                partition = select_partition(request_mem = as.numeric(gsub("[^0-9.-]", "", memory)), 
                                             request_cpu = cores, 
                                             limits=limits)
            }
            
            slurm_job_singularity(jobdir = jobdir,
                                  logdir = logdir, 
                                  sptogo = SAtogo, 
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
}


# check if I got velocities for all SA
# list of SA we got data
SA_got <- list.files(velocity_SA_dir,pattern = "csv")
SA_got <- gsub(".csv","",SA_got)
missing_SA <- Bioshifts_DB$Name[!Bioshifts_DB$Name %in% SA_got]
missing_SA

v3_polygons <- missing_SA
Bioshifts_DB <- Bioshifts_DB[Bioshifts_DB$Name %in% v3_polygons,]

head(Bioshifts_DB)

# Check error file
error_f <- lapply(missing_SA, function(x){
    tmp <- read.csv(here::here(velocity_script_dir,"slurm-log",paste0(x,".err")))
    tmp[nrow(tmp),]
})
error_f
error_f[1]

# A100_P1, A171_P1, A191_P1 = too small area to calculate velocities (not enough grid-cells)
# A220_P1 = period outside the range of the environmental data