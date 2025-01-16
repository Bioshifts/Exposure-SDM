
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
    work_dir <- getwd()
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
    work_dir <- getwd()
}

# source settings
source("R/settings.R")
source("R/my_functions.R")

res_raster <- "1km"
# res_raster <- "25km"

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


any(!Bioshifts_DB$ID %in% v3_polygons)

# Filter Polygons in Study areas v3
v3_polygons <- v3_polygons[v3_polygons %in% Bioshifts_DB$ID]


Bioshifts_DB_mar <- Bioshifts_DB[which(Bioshifts_DB$Eco=="Mar"),]
Bioshifts_DB_ter <- Bioshifts_DB[which(Bioshifts_DB$Eco=="Ter"),]

##############################
# make a SA x Eco x resolution files
Bioshifts_DB_mar <- expand.grid(ID = Bioshifts_DB_mar$ID,
                                resolution = c("25km","50km","110km"))
Bioshifts_DB_mar$Eco <- "Mar"

Bioshifts_DB_ter <- expand.grid(ID = Bioshifts_DB_ter$ID,
                                resolution = c("1km","25km","50km","110km"))
Bioshifts_DB_ter$Eco <- "Ter"

# combine marine and terrestrial
Bioshifts_DB <- rbind(Bioshifts_DB_mar,
                      Bioshifts_DB_ter)

##############################
# Run just for marines
# Bioshifts_DB <- Bioshifts_DB_mar
##############################
# Rerun just for terrestrials
# Bioshifts_DB <- Bioshifts_DB_ter
##############################

########################
# submit jobs

N_jobs_at_a_time = 50
N_Nodes = 1
tasks_per_core = 1

# Check if file exists
check_if_file_exists <- TRUE

for(i in 1:nrow(Bioshifts_DB)){
    
    SAtogo <- as.character(Bioshifts_DB$ID[i])
    ECO <- as.character(Bioshifts_DB$Eco[i])
    resolution <- as.character(Bioshifts_DB$resolution[i])
    
    # SAtogo <- "A1_P1"
    # ECO <- "Ter"
    # resolution <- "1km"
    
    args = c(SAtogo, ECO, resolution)
    job_name <- paste(args,collapse = "_")
    
    ########################
    # Check if file exists
    
    if(check_if_file_exists){
        file.test <- here::here(velocity_SA_dir, paste0(SAtogo,"_",resolution,".csv"))
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
                memory = "200G"
                time = "1-24:00:00"
                partition = select_partition(request_mem = as.numeric(gsub("[^0-9.-]", "", memory)), 
                                             request_cpu = cores, 
                                             limits=limits)
            }
            
            slurm_job_singularity(jobdir = jobdir,
                                  logdir = logdir, 
                                  sptogo = job_name, 
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


#################################
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

#################################
# check of velocity lat is greater than undirectional velocity

# run file 1_1_check_for_errors_velocity_SA.R
missing_SA <- read.csv("errors_SAs.csv")

v3_polygons <- missing_SA$SA
Bioshifts_DB <- Bioshifts_DB[Bioshifts_DB$ID %in% v3_polygons,]

head(Bioshifts_DB)

#################################
# check N cols
all_SA <- list.files(velocity_SA_dir, pattern = '.csv')

# Ncols
# x = all_SA[459]

ncols_files <- lapply(all_SA, function(x){
    SA_data <- read.csv(here(velocity_SA_dir,x))
    
    SA_i_res <- strsplit(x,"_")[[1]][3]
    SA_i_res <- strsplit(SA_i_res,"[.]")[[1]][1]
    
    # ter or mar?
    if(any(grepl("mat",names(SA_data)))){
        Eco <- "Ter"
    } else {
        Eco <- "Mar"
    }
    data.frame(Eco = Eco, res = SA_i_res, ncol = ncol(SA_data))
})
ncols_files <- do.call(rbind,ncols_files)

ncols_files %>% group_by(Eco,res,ncol) %>% tally
