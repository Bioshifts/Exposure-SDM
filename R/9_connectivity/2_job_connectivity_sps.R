
### P.S.: Cannot run this inside singularity container
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)
library(dplyr)
library(data.table)

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
logdir <- here::here(connectivity_script_dir,"slurm-log")
if(!dir.exists(logdir)){
    dir.create(logdir,recursive = TRUE)
}

# create dir for job files
jobdir <- here::here(connectivity_script_dir,"job")
if(!dir.exists(jobdir)){
    dir.create(jobdir,recursive = TRUE)
}

Rscript_file = here::here(connectivity_script_dir,"2_get_connectivity_sps.R")

# get all shifts
shifts_ens_files <- list.files(shift_dir("Mar"), pattern = "_SA.csv")
length(shifts_ens_files)

shifts_ens <- lapply(shifts_ens_files, function(x) {
    info <- strsplit(x,"/")[[1]]
    info <- info[length(info)]
    info <- strsplit(info,"_")[[1]]
    Species <- paste(info[1],info[2],sep = "_")
    ID <- paste(info[3],info[4],sep = "_")
    Start <- as.numeric(strsplit(info[5],"-")[[1]][1])
    End <- as.numeric(strsplit(info[5],"-")[[1]][2])
    
    tmp <- data.frame(Species = Species, ID = ID, Start = Start, End = End)
    
    return(tmp)
})

shifts_ens <- do.call(rbind,shifts_ens)
shifts_ens_mar <- data.frame(shifts_ens)
shifts_ens_mar$ECO = "Mar"

length(unique(shifts_ens_mar$Species))
nrow(shifts_ens_mar)

shifts_ens_files <- list.files(shift_dir("Ter"), pattern = "_SA.csv")
length(shifts_ens_files)

shifts_ens <- lapply(shifts_ens_files, function(x){
    info <- strsplit(x,"/")[[1]]
    info <- info[length(info)]
    info <- strsplit(info,"_")[[1]]
    Species <- paste(info[1],info[2],sep = "_")
    ID <- paste(info[3],info[4],sep = "_")
    Start <- as.numeric(strsplit(info[5],"-")[[1]][1])
    End <- as.numeric(strsplit(info[5],"-")[[1]][2])
    
    tmp <- data.frame(Species = Species, ID = ID, Start = Start, End = End)
    
    return(tmp)
})

shifts_ens <- do.call(rbind,shifts_ens)
shifts_ens_ter <- data.frame(shifts_ens)
shifts_ens_ter$ECO = "Ter"

length(unique(shifts_ens_ter$Species))
nrow(shifts_ens_ter)

shifts_ens <- rbind(shifts_ens_mar, shifts_ens_ter)

length(unique(shifts_ens$Species))
length(unique(shifts_ens$Species[which(shifts_ens$ECO=="Ter")]))
length(unique(shifts_ens$Species[which(shifts_ens$ECO=="Mar")]))

# head(shifts_ens)

shifts_ens %>% filter(Species=="Neottia_nidus-avis")
which(shifts_ens$Species=="Neottia_nidus-avis")

# run only for terrestrials because we have no connectivity data for marine
shifts_ens <- shifts_ens %>% filter(ECO=="Ter")
nrow(shifts_ens)

# get all possible
all_possible <- sapply(1:nrow(shifts_ens), function(i){
    output_dir <- here(work_dir,paste0("Data/Connectivity_SA_edges/",shifts_ens$ECO[i]))
    time_period <- paste(shifts_ens$Start[i],shifts_ens$End[i],sep="-")
    shift_ID <- paste(shifts_ens$Species[i], shifts_ens$ID[i], time_period, shifts_ens$ECO[i], "Conn", sep = "_")
    filetosave <- here(output_dir, paste0(shift_ID,".csv"))
    return(filetosave)
})
I_have <- list.files(here(work_dir,"Data/Connectivity_SA_edges"), full.names = TRUE, recursive = TRUE, pattern = ".csv")

# get missing
shifts_ens <- shifts_ens[which(!all_possible %in% I_have),]
nrow(shifts_ens)

################################################TRUE################################################TRUETRUE
# submit jobs

N_jobs_at_a_time = 200
N_Nodes = 1
tasks_per_core = 1

# Check if file exists
check_if_file_exists <- TRUE

for(i in 1:nrow(shifts_ens)){
    
    spstogo <- shifts_ens$Species[i]
    SAtogo <- shifts_ens$ID[i]
    ECO <- shifts_ens$ECO[i]
    Start <- shifts_ens$Start[i]
    End <- shifts_ens$End[i]
    
    args = c(spstogo, SAtogo, ECO, Start, End)
    jobname <- paste(c(spstogo, SAtogo), collapse = "_")
    
    ########################
    # Check if file exists
    time_period <- paste(Start,End,sep="-")
    
    shift_ID <- paste(spstogo, SAtogo, time_period, ECO, "Conn", sep = "_")
    
    output_dir <- here(work_dir,paste0("Data/Connectivity_SA_edges/",ECO))
    file.test <- here::here(output_dir, paste0(shift_ID,".csv"))
    
    if(check_if_file_exists){
        RUN <- !file.exists(file.test)
    } else {
        RUN <- TRUE
    }
    if(RUN){
        
        # check if job for this SA is running
        test <- system("squeue --format='%.50j' --me", intern = TRUE)
        test <- gsub(" ","",test)
        
        if(!jobname %in% test){
            
            cores = 5 # reduce N cores because of out-of-memory issue
            time = "1-24:00:00"
            memory = "125G"
            partition = select_partition(request_mem = as.numeric(gsub("[^0-9.-]", "", memory)), 
                                         request_cpu = cores, 
                                         limits=limits)
            
            slurm_job_singularity(jobdir = jobdir,
                                  logdir = logdir, 
                                  sptogo = jobname, 
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
