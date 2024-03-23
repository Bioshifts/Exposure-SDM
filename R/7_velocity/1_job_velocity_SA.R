
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
    work_dir <- getwd()
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
    work_dir <- getwd()
}

# source settings
source("R/settings.R")

script.dir <- velocity_SA_scrit_dir

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

Rscript_file = here::here(script.dir,"1_get_velocity_SA.R")

# Load study areas v3
v3_polygons <- gsub(".shp","",list.files(SA_shps_dir,pattern = ".shp"))


Bioshifts_DB_v1 <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v1))
Bioshifts_DB_v2 <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v2))
Bioshifts_DB_v2$ID <- paste0("B",Bioshifts_DB_v2$Paper.ID,"_",Bioshifts_DB_v2$Study.Period)
Bioshifts_DB_v2$ECO <- ifelse(Bioshifts_DB_v2$Ecosystem.Type=="marine","M","T")

Bioshifts_DB_v1 <- Bioshifts_DB_v1[,c("ID","ECO")]
Bioshifts_DB_v2 <- Bioshifts_DB_v2[,c("ID","ECO")]

Bioshifts_DB <- rbind(Bioshifts_DB_v1,
                      Bioshifts_DB_v2)
Bioshifts_DB <- Bioshifts_DB[-which(duplicated(Bioshifts_DB$ID)),]

##############################
# Rerun just for marines
# Bioshifts_DB <- Bioshifts_DB[which(Bioshifts_DB$ECO=="M"),]
##############################
# Rerun just for terrestrials
# Bioshifts_DB <- Bioshifts_DB[which(Bioshifts_DB$ECO=="T"),]
##############################

# Filter Polygons in Study areas v3
v3_polygons <- v3_polygons[v3_polygons %in% Bioshifts_DB$ID]

# Get polygons metadata
v3_polygons_metadata <- read.csv(here::here(Bioshifts_dir,"Geodatabase_Bioshiftsv3_metadata.csv"))
v3_polygons_metadata <- v3_polygons_metadata[v3_polygons_metadata$Name %in% v3_polygons,]
v3_polygons_metadata$ECO <- ifelse(is.na(v3_polygons_metadata$EleExtentm),"Mar","Ter")

########################
# submit jobs

N_jobs_at_a_time = 50

# Check if file exists
check_if_file_exists <- FALSE

for(i in 1:nrow(v3_polygons_metadata)){
    
    SAtogo <- v3_polygons_metadata$Name[i]
    ECO <- v3_polygons_metadata$ECO[i]
    
    args = SAtogo
    
    ########################
    # Check if file exists
    
    file.test <- here::here(velocity_SA_dir, paste0(SAtogo,".csv"))
    
    if(check_if_file_exists){
        RUN <- !file.exists(file.test)
    } else {
        RUN <- TRUE
    }
    if(RUN){
        # Start writing to this file
        sink(here::here(jobdir,paste0(SAtogo,'.sh')))
        
        # the basic job submission script is a bash script
        cat("#!/bin/bash\n")
        
        cat("#SBATCH -N 1\n")
        cat("#SBATCH -n 1\n")
        
        if(ECO=="Ter"){
            cat("#SBATCH --time=2-20:00:00\n")
            cat("#SBATCH --partition=bigmem-amd\n")
            cat("#SBATCH -c 5\n")
            # cat("#SBATCH --mem=100G\n")
            cat("#SBATCH --mem=500G\n") # for big jobs that crash
        } else {
            cat("#SBATCH --partition=normal\n")
            cat("#SBATCH --mem=50G\n")
            cat("#SBATCH --time=1-24:00:00\n")
            cat("#SBATCH -c 1\n")
        }
        
        cat("#SBATCH --job-name=",SAtogo,"\n", sep="")
        cat("#SBATCH --output=",here::here(logdir,paste0(SAtogo,".out")),"\n", sep="")
        cat("#SBATCH --error=",here::here(logdir,paste0(SAtogo,".err")),"\n", sep="")
        # cat("#SBATCH --mail-type=ALL\n")
        # cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
        
        cat(paste0("IMG_DIR='",singularity_image,"'\n"))
        
        cat("module purge\n")
        cat("module load singularity\n")
        
        cat("singularity exec --disable-cache $IMG_DIR Rscript",Rscript_file, args,"\n", sep=" ")
        
        # Close the sink!
        sink()
        
        # Submit to run on cluster
        system(paste("sbatch", here::here(jobdir, paste0(SAtogo,'.sh'))))
        
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


# check if I got velocities for all SA
# list of SA we got data
SA_got <- list.files(velocity_SA_dir,pattern = "csv")
SA_got <- gsub(".csv","",SA_got)
missing_SA <- v3_polygons_metadata$Name[!v3_polygons_metadata$Name %in% SA_got]
missing_SA

v3_polygons <- missing_SA
v3_polygons_metadata <- v3_polygons_metadata[v3_polygons_metadata$Name %in% v3_polygons,]

head(v3_polygons_metadata)

# Check error file
error_f <- lapply(missing_SA, function(x){
    tmp <- read.csv(here::here(script.dir,"slurm-log",paste0(x,".err")))
    tmp[nrow(tmp),]
})
error_f


