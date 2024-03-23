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

# source settings
source("R/settings.R")

script.dir <- here::here(work_dir,"R/2_2_bios_SA")

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

Rscript_file = here::here(script.dir,"get_bios_SA.R")

# Load study areas v3
v3_polygons <- gsub(".shp","",list.files(SA_shps_dir,pattern = ".shp"))

Bioshifts_DB_v1 <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v1))
Bioshifts_DB_v2 <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v2))
Bioshifts_DB_v2$ID <- paste0("B",Bioshifts_DB_v2$Paper.ID,"_",Bioshifts_DB_v2$Study.Period)
Bioshifts_DB_v2$ECO <- ifelse(Bioshifts_DB_v2$Ecosystem.Type=="marine","M","T")
Bioshifts_DB_v2$START <- Bioshifts_DB_v2$Start.Year
Bioshifts_DB_v2$END <- Bioshifts_DB_v2$End.Year

Bioshifts_DB_v1 <- Bioshifts_DB_v1[,c("ID","ECO","START","END")]
Bioshifts_DB_v2 <- Bioshifts_DB_v2[,c("ID","ECO","START","END")]

Bioshifts_DB <- rbind(Bioshifts_DB_v1,
                      Bioshifts_DB_v2)
Bioshifts_DB$START <- round(Bioshifts_DB$START,0)
Bioshifts_DB$END <- round(Bioshifts_DB$END,0)

Bioshifts_DB <- Bioshifts_DB[-which(duplicated(Bioshifts_DB$ID)),]
# head(Bioshifts_DB)

# all(v3_polygons %in% Bioshifts_DB$ID)
# v3_polygons[!v3_polygons %in% Bioshifts_DB$ID]
# # these are new polygons added

##############################
# # Rerun just for terrestrials
# 
# Bioshifts_DB <- Bioshifts_DB[which(Bioshifts_DB$ECO=="T"),]
# 
# # Filter Polygons in Study areas v3
# v3_polygons <- v3_polygons[v3_polygons %in% Bioshifts_DB$ID]

########################
# Check what is missing

# did we get data for all polygons?
got_mar <- list.files(bios_SA_dir("Mar"))
got_ter <- list.files(bios_SA_dir("Ter"))
got <- rbind(data.frame(realm="T",poly=got_ter),
             data.frame(realm="M",poly=got_mar))

I_have <- c()
for(i in 1:nrow(got)) { cat(i, "from", nrow(got),"\r")
    poly_i <- got$poly[i]
    ECO_i <-got$realm[i]
    
    if(ECO_i=="T"){ 
        got_from_poly_i <- list.files(here::here(bios_SA_dir("Ter"),poly_i))
    } else {
        got_from_poly_i <- list.files(here::here(bios_SA_dir("Mar"),poly_i))
    }
    
    years_from_poly_i <- Bioshifts_DB$START[Bioshifts_DB$ID==poly_i]:Bioshifts_DB$END[Bioshifts_DB$ID==poly_i]
    years_I_have <- sapply(years_from_poly_i, function(x){
        any(grepl(x,got_from_poly_i))
    })
    I_have[i] <- all(years_I_have)
}
# I_have

I_have <- got$poly[I_have]

# to go 
v3_polygons <- v3_polygons[!v3_polygons %in% I_have]
length(v3_polygons)

########################
# submit jobs

N_jobs_at_a_time = 100
N_Nodes = 1
tasks_per_core = 1
cores = 1
time = "20:00:00"
memory = "32G"

for(i in 1:length(v3_polygons)){
    
    SAtogo <- v3_polygons[i]
    args = SAtogo
    
    # check the job is running
    test <- system("squeue --format='%.50j' --me", intern = TRUE)
    test <- gsub(" ","",test)
    
    if(!SAtogo %in% test){
        
        # Start writing to this file
        sink(here::here(jobdir,paste0(SAtogo,'.sh')))
        
        # the basic job submission script is a bash script
        cat("#!/bin/bash\n")
        
        cat("#SBATCH -N",N_Nodes,"\n")
        cat("#SBATCH -n",tasks_per_core,"\n")
        cat("#SBATCH -c",cores,"\n")
        
        cat("#SBATCH --job-name=",SAtogo,"\n", sep="")
        cat("#SBATCH --output=",here::here(logdir,paste0(SAtogo,".out")),"\n", sep="")
        cat("#SBATCH --error=",here::here(logdir,paste0(SAtogo,".err")),"\n", sep="")
        cat("#SBATCH --time=",time,"\n", sep="")
        cat("#SBATCH --mem=",memory,"\n", sep="")
        # cat("#SBATCH --mem=","100G","\n", sep="") # for large jobs that fail
        # cat("#SBATCH --mail-type=ALL\n")
        # cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
        if(computer == "matrics"){
            # cat("#SBATCH --partition=normal\n")
            cat("#SBATCH --partition=bigmem\n")# for large jobs that fail
        }
        
        cat("IMG_DIR='/storage/simple/projects/t_cesab/brunno'\n")
        
        if(computer == "muse"){
            cat("IMG_DIR='/storage/simple/projects/t_cesab/brunno'\n")
            
            cat("module purge\n")
            cat("module load singularity/3.5\n")
            
            cat("singularity exec --disable-cache $IMG_DIR/brunnospatial.sif Rscript",Rscript_file, args,"\n", sep=" ")
        }
        if(computer == "matrics"){
            cat(paste0("IMG_DIR='",singularity_image,"'\n"))
            
            cat("module purge\n")
            cat("module load singularity\n")
            
            cat("singularity exec --disable-cache $IMG_DIR Rscript",Rscript_file, args,"\n", sep=" ")
        }
        
        # Close the sink!
        sink()
        
        # Submit to run on cluster
        system(paste("sbatch", here::here(jobdir,paste0(SAtogo,'.sh'))))
        
        # check how many jobs in progress
        tmp <- system("squeue -u $USER",intern = T)
        
        while(length(tmp)>N_jobs_at_a_time){
            
            Sys.sleep(10)
            tmp <- system("squeue -u $USER",intern = T)
            
        }
    }
}


