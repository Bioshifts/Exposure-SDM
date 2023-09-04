### P.S.: Cannot run this inside singularity container
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()


library(tictoc)
library(here)

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
print(command_args)
polygontogo <- command_args
# polygontogo <- "A112_P1" # Mar # South
# polygontogo <- "A43_P1" # Mar # North

print(polygontogo)

cat("\rrunning polygon", polygontogo)

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

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

########################
# submit jobs

N_Nodes = 1
tasks_per_core = 1
cores = 5
time = "20:00:00"
memory = "64G"

for(i in 1:length(v3_polygons)){
    
    SAtogo <- v3_polygons[i]
    
    args = SAtogo
    
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
    cat("#SBATCH --mail-type=ALL\n")
    cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
    
    cat("IMG_DIR='/storage/simple/projects/t_cesab/brunno'\n")
    
    cat("module purge\n")
    cat("module load singularity/3.5\n")
    
    cat("singularity exec --disable-cache $IMG_DIR/brunnospatial.sif Rscript",Rscript_file,args,"\n", sep=" ")
    
    # Close the sink!
    sink()
    
    # Submit to run on cluster
    system(paste("sbatch", here::here(jobdir,paste0(SAtogo,'.sh'))))
    
}

# Check if everything went well

# did we get data for all polygons?
got_mar <- list.files(bios_SA_dir("Mar"))
got_ter <- list.files(bios_SA_dir("Ter"))
got <- c(got_ter,got_mar)

all(v3_polygons %in% got)
