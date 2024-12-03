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
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}
work_dir <- getwd()

# source settings
source("R/settings.R")
source("R/my_functions.R")

script.dir <- here::here(work_dir,"R/2_2_bioclimatic_SA")

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

# Load bioshifts
bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_sdms_selection(bioshifts)

##############################
# # Rerun just for terrestrials
# 
# bioshifts <- bioshifts[which(bioshifts$Eco=="Ter"),]
# 
# # Filter Polygons in Study areas v3
bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_fix_columns(bioshifts)
bioshifts_IDs <- bioshifts %>%
    select(ID, Eco, Start, End) %>%
    unique %>%
    mutate(Start = round(Start,0),
           End = round(End,0))

v3_polygons <- unique(bioshifts_IDs$ID)

########################
# Check what is missing

# did we get data for all polygons?
got_mar <- list.files(bios_SA_dir("Mar"))
got_ter <- list.files(bios_SA_dir("Ter"))

got <- rbind(data.frame(realm="Ter",poly=got_ter),
             data.frame(realm="Mar",poly=got_mar))

I_have <- c()
for(i in 1:nrow(bioshifts_IDs)) { cat(i, "from", nrow(got),"\r")
    
    poly_i <- bioshifts_IDs$ID[i]
    Eco_i <-bioshifts_IDs$Eco[i]
    
    if(Eco_i=="Ter"){ 
        got_from_poly_i <- list.files(here::here(bios_SA_dir("Ter"),poly_i))
    } else {
        got_from_poly_i <- list.files(here::here(bios_SA_dir("Mar"),poly_i))
    }
    
    years_from_poly_i <- bioshifts_IDs$Start[which(bioshifts_IDs$ID==poly_i)]:bioshifts_IDs$End[which(bioshifts_IDs$ID==poly_i)]
    years_I_have <- sapply(years_from_poly_i, function(x){
        any(grepl(x,got_from_poly_i))
    })
    I_have[i] <- length(which(years_I_have)) == length(years_from_poly_i)
}
# I_have

bioshifts_IDs_I_have <- bioshifts_IDs$ID[which(I_have)]

# to go 
v3_polygons <- bioshifts_IDs %>% filter(!ID %in% bioshifts_IDs_I_have)

########################
# submit jobs

N_jobs_at_a_time = 100
N_Nodes = 1
tasks_per_core = 1
cores = 1
time = "20:00:00"
memory = "32G"
partition = "bigmem"

my_res = "1km"

for(i in 1:nrow(v3_polygons)){
    
    SAtogo <- v3_polygons$ID[i]
    Eco <- v3_polygons$Eco[i]
    
    # SAtogo <- "A108_P1"
    # Eco <- "Ter"
    
    args = c(SAtogo, Eco, my_res)
    
    # check the job is running
    test <- system("squeue --format='%.50j' --me", intern = TRUE)
    test <- gsub(" ","",test)
    
    if(!SAtogo %in% test){
        
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
        tmp <- system("squeue -u $USER",intern = T)
        
        while(length(tmp)>N_jobs_at_a_time){
            
            Sys.sleep(10)
            tmp <- system("squeue -u $USER",intern = T)
            
        }
    }
}


