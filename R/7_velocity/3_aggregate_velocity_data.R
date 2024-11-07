
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

# get useful columns
Bioshifts_DB <- Bioshifts_DB %>%
    dplyr::select(ID, Eco, Areakm2, Type) %>%
    unique()

# load velocities
vel_SA <- list.files(velocity_SA_dir, pattern = ".csv", full.names = TRUE)
# separate 1km ter from 25km
pos <- grep("25km",vel_SA)
vel_SA_25 <- vel_SA[pos]
vel_SA <- vel_SA[-pos]

vel_SA <- lapply(vel_SA, read.csv)
vel_SA_25 <- lapply(vel_SA_25, read.csv)

vel_SA <- data.table::rbindlist(vel_SA, fill = TRUE)
vel_SA_25 <- data.table::rbindlist(vel_SA_25, fill = TRUE)

var_names <- c("baseline.mat","trend.mean.mat","trend.sd.mat",
               "v.median.mat","v.mean.mat","v.sd.mat",
               "v.lat.median.mat","v.lat.mean.mat","v.lat.sd.mat",
               "v.ele.median.mat","v.ele.mean.mat","v.ele.sd.mat",
               
               "baseline.map","trend.mean.map","trend.sd.map",
               "v.median.map","v.mean.map","v.sd.map",
               "v.lat.median.map","v.lat.mean.map","v.lat.sd.map",
               "v.ele.median.map","v.ele.mean.map","v.ele.sd.map",
               
               "baseline.sst","trend.mean.sst","trend.sd.sst",
               "v.mean.sst","v.median.sst","v.sd.sst",
               "v.lat.mean.sst","v.lat.median.sst","v.lat.sd.sst")
var_names_25 <- paste0(var_names,".25km")

var_names <- c("ID",var_names)
var_names_25 <- c("ID",var_names_25)

vel_SA <- vel_SA %>% select(var_names)
names(vel_SA)

pos <- grep("ele|sst",var_names)
vel_SA_25 <- vel_SA_25 %>% select(var_names[-pos])
colnames(vel_SA_25) <- var_names_25[-pos]

# get vel.mat 1km
pos <- which(is.na(vel_SA$baseline.mat))
vel_SA$baseline.mat[pos] <- vel_SA$baseline.sst[pos]

pos <- which(is.na(vel_SA$trend.mean.mat))
vel_SA$trend.mean.mat[pos] <- vel_SA$trend.mean.sst[pos]

pos <- which(is.na(vel_SA$trend.sd.mat))
vel_SA$trend.sd.mat[pos] <- vel_SA$trend.sd.sst[pos]

pos <- which(is.na(vel_SA$v.mean.mat))
vel_SA$v.mean.mat[pos] <- vel_SA$v.mean.sst[pos]

pos <- which(is.na(vel_SA$v.median.mat))
vel_SA$v.median.mat[pos] <- vel_SA$v.median.sst[pos]

pos <- which(is.na(vel_SA$v.sd.mat))
vel_SA$v.sd.mat[pos] <- vel_SA$v.sd.sst[pos]

pos <- which(is.na(vel_SA$v.lat.mean.mat))
vel_SA$v.lat.mean.mat[pos] <- vel_SA$v.lat.mean.sst[pos]

pos <- which(is.na(vel_SA$v.lat.median.mat))
vel_SA$v.lat.median.mat[pos] <- vel_SA$v.lat.median.sst[pos]

pos <- which(is.na(vel_SA$v.lat.sd.mat))
vel_SA$v.lat.sd.mat[pos] <- vel_SA$v.lat.sd.sst[pos]


# get vel.map 1km
pos <- which(is.na(vel_SA$baseline.map))
vel_SA$baseline.map[pos] <- vel_SA$baseline.sst[pos]

pos <- which(is.na(vel_SA$trend.mean.map))
vel_SA$trend.mean.map[pos] <- vel_SA$trend.mean.sst[pos]

pos <- which(is.na(vel_SA$trend.sd.map))
vel_SA$trend.sd.map[pos] <- vel_SA$trend.sd.sst[pos]

pos <- which(is.na(vel_SA$v.mean.map))
vel_SA$v.mean.map[pos] <- vel_SA$v.mean.sst[pos]

pos <- which(is.na(vel_SA$v.median.map))
vel_SA$v.median.map[pos] <- vel_SA$v.median.sst[pos]

pos <- which(is.na(vel_SA$v.sd.map))
vel_SA$v.sd.map[pos] <- vel_SA$v.sd.sst[pos]

pos <- which(is.na(vel_SA$v.lat.mean.map))
vel_SA$v.lat.mean.map[pos] <- vel_SA$v.lat.mean.sst[pos]

pos <- which(is.na(vel_SA$v.lat.median.map))
vel_SA$v.lat.median.map[pos] <- vel_SA$v.lat.median.sst[pos]

pos <- which(is.na(vel_SA$v.lat.sd.map))
vel_SA$v.lat.sd.map[pos] <- vel_SA$v.lat.sd.sst[pos]

vel_SA <- data.frame(vel_SA)[,-grep("sst",names(vel_SA))]

Bioshifts_DB <- merge(Bioshifts_DB, 
                      vel_SA, 
                      by = "ID", 
                      all.x = TRUE)

Bioshifts_DB <- merge(Bioshifts_DB, 
                      vel_SA_25, 
                      by = "ID", 
                      all.x = TRUE)

write.csv(Bioshifts_DB, 
          here(Bioshifts_dir,"vel_SA_all.csv"), 
          row.names = FALSE)


