# --------------------------------------------------------
# title: "Calculate climate change velocity for each study area in the Bioshifts database"
# author: "Brunno F Oliveira & Bioshifts group"
# --------------------------------------------------------


# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","elevatr","raster","sf","rgdal","Hmisc","dplyr", "data.table","enmSdmX")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)


# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

# source settings
source("R/settings.R")
source("R/bioshiftsFunction.R")

# detect slurm cores
N_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

# Eckert 4 equal-area projection
Eckt <- enmSdmX::getCRS('Eckert 4')

# create dir to store results
if(!dir.exists(velocity_SA_dir)){
    dir.create(velocity_SA_dir)
}

# Load study areas v3
v3_polygons <- list.files(SA_shps_dir,pattern = ".shp",full.names = TRUE)
v3_polygons_names <- gsub(".shp","",list.files(SA_shps_dir,pattern = ".shp"))

# Load Bioshifts v3
Bioshifts_DB <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB))

# Create Polygon IDs
Bioshifts_DB$Polygon <- paste0("A",Bioshifts_DB$ID_v1,"_",Bioshifts_DB$Study_ID_v1)
pos_v2 <- which(is.na(Bioshifts_DB$ID_v1))
Bioshifts_DB$Polygon[pos_v2] <- paste0("B",Bioshifts_DB$ID_v2[pos_v2],"_",Bioshifts_DB$Study_ID_v2[pos_v2])

# Filter Polygons in Study areas v3
Bioshifts_DB <- Bioshifts_DB %>%
    filter(Polygon %in% v3_polygons_names)

keep_poly <- which(v3_polygons_names %in% Bioshifts_DB$Polygon)
v3_polygons_names <- v3_polygons_names[keep_poly]
v3_polygons <- v3_polygons[keep_poly]

dim(Bioshifts_DB)
length(unique(Bioshifts_DB$Polygon))
length(unique(v3_polygons_names))

########################
# Loop through shape files of study areas

# i =1 
parallel::mclapply(1:length(v3_polygons), function(i) { 
    
    SA_i <- terra::vect(v3_polygons[i])
    
    Bioshifts_i <- Bioshifts_DB %>%
        filter(Polygon == v3_polygons_names[i])
    
    # Is it Terrestrial or Marine?
    # Check then load temperature data
    my_test <- if(any(is.na(SA_i$EleExtentm))){ # it is terrestrial if it has elevation data 
        ECO = "Mar"
    } else {
        ECO = "Ter"
    }
    
    # shifts start and end
    if(grepl("A",v3_polygons_names[i])){
        S_start <- round(unique(Bioshifts_i$START_v1),0)
        S_end <- round(unique(Bioshifts_i$END_v1),0)
        
    } else {
        S_start <- Bioshifts_i$START_v2
        S_end <- Bioshifts_i$END_v2
    }
    S_time <- S_start:S_end
    
    # get layers within time period of shift
    temperature_layers <- list.files(here::here(vars_dir(ECO),paste0("bio_proj_",my_res)))
    temperature_layers_names <- list.files(here::here(vars_dir(ECO),paste0("bio_proj_",my_res)))
    # temperature_layers <- list.files(here::here(vars_dir(ECO),"bio_proj_5km"),full.names = TRUE)
    # temperature_layers_names <- list.files(here::here(vars_dir(ECO),"bio_proj_5km"))
    temperature_layers_pos <- grepl(paste(S_time,collapse = "|"),temperature_layers_names)
    temperature_layers_names <- temperature_layers_names[temperature_layers_pos]
    temperature_layers <- temperature_layers[temperature_layers_pos]
    temperature_layers <- terra::rast(temperature_layers)
    # get MAT only
    temperature_layers <- temperature_layers[[which(names(temperature_layers)=="mat")]]
    names(temperature_layers) <- S_time
    # crop layers to the study area
    temperature_layers <- terra::mask(terra::crop(temperature_layers, SA_i), SA_i)
    
    # project to equal-area
    temperature_layers <- terra::project(temperature_layers,Eckt)
    
    # range 0-1
    temperature_layers <- rast(lapply(temperature_layers,range01raster))
    
    # calculate velocity
    sa_vel <- bioshifts(
        x = temperature_layers,
        times = S_time,
        cores = 2,
        metrics = c("centroid","nsCentroid","nsQuants")
    )
    
    write.csv(sa_vel, here::here(velocity_SA_dir, paste0(v3_polygons_names[i],".csv")), row.names = FALSE)
    
}, mc.cores = N_cores)
