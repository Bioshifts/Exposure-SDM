# --------------------------------------------------------
# title: "Extract average connectivity indices for each study area in the Bioshifts database"
# author: "Brunno F Oliveira & Bioshifts group"
# --------------------------------------------------------

# Calculates gradient-based climate velocities

# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","raster","sf","rgdal","Hmisc","dplyr", "data.table", "readxl")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
print(command_args)
polygontogo <- as.character(paste(command_args[1], collapse = " "))
Eco <- as.character(paste(command_args[2], collapse = " "))
# Eco <- "Mar"

# polygontogo <- "A112_P1" # Mar # South
# polygontogo <- "A43_P1" # Mar # North
# polygontogo <- "A1_P1" # Ter # North
# polygontogo <- "A10_P1" # Ter # Big # North
# polygontogo <- "A31_P1" # Ter # North # Elevation
# polygontogo <- "A67_P1" # Ter # North
# polygontogo <- "A25_P1" # Ter # North # Not big
# polygontogo <-"A116_P1"# Ter # very Big # North
# polygontogo <- "A79_P1"
# polygontogo <- "A102_P1"
cat("\rrunning polygon", polygontogo)

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
source("R/velocity_functions.R")

# create dir to store results
if(!dir.exists(connectivity_results_dir(Eco))){
    dir.create(connectivity_results_dir(Eco), recursive = TRUE)
}
if(!dir.exists(connectivity_data_dir(Eco))){
    dir.create(connectivity_data_dir(Eco), recursive = TRUE)
}
# create dir to store connectivities SA
conn_SA_dir <- here(scratch_dir,"Data/Connectivity/Connectivity_SA")
if(!dir.exists(conn_SA_dir)){
    dir.create(conn_SA_dir)
}

# N cores
ncores <- parallelly::availableCores()

########################
## Get avg connectivity at the SA

cat(Eco,"\n")

# load SA polygon
SA_i <- terra::vect(here::here(SA_shps_dir,paste0(polygontogo,".shp")))

# If terrestrial
if(Eco=="Ter"){
    
    # connectivity variables
    norm_curr <- terra::rast(here(connectivity_data_dir(Eco),"norm_curr.tif"))
    norm_curr_rc <- terra::rast(here(connectivity_data_dir(Eco),"norm_curr_rc.tif"))
    
    # Get connectivity at the study area
    norm_curr_file <- here(conn_SA_dir,paste0(polygontogo,"_norm_curr.tif"))
    if(!file.exists(norm_curr_file)){
        terra::window(norm_curr) <- ext(SA_i)
        norm_curr <- terra::mask(norm_curr, SA_i, filename = norm_curr_file)
    } else {
        norm_curr <- terra::rast(norm_curr_file)
    }
    
    norm_curr_rc_file <- here(conn_SA_dir,paste0(polygontogo,"_norm_curr_rc.tif"))
    if(!file.exists(norm_curr_rc_file)){
        terra::window(norm_curr_rc) <- ext(SA_i)
        norm_curr_rc <- terra::mask(norm_curr_rc, SA_i, filename = norm_curr_rc_file)
    } else {
        norm_curr_rc <- terra::rast(norm_curr_rc_file)
    }
    
    #######
    ## Project to equal area for more accurate statistics
    norm_curr <- terra::project(norm_curr, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    norm_curr_rc <- terra::project(norm_curr_rc, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    
    #######
    # summary connectivity over study area
    conn.mean <- as.numeric(terra::global(norm_curr, mean, na.rm=TRUE)[,1])
    conn.median <- as.numeric(terra::global(norm_curr, median, na.rm=TRUE)[,1])
    conn.sd <- as.numeric(terra::global(norm_curr, sd, na.rm=TRUE)[,1])
    
    conn.class <- table(terra::values(round(norm_curr_rc,0)))
    # 1) impeded--no movement, completely blocked by barriers or resistance
    # 2) diffuse--movement is unimpeded (good thing)
    # 3) intensified--movement is partly bottlenecked
    # 4) channelized--movement is completely bottlenecked.
    Ncell_classes <- data.frame(impeded = conn.class[1],
                                diffuse = conn.class[2],
                                intensified = conn.class[3],
                                channelized = conn.class[4])
    Ncell_classes[is.na(Ncell_classes)] <- 0
    Ncell <- sum(Ncell_classes)
    
    impeded_Ncell <- Ncell_classes$impeded
    diffuse_Ncell <- Ncell_classes$diffuse
    intensified_Ncell <- Ncell_classes$intensified
    channelized_Ncell <- Ncell_classes$channelized
    
    impeded_prop <- Ncell_classes$impeded/Ncell
    diffuse_prop <- Ncell_classes$diffuse/Ncell
    intensified_prop <- Ncell_classes$intensified/Ncell
    channelized_prop <- Ncell_classes$channelized/Ncell
    
    conn_j <- data.frame(SA = polygontogo,
                         conn.mean,
                         conn.median,
                         conn.sd,
                         impeded_Ncell,
                         diffuse_Ncell,
                         intensified_Ncell,
                         channelized_Ncell,
                         Ncell,
                         impeded_prop,
                         diffuse_prop,
                         intensified_prop,
                         channelized_prop)
    
    #######
    # save connectivity results
    write.csv(conn_j, here::here(connectivity_results_dir(Eco), paste0(polygontogo,".csv")), row.names = FALSE)
    
}

if(Eco=="Mar"){
    
    # connectivity variables
    # make raster
    if(!file.exists(here(connectivity_data_dir(Eco),"marine_connectivity.tif"))){
        norm_curr <- read_excel(here(work_dir,"Data","connectivity_marine.xlsx"))
        norm_curr <- rast(
            data.frame(x = marine_connectivity$x,
                       y = marine_connectivity$y,
                       z = marine_connectivity$current_flow),
            type="xyz")
        writeRaster(norm_curr, here(connectivity_data_dir(Eco),"marine_connectivity.tif"))
    } else {
        
        norm_curr <- rast(here(connectivity_data_dir(Eco),"marine_connectivity.tif"))
        
    }
    
    # Get connectivity at the study area
    norm_curr_file <- here(conn_SA_dir,paste0(polygontogo,"_norm_curr.tif"))
    if(!file.exists(norm_curr_file)){
        terra::window(norm_curr) <- ext(SA_i)
        norm_curr <- terra::mask(norm_curr, SA_i, filename = norm_curr_file)
    } else {
        norm_curr <- terra::rast(norm_curr_file)
    }
    
    #######
    ## Project to equal area for more accurate statistics
    norm_curr <- terra::project(norm_curr, Eckt)
    
    #######
    # summary connectivity over study area
    conn.mean <- as.numeric(terra::global(norm_curr, mean, na.rm=TRUE)[,1])
    conn.median <- as.numeric(terra::global(norm_curr, median, na.rm=TRUE)[,1])
    conn.sd <- as.numeric(terra::global(norm_curr, sd, na.rm=TRUE)[,1])
    
    conn_j <- data.frame(SA = polygontogo,
                         conn.mean,
                         conn.median,
                         conn.sd)
    
    #######
    # save connectivity results
    write.csv(conn_j, here::here(connectivity_results_dir(Eco), paste0(polygontogo,".csv")), row.names = FALSE)
    
}