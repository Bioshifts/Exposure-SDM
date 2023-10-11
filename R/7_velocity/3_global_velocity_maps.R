# --------------------------------------------------------
# title: "Create global velocity maps"
# author: "Brunno F Oliveira & Bioshifts group"
# --------------------------------------------------------

# Calculates gradient-based climate velocities

# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","raster","sf","rgdal","Hmisc","dplyr", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

# source settings
source("R/settings.R")
source("R/velocity_functions.R")

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
ECO <- as.character(paste(command_args[1], collapse = " "))
velocity_variable <- as.character(paste(command_args[2], collapse = " "))

# ECO="Ter"
# velocity_variable="mat"
# ECO="Mar"
# velocity_variable="mean"

# set time period
S_time <- 1960:2009
S_time_name <- paste(range(S_time), collapse = "-")

########################
## Calculate velocity

cat("Calculating for", ECO, velocity_variable, "\n")

# get climate layers
if(ECO=="Ter"){
    vars_dir <- here::here(vars_dir(ECO),paste0("bio_proj_",my_res))
}
if(ECO=="Mar"){
    vars_dir <- here::here(vars_dir(ECO),"bio_proj")
}
climate_layers <- list.files(vars_dir)
climate_layers_pos <- grepl(paste(S_time,collapse = "|"),climate_layers)
climate_layers <- climate_layers[climate_layers_pos]
climate_layers <- terra::rast(here::here(vars_dir,climate_layers))


# If terrestrial
# Calculate velocity for LAT and ELE and for each climatic variable
if(ECO=="Ter"){
    
    # get elevation
    elevation <- terra::rast(here::here(work_dir,"Data/elevation_1km.tif"))
    
    # select the climate variable
    climate_layers <- climate_layers[[which(names(climate_layers)==velocity_variable)]]
    
    # force rasters to pair
    elevation <- terra::project(elevation,climate_layers)
    
    #######
    # calculate the trend (C/year)
    cat("calculate the trend\n")
    
    ttrend_file <- here::here(work_dir,"Data",
                              paste(ECO,velocity_variable,"trend.tif", sep = "_"))
    ttrend <- try(terra::rast(ttrend_file),silent = TRUE)
    if(class(ttrend)=="try-error"){
        ttrend = temp_grad(
            climate_layers,
            th = 0.25*nlyr(climate_layers), ## set minimum N obs. to 1/4 time series length
            tempfile = ttrend_file,
            overwrite = TRUE)
    }
    
    
    #######
    # calculate the spatial gradient (C/km)
    cat("calculate the spatial gradient\n")
    
    avg_climate_layers_file <- here::here(work_dir,"Data",
                                          paste(ECO,velocity_variable,"avg_climate_layers.tif", sep = "_"))
    avg_climate_layers <- try(terra::rast(avg_climate_layers_file),silent = TRUE)
    if(class(avg_climate_layers)=="try-error"){
        avg_climate_layers <- terra::app(
            climate_layers, 
            fun = mean, na.rm = TRUE, 
            filename=avg_climate_layers_file,
            overwrite=TRUE)
    }
    
    spgrad = spatial_grad(avg_climate_layers)
    
    #######
    # Get the spatial gradient up slope 
    cat("Get the spatial gradient up slope \n")
    spgrad_ele = spatial_grad(elevation)
    
    # Convert angle to radians
    initial_angle_rad <- .rad(spgrad$angle) # angle of the spatial gradient 
    target_angle_rad <- .rad(spgrad_ele$angle) # angle of the elevation up slope
    
    # Apply conversion >> What is the environmental gradient up slope? 
    spgrad$angle_ele <- spgrad_ele$angle
    spgrad$Grad_ele <- spgrad$Grad * cos(initial_angle_rad - target_angle_rad) 
    
    #######
    ## calculate gradient-based climate velocity:
    ## Across latitude
    cat("calculate velocity across latitude\n")
    gVelLat <- gVelocity(grad = spgrad, slope = ttrend, 
                         grad_col = "NS", truncate = TRUE)
    
    ## Across elevation
    cat("calculate velocity across elevation\n")
    gVelEle <- gVelocity(grad = spgrad, slope = ttrend, 
                         grad_col = "Grad_ele", truncate = TRUE)
    
    ## Undirectional
    cat("calculate velocity undirectional\n")
    gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
    
    ## Angle
    gVelAngle <- gVel
    gVelAngle[spgrad$icell] <- spgrad$angle
    
    ## Angle upslope
    gVelAngleEle <- gVelEle
    gVelAngleEle[spgrad$icell] <- spgrad$angle_ele
    
    #######
    ## change sign of gVelLat if in the south hemisphere to reflect a velocity away of the tropics
    SouthCells <- terra::as.data.frame(gVelLat, xy = TRUE, cell = TRUE) 
    SouthCells <- SouthCells %>% filter(y<0)
    SouthCells <- SouthCells$cell
    if(length(SouthCells)>0){
        gVelLat[SouthCells] <- gVelLat[SouthCells] * -1
    }
    
    # in CHELSA, temperature is *10
    if(velocity_variable=="mat"){
        gVelLat <- gVelLat/10
        gVelEle <- gVelEle/10
        gVel <- gVel/10
    }
    
    #######
    # Save rasters
    cat("Save rasters\n")
    
    terra::writeRaster(
        gVelLat,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelLat",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVelEle,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelEle",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVel,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVel",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVelAngle,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelAngle",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVelAngleEle,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelAngleEle",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    
    
} else {
    
    # If marine
    # Calculate only for SST
    
    # select the climate variable
    climate_layers <- climate_layers[[which(names(climate_layers)==velocity_variable)]]
    
    #######
    # calculate the trend (C/year)
    cat("calculate the trend\n")
    
    ttrend_file <- here::here(work_dir,"Data",
                              paste(ECO,velocity_variable,"trend",paste0(S_time_name,".tif"), sep = "_"))
    ttrend <- try(terra::rast(ttrend_file),silent = TRUE)
    if(class(ttrend)=="try-error"){
        ttrend = temp_grad(
            climate_layers,
            th = 0.25*nlyr(climate_layers), ## set minimum # obs. to 1/4 time series length
            tempfile = ttrend_file,
            overwrite = TRUE)
    }
    
    
    #######
    # calculate the spatial gradient (C/km)
    cat("calculate the spatial gradient\n")
    
    avg_climate_layers_file <- here::here(work_dir,"Data",
                                          paste(ECO,velocity_variable,"avg_climate_layers",paste0(S_time_name,".tif"), sep = "_"))
    avg_climate_layers <- try(terra::rast(avg_climate_layers_file),silent = TRUE)
    if(class(avg_climate_layers)=="try-error"){
        avg_climate_layers <- terra::app(
            climate_layers, 
            fun = mean, na.rm = TRUE, 
            filename=avg_climate_layers_file,
            overwrite=TRUE)
    }
    
    spgrad = spatial_grad(avg_climate_layers)
    
    #######
    ## calculate gradient-based climate velocity:
    ## Across latitude
    cat("calculate velocity across latitude\n")
    gVelLat <- gVelocity(grad = spgrad, slope = ttrend, 
                         grad_col = "NS", truncate = TRUE)
    
    ## Undirectional
    cat("calculate velocity undirectional\n")
    gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
    
    ## Angle
    gVelAngle <- gVel
    gVelAngle[spgrad$icell] <- spgrad$angle
    
    #######
    ## change sign of gVelLat if in the south hemisphere to reflect a velocity away of the tropics
    SouthCells <- terra::as.data.frame(gVelLat, xy = TRUE, cell = TRUE) 
    SouthCells <- SouthCells %>% filter(y<0)
    SouthCells <- SouthCells$cell
    if(length(SouthCells)>0){
        gVelLat[SouthCells] <- gVelLat[SouthCells] * -1
    }
    
    #######
    # Save rasters
    cat("Save rasters\n")
    terra::writeRaster(
        gVelLat,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelLat",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVel,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVel",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVelAngle,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelAngle",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    
}

