# --------------------------------------------------------
# title: "Create global velocity maps"
# author: "Brunno F Oliveira & Bioshifts group"
# --------------------------------------------------------

# Calculates gradient-based climate velocities

# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","raster","sf","rgdal","Hmisc","dplyr", "data.table","parallelly","qs","rnaturalearth","parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
Eco <- as.character(paste(command_args[1], collapse = " "))
velocity_variable <- as.character(paste(command_args[2], collapse = " "))
res_raster <- as.character(paste(command_args[3], collapse = " "))

# Eco="Ter"
# velocity_variable="mat"
# velocity_variable="map"
# res_raster <- "1km"

# Eco="Mar"
# velocity_variable="sst"

# set time period
# period <- 2000:2015
# period <- 1979:2018
period <- 1960:2009
period_name <- paste(range(period), collapse = "-")

# N cores
ncores <- parallelly::availableCores()

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

if(Eco == "Ter"){
    my_res = res_raster
} else {
    my_res = "25km"
}


# create dir to store temporary files
tmp_dir <- here(scratch_dir,"tmp/vel_global")
if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir, recursive = TRUE)
}
# create dir to store results
vel_dir <- here(work_dir,"Data/Velocity_global")
if(!dir.exists(vel_dir)){
    dir.create(vel_dir, recursive = TRUE)
}


########################
# select environmental variable for period
climate_layers <- list.files(bios_dir(Eco), full.names = TRUE)
climate_layers <- climate_layers[grep(paste0(period,collapse = "|"),climate_layers)]
climate_layers <- rast(climate_layers)

# select the environmental variable of interest
if(Eco == "Ter"){
    climate_layers <- climate_layers[[which(names(climate_layers) %in% velocity_variable)]]
} else {
    climate_layers <- climate_layers["mean"]
    names(climate_layers) <- rep(velocity_variable, nlyr(climate_layers))
}

########################
# Get elevation data if running for terrestrial environment
if(Eco == "Ter"){
    # get elevation
    elevation <- terra::rast(here::here(work_dir,"Data/elevation_1km.tif"))
    # force raster to pair
    elevation <- terra::project(elevation, climate_layers, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
}


#######
# Calculate velocity

cat("Calculating for", Eco, velocity_variable, my_res, "\n")

variable_name <- paste(Eco,velocity_variable,my_res,period_name,"global",sep = "_")

# If terrestrial
# Calculate velocity for LAT and ELE and for each climatic variable

# select the climate variable
climate_layers <- climate_layers[velocity_variable]

#######
# calculate the trend (C/year)
cat("Trend\n")

ttrend_file_name <- paste(variable_name, "trend.tif",sep = "_")

ttrend_file <- here::here(tmp_dir,ttrend_file_name)

ttrend <- try(terra::rast(ttrend_file),silent = TRUE)

if(class(ttrend)=="try-error"){
    ttrend = temp_grad(
        climate_layers,
        th = 0.25*nlyr(climate_layers), ## set minimum N obs. to 1/4 time series length
        file_name = ttrend_file,
        overwrite = TRUE,
        ncores = ncores)
}

#######
# Get averaged climate layers
avg_climate_layers_file_name <- paste(variable_name, "avg_climate_layers.tif",sep = "_")

avg_climate_layers_file <- here::here(tmp_dir, avg_climate_layers_file_name)

avg_climate_layers <- try(terra::rast(avg_climate_layers_file),silent = TRUE)

if(class(avg_climate_layers)=="try-error"){
    avg_climate_layers <- terra::app(
        climate_layers, 
        fun = mean, na.rm = TRUE, 
        filename=avg_climate_layers_file,
        overwrite=TRUE,
        cores = ncores)
}

rm(climate_layers);gc()

#######
# Get the spatial gradient (C/km)
cat("Spatial gradient\n")

spgrad_file_name <- paste(variable_name, "spgrad.tif",sep = "_")

spgrad_file <- here::here(tmp_dir,spgrad_file_name)

spgrad <- try(terra::rast(spgrad_file),silent = TRUE)

if(any(class(spgrad)=="try-error")){
    
    spgrad = spatial_grad_big(rx = avg_climate_layers,
                              tmp_dir = here::here(tmp_dir),
                              filename_tiles = paste(variable_name, "spgrad_tile.tif",sep = "_"),
                              filename_final = spgrad_file,
                              ncores = ncores)
    
}

# force raster to pair
ttrend <- terra::project(ttrend, spgrad, threads=TRUE, use_gdal=TRUE, gdal=TRUE)

#######
## calculate gradient-based climate velocity:
cat("Calculate Velocities\n")

## Unprojected
cat("Velocity Unprojected\n")
gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)

## Across latitude
cat("Velocity across latitude\n")
gVelLat <- gVel$Vel * cos(deg_to_rad(gVel$Ang))

# change sign of gVelLat if in the south hemisphere to reflect a velocity away of the tropics
SouthCells <- terra::as.data.frame(gVelLat$Vel, xy = TRUE, cell = TRUE) 
SouthCells <- SouthCells %>% filter(y<0)
SouthCells <- SouthCells$cell
if(length(SouthCells)>0){
    gVelLat[SouthCells] <- gVelLat[SouthCells] * -1
}



#######
cat("Save velocity maps\n")
# save velocity maps

# velocity
terra::writeRaster(gVel, 
                   here::here(vel_dir,paste(variable_name, "gVel.tif",sep = "_")),
                   overwrite=TRUE)

# velocity latitude
terra::writeRaster(gVelLat$Vel, 
                   here::here(vel_dir,paste(variable_name, "gVelLat.tif",sep = "_")),
                   overwrite=TRUE)

# trend
terra::writeRaster(ttrend, 
                   here::here(vel_dir,paste(variable_name, "trend.tif",sep = "_")),
                   overwrite=TRUE)

# spatial gradient
terra::writeRaster(spgrad, 
                   here::here(vel_dir,paste(variable_name, "spatgrad.tif",sep = "_")),
                   overwrite=TRUE)



# Calculate elevation velocities for terrestrial study areas
if(Eco == "Ter" & res_raster == "1km"){
    
    
    ########
    # Get the spatial gradient up slope (elevation/km)
    cat("Spatial gradient up slope\n")
    
    spgrad_ele_file_name <- paste(variable_name, "spgrad_ele.tif",sep = "_")
    
    spgrad_ele_file <- here::here(tmp_dir,spgrad_ele_file_name)
    
    spgrad_ele <- try(terra::rast(spgrad_ele_file), silent = TRUE)
    
    spgrad_ele = spatial_grad_big(rx = elevation,
                                  tmp_dir = here::here(tmp_dir),
                                  filename_tiles = paste(variable_name, "spgradEle_tile.tif",sep = "_"),
                                  filename_final = spgrad_ele_file,
                                  ncores = ncores)
    
    
    ########
    # force rasters to pair
    if(!ext(spgrad_ele)==ext(avg_climate_layers)){
        spgrad_ele <- terra::project(spgrad_ele,avg_climate_layers)
    }
    
    # Get the environmental gradient up slope (C/Elev)
    # Divide the environmental gradient (C/km) by the elevation gradient (Elev/km)
    spgrad_ele$Grad_ele <- spgrad$Grad / spgrad_ele$Grad
    
    #######
    ## calculate gradient-based climate velocity:
    cat("Calculate Velocities\n")
    
    ## Across elevation 
    # cat("Velocity across elevation \n")
    gVelEle <- gVelocity(grad = spgrad_ele, slope = ttrend,
                         grad_col = "Grad_ele", truncate = TRUE)
    
    #######
    ## Project to equal area for more accurate statistics
    gVelEle <- terra::project(gVelEle, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    
    
    #######
    cat("Save velocity maps\n")
    # save velocity maps
    
    # velocity elevation
    terra::writeRaster(spgrad, 
                       here::here(vel_dir,paste(variable_name, "gVelEle.tif",sep = "_")),
                       overwrite=TRUE)
    
    # spatial gradient elevation
    terra::writeRaster(spgrad_ele, 
                       here::here(vel_dir,paste(variable_name, "spatgradEle.tif",sep = "_")),
                       overwrite=TRUE)
    
} 


# delete temporary files
unlink(list.files(tmp_dir,pattern = "global",full.names = TRUE))
