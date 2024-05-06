# --------------------------------------------------------
# title: "Create global velocity maps"
# author: "Brunno F Oliveira & Bioshifts group"
# --------------------------------------------------------

# Calculates gradient-based climate velocities

# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","raster","sf","rgdal","Hmisc","dplyr", "data.table","parallelly","qs","rnaturalearth")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

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
source("R/velocity_functions.R")

# create dir to store temporary files
tmp_dir <- here(scratch_dir,"Data/vel_global")
if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir, recursive = TRUE)
}
# create dir to store results
vel_dir <- here(work_dir,"Data/Velocity_global")
if(!dir.exists(vel_dir)){
    dir.create(vel_dir, recursive = TRUE)
}

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
ECO <- as.character(paste(command_args[1], collapse = " "))
velocity_variable <- as.character(paste(command_args[2], collapse = " "))

# ECO="Ter"
# velocity_variable="mat"
# velocity_variable="map"
# ECO="Mar"
# velocity_variable="mean"

# set time period
S_time <- 1960:2009
S_time_name <- paste(range(S_time), collapse = "-")

ncores <- parallelly::availableCores()
########################
## Calculate velocity

cat("Calculating for", ECO, velocity_variable, "\n")

# If terrestrial
# Calculate velocity for LAT and ELE and for each climatic variable
if(ECO=="Ter"){
    
    ttrend_file <- here::here(vel_dir,paste(ECO,velocity_variable,"trend",paste0(S_time_name,".tif"), sep = "_"))
    avg_climate_layers_file <- here::here(vel_dir,paste(ECO,velocity_variable,"avg_climate_layers",paste0(S_time_name,".tif"), sep = "_"))
    spgrad_file <- here::here(vel_dir,paste(ECO,velocity_variable,"spgrad",paste0(S_time_name,".tif"), sep = "_"))
    spgrad_ele_file <- here::here(vel_dir,paste(ECO,velocity_variable,"spgrad_ele",paste0(S_time_name,".tif"), sep = "_"))
    
    test_if_exists <- all(sapply(c(ttrend_file,avg_climate_layers_file,spgrad_file,spgrad_ele_file),
                                 file.exists))
    if(!test_if_exists){
        # get elevation
        elevation <- terra::rast(here::here(work_dir,"Data/elevation_1km.tif"))
        
        # select the climate variable
        climate_layers <- list.files(bios_dir(ECO))
        # temporal range
        climate_layers_pos <- grepl(paste(S_time,collapse = "|"),climate_layers)
        climate_layers <- climate_layers[climate_layers_pos]
        # climate layer
        climate_layers <- terra::rast(here::here(bios_dir(ECO),climate_layers))
        climate_layers <- climate_layers[[which(names(climate_layers)==velocity_variable)]]
        
        # force raster to pair
        cat("force raster to pair\n")
        elevation <- terra::project(elevation, climate_layers, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    }
    
    
    #######
    # calculate the trend (C/year)
    cat("calculate the trend\n")
    
    ttrend_file <- here::here(vel_dir,
                              paste(ECO,velocity_variable,"trend",paste0(S_time_name,".tif"), sep = "_"))
    ttrend <- try(terra::rast(ttrend_file),silent = TRUE)
    if(class(ttrend)=="try-error"){
        ttrend = temp_grad(
            climate_layers,
            th = 0.25*nlyr(climate_layers), ## set minimum N obs. to 1/4 time series length
            tempfile = ttrend_file,
            overwrite = TRUE,
            ncores = ncores)
    }
    
    #######
    # Get averaged climate layers
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
    cat("calculate the spatial gradient\n")
    
    spgrad <- try(terra::rast(spgrad_file),silent = TRUE)
    if(any(class(spgrad)=="try-error")){
        
        avg_climate_layers_tiles <- avg_climate_layers
        # define tile resolution 
        nrow(avg_climate_layers_tiles) <- round(nrow(avg_climate_layers)/100,0)
        ncol(avg_climate_layers_tiles) <- round(ncol(avg_climate_layers)/100,0)
        avg_climate_layers_tiles <- terra::makeTiles(
            avg_climate_layers, 
            avg_climate_layers_tiles, 
            filename = here::here(tmp_dir,paste(Eco,velocity_variable,S_time_name,"tile_.tif",sep="_")),
            extend = TRUE, na.rm = TRUE, overwrite=TRUE)
        
        # x=avg_climate_layers_tiles[227]
        test <- parallel::mclapply(avg_climate_layers_tiles, function(x) {
            tmp_file <- gsub("tile","spgrad",x)
            test <- try(rast(tmp_file),silent = TRUE)
            if(any(class(test)=="try-error")){
                tmp_tile <- terra::rast(x)
                tmp_data <- spatial_grad(tmp_tile)
                terra::writeRaster(tmp_data,tmp_file,overwrite = TRUE)
            }
        }, mc.cores = ncores)
        
        test <- sapply(test,class)
        remove_these <- which(test=="try-error")
        if(length(remove_these)>0){
            tiles_good <- avg_climate_layers_tiles[-remove_these]
        } else {
            tiles_good <- avg_climate_layers_tiles
        }
        
        spgrad <- terra::vrt(gsub("tile","spgrad",tiles_good),
                             set_names = TRUE,
                             overwrite = TRUE)
        
        terra::writeRaster(spgrad,
                           spgrad_file, 
                           overwrite = TRUE)
        
        spgrad <- terra::rast(spgrad_file)
        
        # delete temporary files
        unlink(avg_climate_layers_tiles)
        unlink(gsub("tile","spgrad",avg_climate_layers_tiles))
        
    }
    
    #######
    # Get the spatial gradient up slope (elevation/km)
    cat("Get the spatial gradient up slope \n")
    
    spgrad_ele <- try(terra::rast(spgrad_ele_file),silent = TRUE)
    if(any(class(spgrad_ele)=="try-error")){
        
        elevation_tiles <- elevation
        # define tile resolution 
        nrow(elevation_tiles) <- round(nrow(elevation)/100,0)
        ncol(elevation_tiles) <- round(ncol(elevation)/100,0)
        elevation_tiles <- terra::makeTiles(
            elevation, 
            elevation_tiles, 
            filename = here::here(tmp_dir,paste(Eco,velocity_variable,S_time_name,"tile_.tif",sep="_")),
            na.rm = TRUE, overwrite=TRUE)
        
        # x=elevation_tiles[1]
        test <- parallel::mclapply(elevation_tiles, function(x) {
            tmp_file <- gsub("tile","spgrad_ele",x)
            test <- try(rast(tmp_file),silent = TRUE)
            if(any(class(test)=="try-error")){
                tmp_tile <- terra::rast(x)
                tmp_data <- spatial_grad(tmp_tile)
                terra::writeRaster(tmp_data,tmp_file,overwrite = TRUE)
            }
        }, mc.cores = ncores)
        
        test <- sapply(test,class)
        remove_these <- which(test=="try-error")
        if(length(remove_these)>0){
            tiles_good <- elevation_tiles[-remove_these]
        } else {
            tiles_good <- elevation_tiles
        }
        
        spgrad_ele <- terra::vrt(gsub("tile","spgrad_ele",tiles_good), 
                                 set_names = TRUE,
                                 overwrite = TRUE)
        
        terra::writeRaster(spgrad_ele,
                           spgrad_ele_file, 
                           overwrite = TRUE)
        
        spgrad_ele <- terra::rast(spgrad_ele_file)
        
        # delete temporary files
        unlink(elevation_tiles)
        unlink(gsub("tile","spgrad_ele",elevation_tiles))
        
    }
    
    # force pair
    if(!ext(spgrad_ele)==ext(spgrad)){
        spgrad <- terra::project(spgrad,avg_climate_layers)
        spgrad_ele <- terra::project(spgrad_ele,avg_climate_layers)
    }
    # What is the environmental gradient up slope? (C/Elev)
    # Divide the environmental gradient (C/km) by the elevation gradient (Elev/km)
    spgrad_ele$Grad_ele <- spgrad$Grad / spgrad_ele$Grad
    
    #######
    ## calculate gradient-based climate velocity:
    cat("Calculate Velocities\n")
    
    ## Unprojected
    cat("Velocity Unprojected\n")
    gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
    
    ## Across elevation 
    cat("Velocity across elevation\n")
    gVelEle <- gVelocity(grad = spgrad_ele, slope = ttrend, 
                         grad_col = "Grad_ele", truncate = TRUE)
    
    ## Across latitude
    cat("Velocity across latitude\n")
    gVelLat <- gVel$Vel * cos(deg_to_rad(gVel$Ang))
    
    ## change sign of gVelLat if in the south hemisphere to reflect a velocity away of the tropics
    SouthCells <- terra::as.data.frame(gVelLat$Vel, xy = TRUE, cell = TRUE) 
    SouthCells <- SouthCells %>% filter(y<0)
    SouthCells <- SouthCells$cell
    if(length(SouthCells)>0){
        gVelLat[SouthCells] <- gVelLat[SouthCells] * -1
    }
    
    
    #######
    # Save rasters
    cat("Save velocity maps\n")
    
    #######
    ## Project to equal area before saving
    cat("Project to equal area\n")
    
    model_raster_file <- here(vars_dir(ECO),paste("model_raster_ter",my_res,"Eckt.tif",sep="_"))
    if(file.exists(model_raster_file)){
        model_raster <- rast(model_raster_file)
    } else {
        model_raster <- rast(here(vars_dir(ECO),paste0("model_raster_ter_",my_res,".tif")))
        model_raster <- terra::project(model_raster, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
        writeRaster(model_raster, model_raster_file)
    }
    gVel <- terra::project(gVel, model_raster, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    gVelLat <- terra::project(gVelLat, model_raster, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    gVelEle <- terra::project(gVelEle, model_raster, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    
    terra::writeRaster(
        gVel$Vel,
        filename = here::here(vel_dir,paste(ECO,velocity_variable,"gVel",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVelLat,
        filename = here::here(vel_dir,paste(ECO,velocity_variable,"gVelLat",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVelEle$Vel,
        filename = here::here(vel_dir,paste(ECO,velocity_variable,"gVelEle",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVel$Ang,
        filename = here::here(vel_dir,paste(ECO,velocity_variable,"gVelAngle",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    
    
} else {
    
    # If marine
    # Calculate only for SST
    
    ttrend_file <- here::here(vel_dir,paste(ECO,velocity_variable,"trend",paste0(S_time_name,".tif"), sep = "_"))
    avg_climate_layers_file <- here::here(vel_dir,paste(ECO,velocity_variable,"avg_climate_layers",paste0(S_time_name,".tif"), sep = "_"))
    spgrad_file <- here::here(vel_dir,paste(ECO,velocity_variable,"spgrad",paste0(S_time_name,".tif"), sep = "_"))
    
    test_if_exists <- all(sapply(c(ttrend_file,avg_climate_layers_file,spgrad_file),
                                 file.exists))
    if(!test_if_exists){
        # select the climate variable
        climate_layers <- list.files(bios_dir(ECO))
        # temporal range
        climate_layers_pos <- grepl(paste(S_time,collapse = "|"),climate_layers)
        climate_layers <- climate_layers[climate_layers_pos]
        # climate layer
        climate_layers <- terra::rast(here::here(bios_dir(ECO),climate_layers))
        climate_layers <- climate_layers[[which(names(climate_layers)==velocity_variable)]]
    }
    
    #######
    # calculate the trend (C/year)
    cat("calculate the trend\n")
    
    ttrend_file <- here::here(vel_dir,
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
    
    # calculate averaged climate layers
    avg_climate_layers <- try(terra::rast(avg_climate_layers_file),silent = TRUE)
    if(class(avg_climate_layers)=="try-error"){
        avg_climate_layers <- terra::app(
            climate_layers, 
            fun = mean, na.rm = TRUE, 
            filename=avg_climate_layers_file,
            overwrite=TRUE)
    }
    
    spgrad <- try(terra::rast(spgrad_file),silent = TRUE)
    if(any(class(spgrad)=="try-error")){
        spgrad = spatial_grad(avg_climate_layers)
        terra::writeRaster(spgrad, spgrad_file)
        
    }
    
    #######
    ## calculate gradient-based climate velocity:
    ## Unprojected
    cat("Velocity Unprojected\n")
    gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
    
    ## Across latitude
    cat("Velocity across latitude\n")
    gVelLat <- gVel$Vel * cos(deg_to_rad(gVel$Ang))
    
    #######
    ## change sign of gVelLat if in the south hemisphere to reflect a velocity away of the tropics
    SouthCells <- terra::as.data.frame(gVelLat$Vel, xy = TRUE, cell = TRUE) 
    SouthCells <- SouthCells %>% filter(y<0)
    SouthCells <- SouthCells$cell
    if(length(SouthCells)>0){
        gVelLat[SouthCells] <- gVelLat[SouthCells] * -1
    }
    
    #######
    # Save rasters
    cat("Save velocity maps\n")
    
    #######
    ## Project to equal area before saving
    cat("Project to equal area\n")
    
    model_raster_file <- here(vars_dir(ECO),"model_raster_mar_Eckt.tif")
    if(file.exists(model_raster_file)){
        model_raster <- rast(model_raster_file)
    } else {
        model_raster <- rast(here(vars_dir(ECO),"model_raster_mar.tif"))
        model_raster <- terra::project(model_raster, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
        writeRaster(model_raster, model_raster_file)
    }
    gVel <- terra::project(gVel, model_raster, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    gVelLat <- terra::project(gVelLat, model_raster, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    gVelEle <- terra::project(gVelEle, model_raster, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    
    terra::writeRaster(
        gVelLat,
        filename = here::here(vel_dir,paste(ECO,velocity_variable,"gVelLat",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVel$Vel,
        filename = here::here(vel_dir,paste(ECO,velocity_variable,"gVel",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVel$Ang,
        filename = here::here(vel_dir,paste(ECO,velocity_variable,"gVelAngle",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    
}

