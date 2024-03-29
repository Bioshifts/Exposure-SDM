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
if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir)
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
    
    # get elevation slope
    elevation_slope <- terra::rast(here::here(work_dir,"Data/elevation_slp_1km.tif"))
    
    # select the climate variable
    climate_layers <- climate_layers[[which(names(climate_layers)==velocity_variable)]]
    
    # project to equal-area
    climate_layers <- terra::project(climate_layers,Eckt)
    
    # force raster to pair
    elevation <- terra::project(elevation,climate_layers)
    elevation_slope <- terra::project(elevation_slope,climate_layers)
    
    # in CHELSA, temperature is *10
    if(velocity_variable=="mat"){
        climate_layers <- climate_layers/10
    }
    
    #######
    # calculate the trend (C/year)
    cat("calculate the trend\n")
    
    ttrend_file <- here::here(work_dir,"Data",
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
    avg_climate_layers_file <- here::here(work_dir,"Data",
                                          paste(ECO,velocity_variable,"avg_climate_layers",paste0(polygontogo,".tif"), sep = "_"))
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
    
    spgrad_file <- here::here(work_dir,"Data",
                              paste(ECO,velocity_variable,"spgrad",paste0(S_time_name,".qs"), sep = "_"))
    spgrad <- try(qs::qread(spgrad_file,nthreads = ncores),silent = TRUE)
    if(any(class(spgrad)=="try-error")){
        
        avg_climate_layers_tiles <- avg_climate_layers
        # define tile resolution 
        nrow(avg_climate_layers_tiles) <- round(nrow(avg_climate_layers)/100,0)
        ncol(avg_climate_layers_tiles) <- round(ncol(avg_climate_layers)/100,0)
        avg_climate_layers_tiles[] <- 1:ncell(avg_climate_layers_tiles)
        # mask raster
        SA_i_Eck <- terra::project(SA_i,avg_climate_layers_tiles)
        avg_climate_layers_tiles <- terra::mask(avg_climate_layers_tiles, SA_i_Eck)
        
        # plot(avg_climate_layers_tiles);dev.off()
        
        parallel::mclapply(terra::cells(avg_climate_layers_tiles), function(x) {
            tmp_file <- here::here(tmp_dir,
                                   paste(ECO,velocity_variable,"spgrad_tile",x,paste0(polygontogo,".qs"), sep = "_"))
            test <- try(qs::qread(tmp_file),silent = TRUE)
            if(any(class(test)=="try-error")){
                tmp_ext <- ext(avg_climate_layers_tiles, x)
                terra::window(avg_climate_layers) <- tmp_ext
                tmp_data <- spatial_grad(avg_climate_layers)
                terra::window(avg_climate_layers) <- NULL
                # plug in "real" icells
                real_cells <- terra::cells(avg_climate_layers, tmp_ext)
                tmp_data$icell <- real_cells
                qs::qsave(tmp_data, tmp_file)
            }
        }, mc.cores = ncores)
        
        spgrad <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
            tmp_file <- here::here(tmp_dir,
                                   paste(ECO,velocity_variable,"spgrad_tile",x,paste0(polygontogo,".qs"), sep = "_"))
            qs::qread(tmp_file)
        })
        
        spgrad <- data.frame(data.table::rbindlist(spgrad))
        
        # delete temporary files
        del_files <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
            tmp_file <- here::here(tmp_dir,
                                   paste(ECO,velocity_variable,"spgrad_tile",x,paste0(S_time_name,".qs"), sep = "_"))
            unlink(tmp_file)
        })
        
        qs::qsave(spgrad, spgrad_file)
        
    }
    
    #######
    # Get the spatial gradient up slope (elevation/km)
    cat("Get the spatial gradient up slope \n")
    
    spgrad_ele_file <- here::here(work_dir,"Data",
                                  paste(ECO,velocity_variable,"spgrad_ele",paste0(S_time_name,".qs"), sep = "_"))
    spgrad_ele <- try(qs::qread(spgrad_ele_file,nthreads = ncores),silent = TRUE)
    if(any(class(spgrad_ele)=="try-error")){
        
        # force rasters to pair
        elevation <- terra::project(elevation,avg_climate_layers)
        elevation_slope <- terra::project(elevation_slope,avg_climate_layers)
        # remove oceans
        elevation[is.na(avg_climate_layers)] <- NA
        elevation_slope[is.na(avg_climate_layers)] <- NA
        
        # remove Antarctica
        avg_climate_layers_tiles <- avg_climate_layers
        ext_crop <- ext(avg_climate_layers_tiles)
        ext_crop[3] <- -62
        avg_climate_layers_tiles <- crop(avg_climate_layers_tiles,ext_crop)
        # define tile resolution 
        nrow(avg_climate_layers_tiles) <- round(nrow(avg_climate_layers)/100,0)
        ncol(avg_climate_layers_tiles) <- round(ncol(avg_climate_layers)/100,0)
        avg_climate_layers_tiles[] <- 1:ncell(avg_climate_layers_tiles)
        # mask raster
        land <- vect(rnaturalearthdata::coastline110)
        land <- terra::project(land,avg_climate_layers_tiles)
        countries <- vect(rnaturalearthdata::countries110) # this is needed for masking
        countries <- terra::project(countries,avg_climate_layers_tiles)
        avg_climate_layers_tiles <- terra::mask(avg_climate_layers_tiles, countries)
        
        parallel::mclapply(terra::cells(avg_climate_layers_tiles), function(x) {
            tmp_file <- here::here(tmp_dir,
                                   paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(S_time_name,".qs"), sep = "_"))
            test <- try(qs::qread(tmp_file),silent = TRUE)
            if(any(class(test)=="try-error")){
                tmp_ext <- ext(avg_climate_layers_tiles, x)
                terra::window(elevation) <- tmp_ext
                tmp_data <- spatial_grad(elevation)
                terra::window(elevation) <- NULL
                # plug in "real" icells
                real_cells <- terra::cells(elevation, tmp_ext)
                tmp_data$icell <- real_cells
                qs::qsave(tmp_data, tmp_file)
            }
        }, mc.cores = ncores)
        
        spgrad_ele <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
            tmp_file <- here::here(tmp_dir,
                                   paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(S_time_name,".qs"), sep = "_"))
            qs::qread(tmp_file)
        })
        
        spgrad_ele <- data.frame(data.table::rbindlist(spgrad_ele))
        
        qs::qsave(spgrad_ele, spgrad_ele_file)
        
        # delete temporary files
        del_files <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
            tmp_file <- here::here(tmp_dir,
                                   paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(S_time_name,".qs"), sep = "_"))
            unlink(tmp_file)
        })
    }
    
    # What is the environmental gradient up slope? (C/Elev)
    # Divide the environmental gradient (C/km) by the elevation gradient (Elev/km)
    spgrad_ele$Grad_ele <- spgrad$Grad / spgrad_ele$Grad
    
    # What is the distance upslope needed to cover the same distance in the flat terrain? (C/km up slope)
    # Divide the environmental gradient by the cosine of the elevation slope 
    # This is the same as the formula to find the hypotenuse if provided the basis (gradient C/km) and the angle (elevation slope) of the triangle.
    spgrad_ele$Grad_ele_dist <- spgrad$Grad / cos(deg_to_rad(elevation_slope[spgrad_ele$icell]))
    
    #######
    ## calculate gradient-based climate velocity:
    cat("Calculate Velocities\n")
    
    ## Unprojected
    cat("Velocity Unprojected\n")
    gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
    
    ## Across elevation 
    cat("Velocity across elevation (method elevation)\n")
    gVelEleUp <- gVelocity(grad = spgrad_ele, slope = ttrend, 
                           grad_col = "Grad_ele", truncate = TRUE)
    
    cat("Velocity across elevation (method distance)\n")
    gVelEleDist <- gVelocity(grad = spgrad_ele, slope = ttrend, 
                             grad_col = "Grad_ele_dist", truncate = TRUE)
    
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
    
    
    # in CHELSA, temperature is *10
    if(velocity_variable=="mat"){
        gVel$Vel <- gVel$Vel/10
        gVelEleUp$Vel <- gVelEleUp$Vel/10
        gVelEleDist$Vel <- gVelEleDist$Vel/10
        gVelLat <- gVelLat/10
    }
    
    #######
    # Save rasters
    cat("Save velocity maps\n")
    
    terra::writeRaster(
        gVelLat,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelLat",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVelEleUp$Vel,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelEleUp",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVelEleDist$Vel,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelEleDist",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVel$Vel,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVel",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVel$Ang,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelAngle",paste0(S_time_name,".tif"), sep = "_")),
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
    
    # calculate averaged climate layers
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
    
    spgrad_file <- here::here(work_dir,"Data",
                              paste(ECO,velocity_variable,"spgrad",paste0(S_time_name,".qs"), sep = "_"))
    spgrad <- try(qs::qread(spgrad_file,nthreads = ncores),silent = TRUE)
    if(any(class(spgrad)=="try-error")){
        spgrad <- spatial_grad(avg_climate_layers)
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
    terra::writeRaster(
        gVelLat,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelLat",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVel$Vel,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVel",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    terra::writeRaster(
        gVel$Ang,
        filename = here::here("Data",paste(ECO,velocity_variable,"gVelAngle",paste0(S_time_name,".tif"), sep = "_")),
        overwrite = TRUE)
    
}

