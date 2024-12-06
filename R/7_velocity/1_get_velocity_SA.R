# --------------------------------------------------------
# title: "Calculate climate change velocity for each study area in the Bioshifts database"
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
# Args
command_args <- commandArgs(trailingOnly = TRUE)
print(command_args)
polygontogo <- as.character(paste(command_args[1], collapse = " "))
Eco <- as.character(paste(command_args[2], collapse = " "))
res_raster <- as.character(paste(command_args[3], collapse = " "))

# Eco <- "Mar"
# Eco <- "Ter"
# res_raster <- "25km"
# res_raster <- "1km"
# polygontogo <- "B564_P2" # Mar # Mediterranean Sea
# polygontogo <- "A146_P1" # Mar # Southwest Pacific Ocean
# polygontogo <- "A138_P1" # Ter # Scotland
# polygontogo <- "A136_P3" # Ter # North America
# polygontogo <- "A31_P1" # Ter # Southern Scandes, Norway
# polygontogo <- "A34_P1" # Ter # Bavarian Forest National Park
# polygontogo <- "A25_P1" # Ter # Asturias region, Cantabrian Range
# polygontogo <-"A116_P1" # Ter # Worldwide, Morth hemisphere
# polygontogo <- "A79_P1" # Ter # Rhone-Saone Valley; Southeastern France
# polygontogo <- "A134_P1" # Ter # United Kingdom
# polygontogo <- "A1_P1" # Ter # French Alps (giffre Valley)

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

if(Eco == "Ter"){
    my_res = res_raster
} else {
    my_res = "25km"
}


# create dir to store results
if(!dir.exists(velocity_SA_dir)){
    dir.create(velocity_SA_dir)
}

# create dir to store temporary files
tmp_dir <- here(tmp_dir,"vel_SA")
if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir, recursive = TRUE)
}

# velocity variables
if(Eco == "Ter"){
    velocity_variables <- c("mat","map")
} else {
    velocity_variables <- c("sst")
}

# N cores
ncores <- parallelly::availableCores()

########################
# get period of interest
bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_fix_columns(bioshifts)
bioshifts <- bioshifts %>% filter(ID == polygontogo)

period <- bioshifts %>% 
    dplyr::select(Start,End) %>%
    round(digits = 0) %>%
    unique 
period = min(period$Start):max(period$End)

########################
# select environmental variable for period
climate_layers <- list.files(bios_dir(Eco), full.names = TRUE)
climate_layers <- climate_layers[grep(paste0(period,collapse = "|"),climate_layers)]
climate_layers <- rast(climate_layers)

# select the environmental variable of interest
if(Eco == "Ter"){
    climate_layers <- climate_layers[[which(names(climate_layers) %in% velocity_variables)]]
} else {
    climate_layers <- climate_layers["mean"]
    names(climate_layers) <- rep(velocity_variables, nlyr(climate_layers))
}

########################
# load SA polygon
SA_i <- terra::vect(here::here(SA_shps_dir,paste0(polygontogo,".shp")))
# add buffer to SA_i to avoid edge effect when calculating velocities
SA_i <- terra::buffer(SA_i, width = res(climate_layers)[1])
# plot(SA_i);dev.off()

########################
# mask environmental variables to the SA
terra::window(climate_layers) <- terra::ext(SA_i)
climate_layers <- terra::mask(climate_layers, mask = SA_i)
terra::window(climate_layers) <- NULL
# plot(climate_layers[[1]]);dev.off()

########################
# Get elevation data if running for terrestrial environment
if(Eco == "Ter"){
    # get elevation
    elevation <- terra::rast(here::here(work_dir,"Data/elevation_1km.tif"))
    # crop elevation to the study area
    terra::window(elevation) <- ext(SA_i)
    elevation <- terra::mask(elevation, SA_i)
    # force raster to pair
    elevation <- terra::project(elevation, climate_layers, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
}

########################
# Big area test - Code will run in parallel for big areas
if(Eco=="Ter" & my_res=="1km"){
    # Big area test
    area_i <- SA_i$Areakm2
    big_raster <- area_i > 10^4
} else {
    big_raster <- FALSE
}

# writeRaster(climate_layers, "A138_P1_climate.tif")
# writeRaster(elevation, "A138_P1_elevation.tif")

########################
## Calculate velocity at the SA

cat(Eco,"\n")

# If terrestrial
# Calculate velocity for LAT and ELE and for each climatic variable
gVelSA <- list()

for(v in 1:length(velocity_variables)){ # for each climate variable
    
    velocity_variable <- velocity_variables[v]
    
    cat(velocity_variable,"\n")
    
    # select the climate variable
    climate_layers_i <- climate_layers[velocity_variable]
    
    #######
    # calculate the trend (C/year)
    cat("Trend\n")
    
    ttrend_file_name <- paste(Eco,velocity_variable,"trend",paste0(polygontogo,".tif"), sep = "_")
    
    ttrend_file <- here::here(tmp_dir,ttrend_file_name)
    
    ttrend <- try(terra::rast(ttrend_file),silent = TRUE)
    
    if(class(ttrend)=="try-error"){
        ttrend = temp_grad(
            climate_layers_i,
            th = 0.25*nlyr(climate_layers_i), ## set minimum N obs. to 1/4 time series length
            file_name = ttrend_file,
            overwrite = TRUE,
            ncores = switch(big_raster + 0, ncores, NULL))
    }
    
    #######
    # Get averaged climate layers
    avg_climate_layers_file_name <- paste(Eco,velocity_variable,"avg_climate_layers",paste0(polygontogo,".tif"), sep = "_")
    
    avg_climate_layers_file <- here::here(tmp_dir, avg_climate_layers_file_name)
    
    avg_climate_layers <- try(terra::rast(avg_climate_layers_file),silent = TRUE)
    
    if(class(avg_climate_layers)=="try-error"){
        avg_climate_layers <- terra::app(
            climate_layers_i, 
            fun = mean, na.rm = TRUE, 
            filename=avg_climate_layers_file,
            overwrite=TRUE,
            cores = switch(big_raster + 0, ncores, NULL))
    }
    
    rm(climate_layers_i);gc()
    
    #######
    # Get the spatial gradient (C/km)
    cat("Spatial gradient\n")
    
    spgrad_file_name <- paste(Eco,velocity_variable,"spgrad",paste0(polygontogo,".tif"), sep = "_")
    
    spgrad_file <- here::here(tmp_dir,spgrad_file_name)
    
    spgrad <- try(terra::rast(spgrad_file),silent = TRUE)
    
    if(any(class(spgrad)=="try-error")){
        
        if(big_raster){
            
            spgrad = spatial_grad_big(rx = avg_climate_layers,
                                      tmp_dir = here::here(tmp_dir),
                                      filename_tiles = paste(Eco,velocity_variable,polygontogo,"tile_.tif",sep="_"),
                                      filename_final = spgrad_file,
                                      ncores = ncores)
            
        } else {
            spgrad = spatial_grad(avg_climate_layers)
            
            terra::writeRaster(spgrad,
                               spgrad_file, 
                               overwrite = TRUE)
            
            spgrad <- terra::rast(spgrad_file)
        }
    }
    
    
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
    ## Project to equal area for more accurate statistics
    gVel <- terra::project(gVel, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    gVelLat <- terra::project(gVelLat, gVel, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    
    # plot(gVelLat);dev.off()
    
    #######
    # summary velocity over study area
    
    baseline <- as.numeric(mean(terra::global(avg_climate_layers,mean,na.rm = TRUE)[,1])) 
    trend.mean <- as.numeric(terra::global(ttrend,mean,na.rm = TRUE)[,1])
    trend.sd <- as.numeric(terra::global(ttrend,sd,na.rm = TRUE)[,1])
    
    # in CHELSA, temperature is *10  
    if(velocity_variable == "mat"){
        baseline <- baseline/10
        trend.mean <- trend.mean/10
        trend.sd <- trend.sd/10
    } 
    
    v.mean <- as.numeric(terra::global(gVel$Vel, mean, na.rm=TRUE)[,1])
    v.median <- as.numeric(terra::global(gVel$Vel, median, na.rm=TRUE)[,1])
    v.sd <- as.numeric(terra::global(gVel$Vel, sd, na.rm=TRUE)[,1])
    
    v.lat.mean <- as.numeric(terra::global(gVelLat$Vel, mean, na.rm=TRUE)[,1])
    v.lat.median <- as.numeric(terra::global(gVelLat$Vel, median, na.rm=TRUE)[,1])
    v.lat.sd <- as.numeric(terra::global(gVelLat$Vel, sd, na.rm=TRUE)[,1])
    
    gVelSA_j <- data.frame(baseline, trend.mean, trend.sd, 
                           v.mean, v.median, v.sd,
                           v.lat.mean, v.lat.median, v.lat.sd)
    names(gVelSA_j) <- paste0(names(gVelSA_j),".",velocity_variable)
    
    gVelSA[[v]] <- gVelSA_j
    
    #######
    cat("Save velocity maps\n")
    # save velocity maps
    
    # velocity
    varname <- paste(polygontogo,velocity_variable,"gVel",my_res,sep="_")
    terra::writeRaster(gVel, 
                       here::here(velocity_SA_dir, paste0(varname,".tif")),
                       overwrite=TRUE)
    
    # velocity latitude
    varname <- paste(polygontogo,velocity_variable,"gVelLat",my_res,sep="_")
    terra::writeRaster(gVelLat$Vel, 
                       here::here(velocity_SA_dir, paste0(varname,".tif")),
                       overwrite=TRUE)
    
    # trend
    varname <- paste(polygontogo,velocity_variable,"trend",my_res,sep="_")
    terra::writeRaster(ttrend, 
                       here::here(velocity_SA_dir, paste0(varname,".tif")),
                       overwrite=TRUE)
    
    # spatial gradient
    varname <- paste(polygontogo,velocity_variable,"spatgrad",my_res,sep="_")
    terra::writeRaster(spgrad, 
                       here::here(velocity_SA_dir, paste0(varname,".tif")),
                       overwrite=TRUE)
} 

gVelSA <- do.call(cbind,gVelSA)
gVelSA <- data.frame(ID=polygontogo,gVelSA)

# Calculate elevation velocities for terrestrial study areas
if(Eco == "Ter" & res_raster == "1km"){
    
    gVelSA_ele <- list()
    
    for(v in 1:length(velocity_variables)){ # for each climate variable
        
        velocity_variable <- velocity_variables[v]
        
        cat(velocity_variable,"\n")
        
        ########
        # Get the spatial gradient up slope (elevation/km)
        cat("Spatial gradient up slope\n")
        
        spgrad_ele_file_name <- paste(Eco,velocity_variable,"spgrad_ele",paste0(polygontogo,".tif"), sep = "_")
        
        spgrad_ele_file <- here::here(tmp_dir,spgrad_ele_file_name)
        
        spgrad_ele <- try(terra::rast(spgrad_ele_file), silent = TRUE)
        
        if(big_raster){
            
            spgrad_ele = spatial_grad_big(rx = elevation,
                                          tmp_dir = here::here(tmp_dir),
                                          filename_tiles = paste(Eco,velocity_variable,polygontogo,"ele_tile_.tif",sep="_"),
                                          filename_final = spgrad_ele_file,
                                          ncores = ncores)
            
        } else {
            spgrad_ele = spatial_grad(elevation)
            
            terra::writeRaster(spgrad_ele, spgrad_ele_file, overwrite = TRUE)
            
            spgrad_ele <- terra::rast(spgrad_ele_file)
        }
        
        
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
        # summary velocity over study area
        
        v.ele.mean <- as.numeric(terra::global(gVelEle$Vel, mean, na.rm=TRUE)[,1])
        v.ele.median <- as.numeric(terra::global(gVelEle$Vel, median, na.rm=TRUE)[,1])
        v.ele.sd <- as.numeric(terra::global(gVelEle$Vel, sd, na.rm=TRUE)[,1])
        
        gVelSA_ele_j <- data.frame(v.ele.mean, v.ele.median, v.ele.sd)
        names(gVelSA_ele_j) <- paste0(names(gVelSA_ele_j),".",velocity_variable)
        
        gVelSA_ele[[v]] <- gVelSA_ele_j
        
        
        #######
        cat("Save velocity maps\n")
        # save velocity maps
        
        # velocity elevation
        varname <- paste(polygontogo,velocity_variable,"gVelEle",my_res,sep="_")
        terra::writeRaster(gVelEle,
                           here::here(velocity_SA_dir, paste0(varname,".tif")),
                           overwrite=TRUE)
        
        # spatial gradient elevation
        varname <- paste(polygontogo,velocity_variable,"spatgradEle",my_res,sep="_")
        terra::writeRaster(spgrad_ele, 
                           here::here(velocity_SA_dir, paste0(varname,".tif")),
                           overwrite=TRUE)
        
    } 
    
    gVelSA_ele <- do.call(cbind,gVelSA_ele)
    
    gVelSA <- cbind(gVelSA, gVelSA_ele)
}


#######
# save velocity results
write.csv(gVelSA, here::here(velocity_SA_dir, paste0(polygontogo,"_",my_res,".csv")), row.names = FALSE)

# delete temporary files
unlink(list.files(tmp_dir,pattern = polygontogo,full.names = TRUE))
