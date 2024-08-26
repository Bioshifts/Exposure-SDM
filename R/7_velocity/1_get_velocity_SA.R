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

# Eco <- "Ter"
# res_raster <- "25km"
# polygontogo <- "A112_P1" # Mar # South
# polygontogo <- "A43_P1" # Mar # North
# polygontogo <- "A1_P1" # Ter # North
# polygontogo <- "A10_P1" # Ter # Big # North
# polygontogo <- "A31_P1" # Ter # North # Elevation
# polygontogo <- "A67_P1" # Ter # North
# polygontogo <- "A25_P1" # Ter # North # Not big
# polygontogo <-"A116_P1"# Ter # very Big # North
# polygontogo <- "A79_P1"

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

my_res <- res_raster
# create dir to store results
if(!dir.exists(velocity_SA_dir)){
    dir.create(velocity_SA_dir)
}

# create dir to store temporary files
tmp_dir <- here(tmp_dir,"vel_SA")
if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir)
}

# velocity variables
velocity_variables <- c("mat","map")

# gradients
gradients <- c("LAT","ELE")

# N cores
ncores <- parallelly::availableCores()

########################
# load SA polygon
SA_i <- terra::vect(here::here(SA_shps_dir,paste0(polygontogo,".shp")))
# plot(SA_i);dev.off()

if(Eco == "Ter"){
    # get elevation
    elevation <- terra::rast(here::here(work_dir,"Data/elevation_1km.tif"))
}

# get layers for SA
climate_layers <- list.files(here(bios_SA_dir(Eco),polygontogo),full.names = TRUE)
climate_layers <- terra::rast(climate_layers)
# plot(climate_layers[[1]]);dev.off()

########################
## Calculate velocity at the SA

cat(Eco,"\n")

# If terrestrial
# Calculate velocity for LAT and ELE and for each climatic variable
if(Eco=="Ter" & my_res=="1km"){
    
    # Big area test
    area_i <- SA_i$Areakm2
    big <- area_i > 10^4
    
    # crop elevation to the study area
    terra::window(elevation) <- ext(SA_i)
    elevation <- terra::mask(elevation, SA_i)
    
    gVelSA <- list()
    
    for(v in 1:length(velocity_variables)){ # for each climate variable
        
        velocity_variable <- velocity_variables[v]
        
        cat(velocity_variable,"\n")
        
        # select the climate variable
        climate_layers_i <- climate_layers[[which(names(climate_layers)==velocity_variable)]]
        # force raster to pair
        elevation <- terra::project(elevation, climate_layers_i, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
        
        #######
        # calculate the trend (C/year)
        cat("Trend\n")
        
        ttrend_file <- here::here(tmp_dir,
                                  paste(Eco,velocity_variable,"trend",paste0(polygontogo,".tif"), sep = "_"))
        ttrend <- try(terra::rast(ttrend_file),silent = TRUE)
        if(class(ttrend)=="try-error"){
            ttrend = temp_grad(
                climate_layers_i,
                th = 0.25*nlyr(climate_layers_i), ## set minimum N obs. to 1/4 time series length
                tempfile = ttrend_file,
                overwrite = TRUE,
                ncores = ncores)
        }
        
        #######
        # Get averaged climate layers
        avg_climate_layers_file <- here::here(tmp_dir,
                                              paste(Eco,velocity_variable,"avg_climate_layers",paste0(polygontogo,".tif"), sep = "_"))
        avg_climate_layers <- try(terra::rast(avg_climate_layers_file),silent = TRUE)
        if(class(avg_climate_layers)=="try-error"){
            avg_climate_layers <- terra::app(
                climate_layers_i, 
                fun = mean, na.rm = TRUE, 
                filename=avg_climate_layers_file,
                overwrite=TRUE,
                cores = ncores)
        }
        
        rm(climate_layers_i);gc()
        
        #######
        # Get the spatial gradient (C/km)
        cat("Spatial gradient\n")
        
        spgrad_file <- here::here(tmp_dir,
                                  paste(Eco,velocity_variable,"spgrad",paste0(polygontogo,".tif"), sep = "_"))
        spgrad <- try(terra::rast(spgrad_file),silent = TRUE)
        if(any(class(spgrad)=="try-error")){
            
            if(big){
                
                avg_climate_layers_tiles <- avg_climate_layers
                # define tile resolution 
                nrow(avg_climate_layers_tiles) <- round(nrow(avg_climate_layers)/1000,0)
                ncol(avg_climate_layers_tiles) <- round(ncol(avg_climate_layers)/1000,0)
                # avg_climate_layers_tiles[] <- 1:ncell(avg_climate_layers_tiles)
                # plot(avg_climate_layers_tiles);dev.off()
                avg_climate_layers_tiles <- terra::makeTiles(
                    avg_climate_layers, 
                    avg_climate_layers_tiles, 
                    filename = here::here(tmp_dir,paste(Eco,velocity_variable,polygontogo,"tile_.tif",sep="_")),
                    na.rm = TRUE, overwrite=TRUE)
                
                # x=avg_climate_layers_tiles[1]
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
                
            } else {
                spgrad = spatial_grad(avg_climate_layers)
                
                terra::writeRaster(spgrad, spgrad_file)
                
                spgrad <- terra::rast(spgrad_file)
            }
        }
        
        
        # #######
        # Get the spatial gradient up slope (elevation/km)
        cat("Spatial gradient up slope\n")
        
        spgrad_ele_file <- here::here(tmp_dir,
                                      paste(Eco,velocity_variable,"spgrad_ele",paste0(polygontogo,".tif"), sep = "_"))
        spgrad_ele <- try(expr = {terra::rast(spgrad_ele_file)}, silent = TRUE)
        if(any(class(spgrad_ele)=="try-error")){
            
            if(big){
                
                elevation_tiles <- elevation
                # define tile resolution 
                nrow(elevation_tiles) <- round(nrow(elevation)/100,0)
                ncol(elevation_tiles) <- round(ncol(elevation)/100,0)
                elevation_tiles <- terra::makeTiles(
                    elevation, 
                    elevation_tiles, 
                    filename = here::here(tmp_dir,paste(Eco,velocity_variable,polygontogo,"tile_.tif",sep="_")),
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
                
            } else {
                spgrad_ele = spatial_grad(elevation)
                
                terra::writeRaster(spgrad_ele, spgrad_ele_file)
                
                spgrad_ele <- terra::rast(spgrad_ele_file)
                
            }
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
        
        # ## Across elevation 
        # cat("Velocity across elevation \n")
        gVelEle <- gVelocity(grad = spgrad_ele, slope = ttrend,
                             grad_col = "Grad_ele", truncate = TRUE)
        
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
        gVelEle <- terra::project(gVelEle, gVel, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
        
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
        
        v.ele.mean <- s.numeric(terra::global(gVelEle$Vel, mean, na.rm=TRUE)[,1])
        v.ele.median <- s.numeric(terra::global(gVelEle$Vel, median, na.rm=TRUE)[,1])
        v.ele.sd <- s.numeric(terra::global(gVelEle$Vel, sd, na.rm=TRUE)[,1])
        
        gVelSA_j <- data.frame(baseline, trend.mean, trend.sd, 
                               v.mean, v.median, v.sd, 
                               v.lat.mean, v.lat.median, v.lat.sd, 
                               v.ele.mean, v.ele.median, v.ele.sd)
        names(gVelSA_j) <- paste0(names(gVelSA_j),".",velocity_variable)
        
        gVelSA[[v]] <- gVelSA_j
        
        #######
        cat("Save velocity maps\n")
        # save velocity maps
        terra::writeRaster(gVel$Vel, 
                           here::here(velocity_SA_dir, 
                                      paste(polygontogo,velocity_variable,"gVel.tif",sep="_")),
                           overwrite=TRUE)
        terra::writeRaster(gVelLat$Vel, 
                           here::here(velocity_SA_dir, 
                                      paste(polygontogo,velocity_variable,"gVelLat.tif",sep="_")),
                           overwrite=TRUE)
        terra::writeRaster(gVelEle$Vel,
                           here::here(velocity_SA_dir,
                                      paste(polygontogo,velocity_variable,"gVelEle.tif",sep="_")),
                           overwrite=TRUE)
        
    } 
    
    gVelSA <- do.call(cbind,gVelSA)
    gVelSA$ID=polygontogo
    
    #######
    # save velocity results
    write.csv(gVelSA, here::here(velocity_SA_dir, paste0(polygontogo,".csv")), row.names = FALSE)
    
    # delete temporary files
    unlink(list.files(tmp_dir,pattern = polygontogo,full.names = TRUE))
    
    
} 
if(Eco=="Mar"){
    
    # If marine
    # Calculate only for SST
    velocity_variable <- "sst"
    
    # select the climate variable
    climate_layers_i <- climate_layers[[which(names(climate_layers)=="mean")]]
    
    #######
    # calculate the trend (C/year)
    cat("Trend\n")
    ttrend = temp_grad(climate_layers_i,
                       th = 0.25*nlyr(climate_layers_i) ## set minimum # obs. to 1/4 time series length
    )
    
    #######
    # calculate the spatial gradient (C/km)
    cat("Spatial gradient\n")
    
    ## calculate averaged climate layers
    avg_climate_layers <- terra::app(
        climate_layers_i, 
        fun = mean, na.rm = TRUE, 
        overwrite=TRUE,
        cores = ncores)
    
    gc()
    
    rm(climate_layers_i)
    
    spgrad = spatial_grad(avg_climate_layers)
    
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
    ## Project to equal area for more accurate statistics
    gVel <- terra::project(gVel, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    gVelLat <- terra::project(gVelLat, gVel, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    
    #######
    # summary velocity over study area
    baseline <- mean(global(avg_climate_layers,mean,na.rm = TRUE)[,1])
    trend.mean <- as.numeric(global(ttrend,mean,na.rm = TRUE)[,1])
    trend.sd <- as.numeric(global(ttrend,sd,na.rm = TRUE)[,1])
    
    v.mean <- terra::global(gVel$Vel, mean, na.rm=TRUE)[,1]
    v.median <- terra::global(gVel$Vel, median, na.rm=TRUE)[,1]
    v.sd <- terra::global(gVel$Vel, sd, na.rm=TRUE)[,1]
    
    v.lat.mean <- as.numeric(global(gVelLat$Vel, mean, na.rm=TRUE))
    v.lat.median <- as.numeric(global(gVelLat$Vel, median, na.rm=TRUE))
    v.lat.sd <- as.numeric(global(gVelLat$Vel, sd, na.rm=TRUE))
    
    gVelSA <- data.frame(baseline, trend.mean, trend.sd, 
                         v.mean, v.median, v.sd, 
                         v.lat.mean, v.lat.median, v.lat.sd)
    names(gVelSA) <- paste0(names(gVelSA),".",velocity_variable)
    gVelSA$ID=polygontogo
    
    #######
    # save velocity results
    write.csv(gVelSA, here::here(velocity_SA_dir, paste0(polygontogo,".csv")), row.names = FALSE)
    
    #######
    # save velocity maps
    cat("Save velocity maps\n")
    terra::writeRaster(gVel$Vel, 
                       here::here(velocity_SA_dir, 
                                  paste(polygontogo,velocity_variable,"gVel.tif",sep="_")),
                       overwrite=TRUE)
    terra::writeRaster(gVelLat, 
                       here::here(velocity_SA_dir, 
                                  paste(polygontogo,velocity_variable,"gVelLat.tif",sep="_")),
                       overwrite=TRUE)
    
}


if((Eco=="Ter" & !my_res=="1km")){
    
    gVelSA <- list()
    
    for(v in 1:length(velocity_variables)){ # for each climate variable
        
        velocity_variable <- velocity_variables[v]
        
        cat(velocity_variable,"\n")
        
        # select the climate variable
        climate_layers_i <- climate_layers[[which(names(climate_layers)==velocity_variable)]]
        
        
        #######
        # calculate the trend (C/year)
        cat("Trend\n")
        ttrend = temp_grad(climate_layers_i,
                           th = 0.25*nlyr(climate_layers_i) ## set minimum # obs. to 1/4 time series length
        )
        
        #######
        # calculate the spatial gradient (C/km)
        cat("Spatial gradient\n")
        
        ## calculate averaged climate layers
        avg_climate_layers <- terra::app(
            climate_layers_i, 
            fun = mean, na.rm = TRUE, 
            overwrite=TRUE,
            cores = ncores)
        
        gc()
        
        rm(climate_layers_i)
        
        spgrad = spatial_grad(avg_climate_layers)
        
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
        ## Project to equal area for more accurate statistics
        gVel <- terra::project(gVel, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
        gVelLat <- terra::project(gVelLat, gVel, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
        
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
        terra::writeRaster(gVel$Vel, 
                           here::here(velocity_SA_dir, 
                                      paste(polygontogo,velocity_variable,my_res,"gVel.tif",sep="_")),
                           overwrite=TRUE)
        terra::writeRaster(gVelLat$Vel, 
                           here::here(velocity_SA_dir, 
                                      paste(polygontogo,velocity_variable,my_res,"gVelLat.tif",sep="_")),
                           overwrite=TRUE)
        
    } 
    
    gVelSA <- do.call(cbind,gVelSA)
    gVelSA$ID=polygontogo
    
    #######
    # save velocity results
    write.csv(gVelSA, here::here(velocity_SA_dir, paste0(paste(polygontogo,my_res,sep = "_"),".csv")), row.names = FALSE)
    
}
