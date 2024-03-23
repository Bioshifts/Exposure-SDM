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
polygontogo <- command_args

# polygontogo <- "A169_P1" # Mar # Small
# polygontogo <- "A112_P1" # Mar # South
# polygontogo <- "A43_P1" # Mar # North
# polygontogo <- "A1_P1" # Ter # North
# polygontogo <- "A10_P1" # Ter # Big # North
# polygontogo <- "A31_P1" # Ter # North # Elevation
# polygontogo <- "A67_P1" # Ter # North
# polygontogo <-"A116_P1"# Ter # very Big # North

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
source("R/velocity_functions.R")

# create dir to store results
if(!dir.exists(velocity_SA_dir)){
    dir.create(velocity_SA_dir)
}

# create dir to store temporary files
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


# all_polis <- list.files(here::here(SA_shps_dir),pattern = ".shp")
# all_polis <- gsub(".shp","",all_polis)

# v1 or v2?
if(grepl("A",polygontogo)){
    # Load Bioshifts version
    Bioshifts_DB <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v1))
}
if(grepl("B",polygontogo)){
    # Load Bioshifts version
    Bioshifts_DB <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v2))
    # Create Polygon IDs
    Bioshifts_DB$ID <- paste0("B",Bioshifts_DB$Paper.ID,"_",Bioshifts_DB$Study.Period)
}

# Filter Polygons in Study areas v3
Bioshifts_DB <- filter(Bioshifts_DB, ID == polygontogo)

# Is it Terrestrial or Marine?
# Check then load temperature data
my_test <- if(any(is.na(SA_i$EleExtentm))){ # it is terrestrial if it has elevation data 
    ECO = "Mar"
} else {
    ECO = "Ter"
    # get elevation
    elevation <- terra::rast(here::here(work_dir,"Data/elevation_1km.tif"))
    # get elevation slope
    elevation_slope <- terra::rast(here::here(work_dir,"Data/elevation_slp_1km.tif"))
    
}

# shifts start and end
if(grepl("A",polygontogo)){
    S_start <- round(unique(Bioshifts_DB$START),0)
    S_end <- round(unique(Bioshifts_DB$END),0)
} else {
    S_start <- round(Bioshifts_DB$Start.Year,0)
    S_end <- round(Bioshifts_DB$End.Year,0)
}
S_time <- min(S_start):max(S_end)


# get layers within time period of shift
if(ECO=="Ter"){
    vars_dir <- here::here(vars_dir(ECO),paste0("bio_proj_",my_res,"_SA"),polygontogo)
}
if(ECO=="Mar"){
    vars_dir <- here::here(vars_dir(ECO),"bio_proj_SA",polygontogo)
}
climate_layers <- list.files(vars_dir)
climate_layers_pos <- grepl(paste(S_time,collapse = "|"),climate_layers)
climate_layers <- climate_layers[climate_layers_pos]
climate_layers <- terra::rast(here::here(vars_dir,climate_layers))

########################
## Calculate velocity at the SA

cat(ECO,"\n")

# If terrestrial
# Calculate velocity for LAT and ELE and for each climatic variable
if(ECO=="Ter"){
    
    # Big area test
    area_i <- SA_i$Areakm2
    big <- area_i > 10^4
    
    # # crop elevation to the study area
    # terra::window(elevation) <- ext(SA_i)
    # elevation <- terra::mask(elevation, SA_i)
    # 
    # # crop elevation slp to the study area
    # terra::window(elevation_slope) <- ext(SA_i)
    # elevation_slope <- terra::mask(elevation_slope, SA_i)
    
    gVelSA <- list()
    
    for(v in 1:length(velocity_variables)){ # for each climate variable
        
        velocity_variable <- velocity_variables[v]
        cat(velocity_variable,"\n")
        
        # select the climate variable
        climate_layers_i <- climate_layers[[which(names(climate_layers)==velocity_variable)]]
        
        # crop climate variables to the study area
        terra::window(climate_layers_i) <- ext(SA_i)
        climate_layers_i <- terra::mask(climate_layers_i, SA_i)
        
        # project to equal-area
        climate_layers_i <- terra::project(climate_layers_i,Eckt)
        
        # if study area is two small (< 8 cells), there is no way climate gradients can be calculated.
        # A possible solution is to disaggregate the raster to a finer resolution
        if(ncell(climate_layers_i) < 12){
            # Force smaller resolution
            climate_layers_i <- terra::disagg(
                climate_layers_i, 
                fact = 12/ncell(climate_layers_i),
                method = "bilinear")
        }
        
        # force raster to pair
        # elevation <- terra::project(elevation,climate_layers_i)
        # elevation_slope <- terra::project(elevation_slope,climate_layers_i)
        
        # in CHELSA, temperature is *10
        if(velocity_variable=="mat"){
            climate_layers_i <- climate_layers_i/10
        }
        
        
        #######
        # calculate the trend (C/year)
        cat("Trend\n")
        
        ttrend_file <- here::here(work_dir,"Data",
                                  paste(ECO,velocity_variable,"trend",paste0(polygontogo,".tif"), sep = "_"))
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
        cat("Spatial gradient\n")
        
        spgrad_file <- here::here(work_dir,"Data",
                                  paste(ECO,velocity_variable,"spgrad",paste0(S_time_name,".qs"), sep = "_"))
        spgrad <- try(qs::qread(spgrad_file,nthreads = ncores),silent = TRUE)
        
        if(any(class(spgrad)=="try-error")){
            
            if(big){
                
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
                
            } else {
                spgrad = spatial_grad(avg_climate_layers)
            }
            
            qs::qsave(spgrad, spgrad_file)
        }
        
        
        # #######
        # # Get the spatial gradient up slope (elevation/km)
        # cat("Spatial gradient up slope\n")
        # 
        # spgrad_ele_file <- here::here(tmp_dir,
        #                               paste(ECO,velocity_variable,"spgrad_ele",paste0(polygontogo,".qs"), sep = "_"))
        # spgrad_ele <- try(qs::qread(spgrad_ele_file,nthreads = ncores),silent = TRUE)
        # if(any(class(spgrad_ele)=="try-error")){
        #     
        #     if(big){
        #         
        #         avg_climate_layers_tiles <- avg_climate_layers
        #         # define tile resolution 
        #         nrow(avg_climate_layers_tiles) <- round(nrow(avg_climate_layers)/100,0)
        #         ncol(avg_climate_layers_tiles) <- round(ncol(avg_climate_layers)/100,0)
        #         avg_climate_layers_tiles[] <- 1:ncell(avg_climate_layers_tiles)
        #         # mask raster
        #         SA_i_Eck <- terra::project(SA_i,avg_climate_layers_tiles)
        #         avg_climate_layers_tiles <- terra::mask(avg_climate_layers_tiles, SA_i_Eck)
        #         
        #         # plot(avg_climate_layers_tiles);dev.off()
        #         
        #         parallel::mclapply(terra::cells(avg_climate_layers_tiles), function(x) {
        #             tmp_file <- here::here(tmp_dir,
        #                                    paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(polygontogo,".qs"), sep = "_"))
        #             test <- try(qs::qread(tmp_file),silent = TRUE)
        #             if(any(class(test)=="try-error")){
        #                 tmp_ext <- terra::ext(avg_climate_layers_tiles, cells = x)
        #                 real_cells <- terra::cells(elevation, tmp_ext)
        #                 terra::window(elevation) <- tmp_ext
        #                 tmp_data <- spatial_grad(elevation)
        #                 terra::window(elevation) <- NULL
        #                 # plug in "real" icells
        #                 tmp_data$icell <- real_cells
        #                 qs::qsave(tmp_data, tmp_file)
        #             }
        #         }, mc.cores = ncores)
        #         
        #         spgrad_ele <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
        #             tmp_file <- here::here(tmp_dir,
        #                                    paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(polygontogo,".qs"), sep = "_"))
        #             try(qs::qread(tmp_file),silent = TRUE)
        #         })
        #         
        #         spgrad_ele <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
        #             tmp_file <- here::here(tmp_dir,
        #                                    paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(polygontogo,".qs"), sep = "_"))
        #             qs::qread(tmp_file)
        #         })
        #         
        #         spgrad_ele <- data.frame(data.table::rbindlist(spgrad_ele))
        #         
        #         qs::qsave(spgrad_ele, spgrad_ele_file)
        #         
        #         # delete temporary files
        #         del_files <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
        #             tmp_file <- here::here(tmp_dir,
        #                                    paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(polygontogo,".qs"), sep = "_"))
        #             unlink(tmp_file)
        #         })
        #     } else {
        #         spgrad_ele = spatial_grad(elevation)
        #     }
        # }
        # 
        # 
        # # What is the environmental gradient up slope? (C/Elev)
        # # Divide the environmental gradient (C/km) by the elevation gradient (Elev/km)
        # spgrad_ele$Grad_ele <- spgrad$Grad / spgrad_ele$Grad
        # 
        # # What is the distance upslope needed to cover the same distance in the flat terrain? (C/km up slope)
        # # Divide the environmental gradient by the cosine of the elevation slope 
        # # This is the same as the formula to find the hypotenuse if provided the basis (gradient C/km) and the angle (elevation slope) of the triangle.
        # spgrad_ele$Grad_ele_dist <- spgrad$Grad / cos(deg_to_rad(elevation_slope[spgrad_ele$icell]))
        
        #######
        ## calculate gradient-based climate velocity:
        cat("Calculate Velocities\n")
        
        ## Unprojected
        cat("Velocity Unprojected\n")
        gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
        
        # ## Across elevation 
        # cat("Velocity across elevation (method elevation)\n")
        # gVelEleUp <- gVelocity(grad = spgrad_ele, slope = ttrend, 
        #                        grad_col = "Grad_ele", truncate = TRUE)
        # 
        # cat("Velocity across elevation (method distance)\n")
        # gVelEleDist <- gVelocity(grad = spgrad_ele, slope = ttrend, 
        #                          grad_col = "Grad_ele_dist", truncate = TRUE)
        
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
        # summary velocity over study area
        baseline <- mean(global(avg_climate_layers,mean,na.rm = TRUE)[,1])
        trend.mean <- as.numeric(global(ttrend,mean,na.rm = TRUE)[,1])
        trend.sd <- as.numeric(global(ttrend,sd,na.rm = TRUE)[,1])
        
        v.mean <- terra::global(gVel$Vel, mean, na.rm=TRUE)[,1]
        v.median <- terra::global(gVel$Vel, median, na.rm=TRUE)[,1]
        v.sd <- terra::global(gVel$Vel, sd, na.rm=TRUE)[,1]
        
        v.lat.mean <- terra::global(gVelLat$Vel, mean, na.rm=TRUE)[,1]
        v.lat.median <- terra::global(gVelLat$Vel, median, na.rm=TRUE)[,1]
        v.lat.sd <- terra::global(gVelLat$Vel, sd, na.rm=TRUE)[,1]
        
        # v.ele.mean <- global(gVelEleUp$Vel, mean, na.rm=TRUE)[,1]
        # v.ele.median <- global(gVelEleUp$Vel, median, na.rm=TRUE)[,1]
        # v.ele.sd <- global(gVelEleUp$Vel, sd, na.rm=TRUE)[,1]
        # 
        # v.ele2.mean <- global(gVelEleDist$Vel, mean, na.rm=TRUE)[,1]
        # v.ele2.median <- global(gVelEleDist$Vel, median, na.rm=TRUE)[,1]
        # v.ele2.sd <- global(gVelEleDist$Vel, sd, na.rm=TRUE)[,1]
        
        gVelSA_j <- data.frame(baseline, trend.mean, trend.sd, 
                               v.mean, v.median, v.sd, 
                               v.lat.mean, v.lat.median, v.lat.sd, 
                               # v.ele.mean, v.ele.median, v.ele.sd,
                               # v.ele2.mean, v.ele2.median, v.ele2.sd
                               )
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
        
        # terra::writeRaster(gVelEleUp$Vel, 
        #                    here::here(velocity_SA_dir, 
        #                               paste(polygontogo,velocity_variable,"gVelEleUp.tif",sep="_")),
        #                    overwrite=TRUE) 
        # terra::writeRaster(gVelEleDist$Vel, 
        #                    here::here(velocity_SA_dir, 
        #                               paste(polygontogo,velocity_variable,"gVelEleDist.tif",sep="_")),
        #                    overwrite=TRUE) 
        
    } 
    
    gVelSA <- do.call(cbind,gVelSA)
    gVelSA$ID=polygontogo
    
    #######
    # save velocity results
    write.csv(gVelSA, here::here(velocity_SA_dir, paste0(polygontogo,".csv")), row.names = FALSE)
    
    # delete temporary files
    unlink(list.files(tmp_dir,pattern = polygontogo,full.names = TRUE))
    
    
    
} else {
    # If marine
    # Calculate only for SST
    
    velocity_variable <- "sst"
    
    gVelSA <- data.frame()
    
    # select the climate variable
    climate_layers_i <- climate_layers[[which(names(climate_layers)=="mean")]]
    
    # crop climate variables to the study area
    terra::window(climate_layers_i) <- ext(SA_i)
    climate_layers_i <- terra::mask(climate_layers_i, SA_i)
    
    # project to equal-area
    climate_layers_i <- terra::project(climate_layers_i,Eckt)
    
    # if study area is two small (< 8 cells), there is no way climate gradients can be calculated.
    # A possible solution is to disaggregate the raster to a finer resolution
    if(ncell(climate_layers_i) < 12){
        # Force smaller resolution
        climate_layers_i <- terra::disagg(
            climate_layers_i, 
            fact = 12/ncell(climate_layers_i),
            method = "bilinear")
    } 
    
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


