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
# polygontogo <- "A30_P1" # random test


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
    # get elevation layer
    if(!file.exists(here::here(work_dir,"Data/elevation_1km.tif"))){
        download.file("https://data.earthenv.org/topography/elevation_1KMmn_GMTEDmn.tif",
                      here::here(work_dir,"Data/elevation_1km.tif"))
    }
    elevation <- terra::rast(here::here(work_dir,"Data/elevation_1km.tif"))
}

# shifts start and end
if(grepl("A",polygontogo)){
    S_start <- round(unique(Bioshifts_DB$START),0)
    S_end <- round(unique(Bioshifts_DB$END),0)
} else {
    S_start <- round(Bioshifts_DB$Start.Year,0)
    S_end <- round(Bioshifts_DB$End.Year,0)
}
S_time <- S_start[1]:S_end[1]


# get layers within time period of shift
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

########################
## Calculate velocity at the SA

cat(ECO,"\n")

# If terrestrial
# Calculate velocity for LAT and ELE and for each climatic variable
if(ECO=="Ter"){
    
    # crop elevation to the study area
    terra::window(elevation) <- ext(SA_i)
    elevation <- terra::mask(elevation, SA_i)
    
    gVelSA <- list()
    
    for(v in 1:length(velocity_variables)){ # for each climate variable
        
        velocity_variable <- velocity_variables[v]
        cat(velocity_variable,"\n")
        
        # select the climate variable
        climate_layers_i <- climate_layers[[which(names(climate_layers)==velocity_variable)]]
        
        # crop layers to the study area
        terra::window(climate_layers_i) <- ext(SA_i)
        climate_layers_i <- terra::mask(climate_layers_i, SA_i)
        terra::window(climate_layers_i) <- NULL
        
        # project to equal-area
        climate_layers_i <- terra::project(climate_layers_i,Eckt)
        
        # in CHELSA, temperature is *10
        if(velocity_variable=="mat"){
            climate_layers_i <- climate_layers_i/10
        }
        
        # if study area is two small (there are less then 8 cells), there is no way climate gradients can be calculated.
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
        # parallelize for big areas
        cat("Trend\n")
        
        ttrend_file <- 
            here::here(tmp_dir,
                       paste(ECO,velocity_variable,"tmp_trend",paste0(polygontogo,".tif"), sep = "_"))
        ttrend <- try(terra::rast(ttrend_file),silent = TRUE)
        if(class(ttrend)=="try-error"){
            
            if(ncores > 1){
                
                ttrend = temp_grad(
                    climate_layers_i,
                    th = 0.25*nlyr(climate_layers_i), ## set minimum N obs. to 1/4 time series length
                    tempfile = ttrend_file,
                    overwrite = TRUE,
                    ncores = ncores)
                
            } else {
                
                ttrend = temp_grad(climate_layers_i,
                                   th = 0.25*nlyr(climate_layers_i)) ## set minimum # obs. to 1/4 time series length
                
            }
        } 
        
        
        #######
        # Get the spatial gradient (C/km)
        cat("spatial gradient\n")
        
        spgrad_file <- 
            here::here(tmp_dir,
                       paste(ECO,velocity_variable,"spgrad",paste0(polygontogo,".qs"), sep = "_"))
        spgrad <- try(qs::qread(spgrad_file,nthreads = ncores),silent = TRUE)
        if(any(class(spgrad)=="try-error")){
            
            # 1) calculate averaged climate layers
            avg_climate_layers_file <- 
                here::here(tmp_dir,
                           paste(ECO,velocity_variable,"tmp_avg_climate_layers",paste0(polygontogo,".tif"), sep = "_"))
            avg_climate_layers <- try(terra::rast(avg_climate_layers_file),silent = TRUE)
            if(any(class(avg_climate_layers)=="try-error")){
                avg_climate_layers <- terra::app(
                    climate_layers_i, 
                    fun = mean, na.rm = TRUE, 
                    filename=avg_climate_layers_file,
                    overwrite=TRUE,
                    cores = ncores)
            }
            
            gc()
            
            # 2) calculate the spatial gradient
            
            # Parallel for big areas
            area_i <- SA_i$Areakm2
            
            if(area_i > 10^4){
                
                avg_climate_layers_tiles <- avg_climate_layers
                # define tile resolution 
                nrow(avg_climate_layers_tiles) <- 10
                ncol(avg_climate_layers_tiles) <- 20
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
                        tmp_ext <- terra::ext(avg_climate_layers_tiles, cells = x)
                        real_cells <- terra::cells(avg_climate_layers, tmp_ext)
                        terra::window(avg_climate_layers) <- tmp_ext
                        tmp_data <- spatial_grad(avg_climate_layers)
                        terra::window(avg_climate_layers) <- NULL
                        # plug in "real" icells
                        tmp_data$icell <- real_cells
                        qs::qsave(tmp_data, tmp_file)
                    }
                }, mc.cores = ncores)
                
                spgrad <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
                    tmp_file <- here::here(tmp_dir,
                                           paste(ECO,velocity_variable,"spgrad_tile",x,paste0(polygontogo,".qs"), sep = "_"))
                    try(qs::qread(tmp_file),silent = TRUE)
                })
                
                spgrad <- data.frame(data.table::rbindlist(spgrad))
                
                qs::qsave(spgrad, spgrad_file)
                
                # delete temporary files
                del_files <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
                    tmp_file <- here::here(tmp_dir,
                                           paste(ECO,velocity_variable,"spgrad_tile",x,paste0(polygontogo,".qs"), sep = "_"))
                    unlink(tmp_file)
                })
            } else {
                spgrad = spatial_grad(avg_climate_layers)
                qs::qsave(spgrad, spgrad_file)
            }
            
        }
        
        
        
        
        
        #######
        # Get the spatial gradient up slope (elevation/km)
        cat("spatial gradient up slope\n")
        
        spgrad_ele_file <- 
            here::here(tmp_dir,
                       paste(ECO,velocity_variable,"spgrad_ele",paste0(polygontogo,".qs"), sep = "_"))
        spgrad_ele <- try(qs::qread(spgrad_ele_file,nthreads = ncores),silent = TRUE)
        if(any(class(spgrad_ele)=="try-error")){
            
            # 1) calculate averaged climate layers
            avg_climate_layers_file <- 
                here::here(tmp_dir,
                           paste(ECO,velocity_variable,"tmp_avg_climate_layers",paste0(polygontogo,".tif"), sep = "_"))
            avg_climate_layers <- try(terra::rast(avg_climate_layers_file),silent = TRUE)
            if(any(class(avg_climate_layers)=="try-error")){
                avg_climate_layers <- terra::app(
                    climate_layers_i, 
                    fun = mean, na.rm = TRUE, 
                    filename=avg_climate_layers_file,
                    overwrite=TRUE,
                    cores = ncores)
            }
            
            gc()
            
            # 2) calculate the spatial gradient
            
            # Parallel for big areas
            area_i <- SA_i$Areakm2
            
            if(area_i > 10^4){
                
                avg_climate_layers_tiles <- avg_climate_layers
                # define tile resolution 
                nrow(avg_climate_layers_tiles) <- 10
                ncol(avg_climate_layers_tiles) <- 20
                avg_climate_layers_tiles[] <- 1:ncell(avg_climate_layers_tiles)
                # mask raster
                SA_i_Eck <- terra::project(SA_i,avg_climate_layers_tiles)
                avg_climate_layers_tiles <- terra::mask(avg_climate_layers_tiles, SA_i_Eck)
                
                # force rasters to pair
                elevation <- terra::project(elevation,avg_climate_layers)
                # remove oceans
                elevation[is.na(avg_climate_layers)] <- NA
                
                # plot(avg_climate_layers_tiles);dev.off()
                
                parallel::mclapply(terra::cells(avg_climate_layers_tiles), function(x) {
                    tmp_file <- here::here(tmp_dir,
                                           paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(polygontogo,".qs"), sep = "_"))
                    test <- try(qs::qread(tmp_file),silent = TRUE)
                    if(any(class(test)=="try-error")){
                        tmp_ext <- terra::ext(avg_climate_layers_tiles, cells = x)
                        real_cells <- terra::cells(elevation, tmp_ext)
                        terra::window(elevation) <- tmp_ext
                        tmp_data <- spatial_grad(elevation)
                        terra::window(elevation) <- NULL
                        # plug in "real" icells
                        tmp_data$icell <- real_cells
                        qs::qsave(tmp_data, tmp_file)
                    }
                }, mc.cores = ncores)
                
                spgrad_ele <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
                    tmp_file <- here::here(tmp_dir,
                                           paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(polygontogo,".qs"), sep = "_"))
                    try(qs::qread(tmp_file),silent = TRUE)
                })
                
                spgrad_ele <- data.frame(data.table::rbindlist(spgrad_ele))
                
                qs::qsave(spgrad_ele, spgrad_ele_file)
                
                # delete temporary files
                del_files <- lapply(terra::cells(avg_climate_layers_tiles), function(x) {
                    tmp_file <- here::here(tmp_dir,
                                           paste(ECO,velocity_variable,"spgrad_ele_tile",x,paste0(polygontogo,".qs"), sep = "_"))
                    unlink(tmp_file)
                })
            } else {
                spgrad_ele = spatial_grad(avg_climate_layers)
                qs::qsave(spgrad_ele, spgrad_ele_file)
            }
            
        }
        
        
        # Convert angle to radians
        initial_angle_rad <- deg_to_rad(spgrad$angle) # angle of the spatial gradient 
        target_angle_rad <- deg_to_rad(spgrad_ele$angle) # angle of the elevation up slope
        conversion_rate <- cos(initial_angle_rad - target_angle_rad) 
        
        # Apply conversion >> What is the environmental gradient up slope? (C/km up slope)
        spgrad_ele$Grad_ele <- spgrad_ele$Grad * conversion_rate
        
        #######
        ## calculate gradient-based climate velocity:
        ## Across latitude
        cat("calculate velocity across latitude\n")
        gVelLat <- gVelocity(grad = spgrad, slope = ttrend, 
                             grad_col = "NS", truncate = TRUE)
        
        ## Across elevation
        cat("calculate velocity across elevation\n")
        gVelEle <- gVelocity(grad = spgrad_ele, slope = ttrend, 
                             grad_col = "Grad_ele", truncate = TRUE)
        
        ## Undirectional
        cat("calculate velocity undirectional\n")
        gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
        
        #######
        ## change sign of gVelLat if in the south hemisphere to reflect a velocity away of the tropics
        SouthCells <- terra::as.data.frame(gVelLat, xy = TRUE, cell = TRUE) 
        SouthCells <- SouthCells %>% filter(y<0)
        SouthCells <- SouthCells$cell
        if(length(SouthCells)>0){
            gVelLat[SouthCells] <- gVelLat[SouthCells] * -1
        }
        
        #######
        # summary velocity over study area
        baseline <- mean(global(climate_layers_i,mean,na.rm = TRUE)[,1])
        trend.mean <- as.numeric(global(ttrend,mean,na.rm = TRUE)[,1])
        trend.sd <- as.numeric(global(ttrend,sd,na.rm = TRUE)[,1])
        
        v.mean <- terra::global(gVel$GradVel, mean, na.rm=TRUE)[,1]
        v.median <- terra::global(gVel$GradVel, median, na.rm=TRUE)[,1]
        v.sd <- terra::global(gVel$GradVel, sd, na.rm=TRUE)[,1]
        
        v.lat.mean <- terra::global(gVelLat$GradVel, mean, na.rm=TRUE)[,1]
        v.lat.median <- terra::global(gVelLat$GradVel, median, na.rm=TRUE)[,1]
        v.lat.sd <- terra::global(gVelLat$GradVel, sd, na.rm=TRUE)[,1]
        
        v.ele.mean <- global(gVelEle$GradVel, mean, na.rm=TRUE)[,1]
        v.ele.median <- global(gVelEle$GradVel, median, na.rm=TRUE)[,1]
        v.ele.sd <- global(gVelEle$GradVel, sd, na.rm=TRUE)[,1]
        
        gVelSA_j <- data.frame(baseline, trend.mean, trend.sd, 
                               v.mean, v.median, v.sd, 
                               v.lat.mean, v.lat.median, v.lat.sd, 
                               v.ele.mean, v.ele.median, v.ele.sd)
        names(gVelSA_j) <- paste0(names(gVelSA_j),".",velocity_variable)
        
        gVelSA[[v]] <- gVelSA_j
    } 
    
    gVelSA <- do.call(cbind,gVelSA)
    gVelSA$ID=polygontogo
    
    
} else {
    # If marine
    # Calculate only for SST
    
    velocity_variable <- "sst"
    
    gVelSA <- data.frame()
    
    # select the climate variable
    climate_layers_i <- climate_layers[[which(names(climate_layers)=="mean")]]
    
    # crop layers to the study area
    climate_layers_i <- terra::mask(terra::crop(climate_layers_i, SA_i), SA_i)
    
    # project to equal-area
    climate_layers_i <- terra::project(climate_layers_i,Eckt)
    
    # if study area is two small (there are less then 8 cells), there is no way climate gradients can be calculated.
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
    cat("spatial gradient\n")
    spgrad = spatial_grad(climate_layers_i)
    
    #######
    ## calculate gradient-based climate velocity:
    cat("calculate gradient-based climate velocity\n")
    gVel <- gVelocity(grad = spgrad, slope = ttrend,
                      truncate = TRUE)
    
    gVelLat <- gVelocity(grad = spgrad, slope = ttrend, 
                         grad_col = "NS", truncate = TRUE)
    
    #######
    ## change sign of gVelLat if in the south hemisphere to reflect a velocity away of the tropics
    SouthCells <- terra::as.data.frame(gVelLat, xy = TRUE, cell = TRUE) 
    SouthCells <- SouthCells %>% filter(y<0)
    SouthCells <- SouthCells$cell
    if(length(SouthCells)>0){
        gVelLat[SouthCells] <- gVelLat[SouthCells] * -1
    }
    
    #######
    # summary velocity over study area
    baseline <- mean(global(climate_layers_i,mean,na.rm = TRUE)[,1])
    trend.mean <- as.numeric(global(ttrend,mean,na.rm = TRUE)[,1])
    trend.sd <- as.numeric(global(ttrend,sd,na.rm = TRUE)[,1])
    
    v.mean <- terra::global(gVel$GradVel, mean, na.rm=TRUE)[,1]
    v.median <- terra::global(gVel$GradVel, median, na.rm=TRUE)[,1]
    v.sd <- terra::global(gVel$GradVel, sd, na.rm=TRUE)[,1]
    
    v.lat.mean <- as.numeric(global(gVelLat$GradVel, mean, na.rm=TRUE))
    v.lat.median <- as.numeric(global(gVelLat$GradVel, median, na.rm=TRUE))
    v.lat.sd <- as.numeric(global(gVelLat$GradVel, sd, na.rm=TRUE))
    
    gVelSA <- data.frame(baseline, trend.mean, trend.sd, 
                         v.mean, v.median, v.sd, 
                         v.lat.mean, v.lat.median, v.lat.sd)
    names(gVelSA) <- paste0(names(gVelSA),".",velocity_variable)
    gVelSA$ID=polygontogo
}

# save velocity results
write.csv(gVelSA, here::here(velocity_SA_dir, paste0(polygontogo,".csv")), row.names = FALSE)

# save velocity maps
terra::writeRaster(gVel$GradVel, here::here(velocity_SA_dir, paste0(polygontogo,"_gVel.tif")),overwrite=TRUE)
terra::writeRaster(gVelLat$GradVel, here::here(velocity_SA_dir, paste0(polygontogo,"_gVelLat.tif")),overwrite=TRUE)
if(ECO=="Ter"){
    terra::writeRaster(gVelEle$GradVel, here::here(velocity_SA_dir, paste0(polygontogo,"_gVelEle.tif")),overwrite=TRUE) 
    
    # delete temporary files
    unlink(list.files(tmp_dir,pattern = polygontogo,full.names = TRUE))
    
}
