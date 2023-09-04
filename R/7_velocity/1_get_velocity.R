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
# polygontogo <- "A112_P1" # Mar # South
# polygontogo <- "A43_P1" # Mar # North

print(polygontogo)

cat("\rrunning polygon", polygontogo)

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

# source settings
source("R/settings.R")
source("R/velocity_functions.R")

# create dir to store results
if(!dir.exists(velocity_SA_dir)){
    dir.create(velocity_SA_dir)
}

# velocity variables
velocity_variables <- c("mat","map")

# gradients
gradients <- c("LAT","ELE")


########################
# load SA polygon
SA_i <- terra::vect(here::here(SA_shps_dir,paste0(polygontogo,".shp")))

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
    Bioshifts_DB$Paper.ID
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
S_time <- S_start:S_end


# get layers within time period of shift
if(ECO=="Ter"){
    vars_dir <- here::here(vars_dir(ECO),paste0("bio_proj_",my_res))
} else {
    vars_dir <- here::here(vars_dir(ECO),"SST")
}
climate_layers <- list.files(vars_dir)
climate_layers_pos <- grepl(paste(S_time,collapse = "|"),climate_layers)
climate_layers <- climate_layers[climate_layers_pos]
climate_layers <- terra::rast(here::here(vars_dir,climate_layers))

########################
## Calculate velocity at the SA

# If terrestrial
# Calculate velocity for LAT and ELE and for each climatic variable
if(ECO=="Ter"){
    
    # crop elevation to the study area
    elevation_i <- terra::mask(terra::crop(elevation, SA_i), SA_i)
    plot(elevation_i)
    rmote::stop_rmote()
    # project to equal-area
    elevation_i <- terra::project(elevation_i,Eckt)
    
    gVelSA <- data.frame()
    
    for(v in 1:length(velocity_variables)){ # for each climate variable
        
        # select the climate variable
        climate_layers_i <- climate_layers[[which(names(climate_layers)==velocity_variables[v])]]
        
        # crop layers to the study area
        climate_layers_i <- terra::mask(terra::crop(climate_layers_i, SA_i), SA_i)
        
        # project to equal-area
        climate_layers_i <- terra::project(climate_layers_i,Eckt)
        
        # force raster to pair
        elevation_i <- terra::project(elevation_i,climate_layers_i)
        
        # in CHELSA, temperature is *10
        if(velocity_variables[v]=="mat"){
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
            
            # Force smaller resolution
            elevation_i <- terra::disagg(
                elevation_i, 
                fact = 12/ncell(elevation_i),
                method = "bilinear")
            
        } 
        
        #######
        # calculate the trend (C/year)
        ttrend = temp_grad(climate_layers_i,
                           th = 0.25*nlyr(climate_layers_i) ## set minimum # obs. to 1/4 time series length
        )
        
        #######
        # calculate the spatial gradient (C/km)
        spgrad = spatial_grad(climate_layers_i)
        # from meters to km
        spgrad$WE <- spgrad$WE / 1000
        spgrad$NS <- spgrad$NS / 1000
        spgrad$Grad <- spgrad$Grad / 1000
        
        #######
        # Get the spatial gradient up slope
        spgrad_ele = spatial_grad(elevation_i)
        
        # Convert angles to radians
        initial_angle_rad <- .rad(spgrad$angle) # angle of the spatial gradient 
        target_angle_rad <- .rad(spgrad_ele$angle) # angle of the elevation up slope
        
        # Apply conversion >> What is the environmental gradient up slope? 
        spgrad$angle_ele <- spgrad_ele$angle
        spgrad$Grad_ele <- spgrad$Grad * cos(initial_angle_rad - target_angle_rad) 
        
        #######
        ## calculate gradient-based climate velocity:
        # gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
        
        gVelLat <- gVelocity(grad = spgrad, slope = ttrend, 
                             grad_col = "NS", truncate = TRUE)
        
        gVelEle <- gVelocity(grad = spgrad, slope = ttrend, 
                             grad_col = "Grad_ele", truncate = TRUE)
        
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
        
        v.lat.mean <- global(gVelLat, mean, na.rm=TRUE)
        v.lat.median <- global(gVelLat, median, na.rm=TRUE)
        v.lat.sd <- global(gVelLat, sd, na.rm=TRUE)
        
        v.ele.mean <- global(gVelLat, mean, na.rm=TRUE)
        v.ele.median <- global(gVelLat, median, na.rm=TRUE)
        v.ele.sd <- global(gVelLat, sd, na.rm=TRUE)
        
        gVelSA_j <- data.frame(baseline, trend.mean, trend.sd, 
                               v.lat.mean, v.lat.median, v.lat.sd, 
                               v.ele.mean, v.ele.median, v.ele.sd)
        names(gVelSA_j) <- paste0(names(gVelSA_j),".",velocity_variables[v])
        gVelSA_j$ID=polygontogo
    } 
    gVelSA <- rbind(gVelSA,gVelSA_j)
    
    
} else {
    # If marine
    # Calculate only for SST
    
    gVelSA <- data.frame()
    
    # crop layers to the study area
    climate_layers_i <- terra::mask(terra::crop(climate_layers, SA_i), SA_i)
    
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
    ttrend = temp_grad(climate_layers_i,
                       th = 0.25*nlyr(climate_layers_i) ## set minimum # obs. to 1/4 time series length
    )
    
    #######
    # calculate the spatial gradient (C/km)
    spgrad = spatial_grad(climate_layers_i)
    # from meters to km
    spgrad$WE <- spgrad$WE / 1000
    spgrad$NS <- spgrad$NS / 1000
    spgrad$Grad <- spgrad$Grad / 1000
    # na.omit(spgrad)
    
    #######
    ## calculate gradient-based climate velocity:
    # gVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
    gVelLat <- gVelocity(grad = spgrad, slope = ttrend, grad_col = "NS", truncate = TRUE)
    
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
    
    v.lat.mean <- as.numeric(global(gVelLat, mean, na.rm=TRUE))
    v.lat.median <- as.numeric(global(gVelLat, median, na.rm=TRUE))
    v.lat.sd <- as.numeric(global(gVelLat, sd, na.rm=TRUE))
    
    gVelSA <- data.frame(baseline = baseline, 
                         trend.mean = trend.mean, 
                         trend.sd = trend.sd, 
                         v.lat.mean = v.lat.mean, 
                         v.lat.median = v.lat.median, 
                         v.lat.sd = v.lat.sd)
    names(gVelSA) <- paste0(names(gVelSA),".sst")
    gVelSA$ID=polygontogo
    gVelSA2=gVelSA
    gVelSA1
    gVelSA2
}

write.csv(gVelSA, here::here(velocity_SA_dir, paste0(polygontogo,".csv")), row.names = FALSE)


gVelocity(spatial_grad(climate_layers_i),
          temp_grad(climate_layers_i,
                    th = 0.25*nlyr(climate_layers_i)))

VoCC::gVoCC(VoCC::tempTrend(raster::stack(climate_layers_i),
                            th = 0.25*nlyr(climate_layers_i)), 
            VoCC::spatGrad(raster::stack(climate_layers_i)))
