# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","elevatr","raster","sf","rgdal","Hmisc","dplyr", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

########################
# load functions
source(here::here("R/range_shift_functions.R"))

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
sptogo <- as.character(paste(command_args[1], collapse = " "))
realm <- as.character(paste(command_args[2], collapse = " "))

# sptogo <- SDMsSpList[1]
# sptogo <- "Abra_alba"
# sptogo <- "Ampelisca_abdita"
# realm = "Mar"

output.dir <- here::here("/media/seagate/boliveira/SDMs/MaxNet",realm)

########################
# Elevation map
if(realm == "Ter"){
    if(!file.exists(here::here("/media/seagate/boliveira/Land","elev_raster_ter.tif"))){
        ter.ras <- terra::rast(here::here("/media/seagate/boliveira/Land","model_raster_ter.tif"))
        Elev <- get_elev_raster(raster(ter.ras), z = 6)
        Elev <- project(rast(Elev), ter.ras)
        # save
        writeRaster(Elev, here::here("/media/seagate/boliveira/Land","elev_raster_ter.tif"))
    } else {
        Elev <- rast(here::here("/media/seagate/boliveira/Land","elev_raster_ter.tif"))
    }
}
# exposure variables
if(realm == "Ter"){
    exp.vars <- c("tas","pr")
} 
if(realm == "Mar"){
    exp.vars <- c("sst")
}

########################
# Load files

# bioshift data
bioshift_info <- read.csv(here::here(output.dir,sptogo,paste(sptogo,"bioshift.csv")))

# shift info
shift_info <- read.csv(here::here(output.dir,sptogo,paste(sptogo,"shift_info.csv")))
Time_periods = shift_info$time_period

# Study area layer
fgdb <- here::here("Data/Study_Areas_v1/Study_Areas.gdb")
StudyID <- unique(bioshift_info$ID)
StudyArea <- lapply(StudyID, function(x) {
    tmp <- readOGR(dsn=fgdb,layer=x)
    vect(tmp)
})
StudyArea <- vect(StudyArea)
StudyArea$NAME = StudyID

# range shift SA
sdms <- list.files(here::here(output.dir,sptogo), pattern = "SA", full.names = TRUE)
sdms <- sdms[grep("tif",sdms)]
sdms <- sdms[-grep("bios",sdms)]

sdms_SA <- lapply(sdms, rast)

# range shift BG
sdms <- list.files(here::here(output.dir,sptogo), pattern = "BG", full.names = TRUE)
sdms <- sdms[grep("tif",sdms)]
sdms <- sdms[-grep("bios",sdms)]

sdms_BG <- rast(sdms)

# shift type
type <- unique(bioshift_info$Type)

#####################################
# Calc SA range shifts 

if(!file.exists(here::here(output.dir,sptogo,paste(sptogo, "shifts_SA.csv")))){
    sui_SA_cells <- list()
    shifts_SA <- list()
    for(i in 1:length(sdms_SA)){
        
        sui <- terra::spatSample(sdms_SA[[i]], 10000, na.rm = TRUE, cells = TRUE, xy = TRUE)
        sui_SA_cells[[i]] <- sui$cell
        
        shifts_SA_tmp <- list()
        for(j in 1:length(type)){
            
            if(type[i] == "LAT"){
                # get type
                var_t <- sui$y
            }
            if(type[i] == "ELE"){
                # crop elevation data to the SA
            }
            
            tmp <- range_shift(sui_t1 = sui[,4], sui_t2 = sui[,5], var_t = var_t)
            info <- names(sdms_SA[[i]])
            info1 <- strsplit(info[1]," ")[[1]]
            info2 <- strsplit(info[2]," ")[[1]]
            mysps <- info1[1]
            myid <- info1[2]
            mystart <- info1[3]
            myend <- info2[3]
            tmp$ID <- myid
            tmp$Type <- type[j]
            tmp$START <- mystart
            tmp$END <- myend
            tmp$Species <- mysps
            
            shifts_SA_tmp[[j]] <- tmp
        }
        shifts_SA_tmp <- rbindlist(shifts_SA_tmp)
        shifts_SA[[i]] <- shifts_SA_tmp
    }
    shifts_SA <- rbindlist(shifts_SA)
    
    write.csv(shifts_SA, 
              here::here(output.dir,sptogo,paste(sptogo, "shifts_SA.csv")),
              row.names = FALSE)
}


#####################################
# Calc BG range shifts 
# get suitability
sui <- terra::spatSample(sdms_BG, 10000, na.rm = TRUE, cells = TRUE, xy = TRUE)

if(!file.exists(here::here(output.dir,sptogo,paste(sptogo, "shifts_BG.csv")))){
    shifts_BG <- list()
    
    for(i in 1:length(type)){
        
        if(type[i] == "LAT"){
            # get type
            var_t <- sui$y
        }
        if(type[i] == "ELE"){
            # crop elevation data to the SA
        }
        
        tmp <- range_shift(sui_t1 = sui[,4], sui_t2 = sui[,5], var_t = var_t)
        info <- names(sdms_BG)
        info1 <- strsplit(info[1]," ")[[1]]
        info2 <- strsplit(info[2]," ")[[1]]
        mysps <- info1[1]
        myid <- info1[2]
        mystart <- info1[3]
        myend <- info2[3]
        tmp$ID <- myid
        tmp$Type <- type[j]
        tmp$START <- mystart
        tmp$END <- myend
        tmp$Species <- mysps
        
        shifts_BG[[j]] <- tmp
    }
    shifts_BG <- rbindlist(shifts_BG)
    
    write.csv(shifts_BG, 
              here::here(output.dir,sptogo,paste(sptogo, "shifts_BG.csv")),
              row.names = FALSE)
}

#####################################
# Range position
# Latitudinal range filling of the study area relative to the BG area

if(!file.exists(here::here(output.dir,sptogo,paste(sptogo, "range_position.csv")))){
    
    range_position <- list()
    for(i in 1:nrow(shift_info)){
        range_position_tmp <- list()
        for(j in 1:length(type)){
            range_position_tmp <- range_pos(shifts_SA[which(shifts_SA$ID==shift_info$ID[i]),],
                                            shifts_BG[which(shifts_SA$Type==type[j]),])
            range_position_tmp$type = type[j]
            range_position_tmp$SA = shift_info$ID[i]
        }
        range_position[[i]] <- range_position_tmp
    }
    range_position <- rbindlist(range_position)
    
    write.csv(range_position, 
              here::here(output.dir,sptogo,paste(sptogo, "range_position.csv")),
              row.names = FALSE)
}
#####################################
# Calc exposure at the BG
# mean weighted env differences from t1 >> t2

if(!file.exists(here::here(output.dir,sptogo,paste(sptogo, "exposure_BG.csv")))){
    
    BG_start <- lapply(1:length(Time_periods), function(i) {
        rast(here::here(output.dir,sptogo,paste(sptogo,Time_periods[i],"BG start bios.tif")))
    })
    BG_start <- BG_start[[1]]
    BG_start <- terra::as.data.frame(BG_start[sui$cell])
    BG_start <- BG_start * sui[,4]
    
    BG_end <- lapply(1:length(Time_periods), function(i) {
        rast(here::here(output.dir,sptogo,paste(sptogo,Time_periods[i],"BG end bios.tif")))
    })
    BG_end <- BG_end[[1]]
    BG_end <- terra::as.data.frame(BG_end[sui$cell])
    BG_end <- BG_end * sui[,5]
    
    # exposure at the background
    BG_exp <- BG_end - BG_start
    BG_exp <- data.frame(apply(BG_exp, 2, mean))
    colnames(BG_exp) <- "Exposure BG"
    
    write.csv(BG_exp, 
              here::here(output.dir,sptogo,paste(sptogo, "exposure_BG.csv")),
              row.names = FALSE)
    
}
#####################################
# Calc exposure at the SA
# mean weighted env differences from t1 >> t2

if(!file.exists(here::here(output.dir,sptogo,paste(sptogo, "exp_SA.csv")))){
    
    exp_SA <- list()
    for(i in 1:nrow(shift_info)){
        
        sui <- sdms_SA[[i]][sui_SA_cells[[i]]]
        
        SA_start <- lapply(1:length(Time_periods), function(i) {
            tmp <- rast(here::here(output.dir,sptogo,paste(sptogo,shift_info$ID[i],Time_periods[i],"SA start bios.tif")))
            tmp <- terra::as.data.frame(tmp[sui_SA_cells[[i]]])
            tmp * sui[,1]
        })
        
        SA_end <- lapply(1:length(Time_periods), function(i) {
            tmp <- rast(here::here(output.dir,sptogo,paste(sptogo,shift_info$ID[i],Time_periods[i],"SA end bios.tif")))
            tmp <- terra::as.data.frame(tmp[sui_SA_cells[[i]]])
            tmp * sui[,1]
        })
        
        # exposure at the background
        tmp <- lapply(1:length(Time_periods), function(i) { 
            tmp <- SA_end[[i]] - SA_start[[i]]
            tmp <- data.frame(apply(tmp, 2, mean))
            colnames(tmp) <- "Exposure SA"
            tmp$bio <- rownames(tmp)
            tmp$Time_periods <- Time_periods[i]
            return(tmp)
        })
        
        exp_SA[[i]] <- rbindlist(tmp)
    }
    exp_SA <- rbindlist(exp_SA)
    
    write.csv(exp_SA, 
              here::here(output.dir,sptogo,paste(sptogo, "exp_SA.csv")),
              row.names = FALSE)
    
}
