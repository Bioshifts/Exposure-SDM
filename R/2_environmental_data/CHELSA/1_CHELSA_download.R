### Download CHELSA climate variables

# Data from Chelsa (https://chelsa-climate.org/downloads/)
# 
# Variables of interest following Vanderwal 2013 NatClimChange (mean, minimum, maximum and standard deviation of monthly temperature, sum and coefficient of variation of precipitation, and sum of wettest and driest quarter for precipitation):  
# * Mean monthly daily air temperature >> tas  *
# * Mean monthly minimum air temperature >> tasmin  * 
# * Mean monthly maximum air temperature >> tasmax  * 
# * Total monthly precipitation >> pr  * 
# 
# From these variables we will create these others variables by integrating over the last 12 months from species occurrence:
# * Mean annual air temperature >> tas  *
# * Mean annual minimum air temperature >> tasmin  * 
# * Mean annual maximum air temperature >> tasmax  * 
# * Temperature seasonality >> SD(tas)  * 
# * Annual precipitation amount >> pr  * 
# * Precipitation seasonality >> CV(pr)  * 
# * Mean monthly precipitation amount of the wettest quarter >> The wettest quarter of the year is determined (to the nearest month)  * 
# * Mean monthly precipitation amount of the driest quarter >> The driest quarter of the year is determined (to the nearest month)  * 
# 
# DO NOT PARALLELIZE THIS WORK! YOU ARE LIMITED BY NETWORK SPEED. IF PARALLELIZING YOU END UP WITH A LOT OF SLOW DOWNLOADS.
# MAYBE IT IS POSSIBLE TO PARALLELIZE IN COMPUTER CLUSTER IF NETWORK SPEED IS SHARED ACROSS CORES (?)

##########################

### SETUP

list.of.packages <- c("terra","here",
                      "reshape2","tidyr","tidyverse",
                      "foreach","doParallel","pbapply")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

# increase default timeout (=60s) prevents downloads to stop before finishing
options(timeout=3000) 

# set computer
computer = "muse"

# working directory
work_dir <- 
    if(computer == "muse"){
        "/storage/simple/projects/t_cesab/brunno/Exposure-SDM"
    }
setwd(work_dir)

# source settings
source("R/settings.R")

# load functions and settings
source("R/getChelsa.R")
source("R/my_functions.R")

# settings
t.period <- temporal_range_env_data(realm = "Ter") # Limits of the CHELSA database

realm = "Ter"

# path to save model raster
nc_path <- vars_dir(realm = realm) 
if(!dir.exists(here::here(nc_path))){
    dir.create(here::here(nc_path), recursive = TRUE)
}

# raw CHELSE data are save here
scratch_dir_1km <- here::here(scratch_dir,"Data/cruts_1km")
scratch_dir_5km <- here::here(scratch_dir,"Data/cruts_5km")

periods <- format(
    seq.Date(from = as.Date(paste0(t.period[1],"/01/01")), 
             to = as.Date(paste0(t.period[2],"/12/01")), 
             by = "month"), 
    "%m_%Y")


vars <- c("tmax", "tmin", "prec")

# create folders to store variables
for (i in 1:length(vars)){
    if(!dir.exists(here::here(scratch_dir_1km,vars[i]))){
        dir.create(here::here(scratch_dir_1km,vars[i]),recursive = TRUE)
    }
    if(!dir.exists(here::here(scratch_dir_5km,vars[i]))){
        dir.create(here::here(scratch_dir_5km,vars[i]),recursive = TRUE)
    }
}

#############################
# Create land mask
# -32768 is the value that represents ocean in cruts rasters

if(!file.exists(here::here(nc_path,"model_raster_ter_1km.tif"))){
    
    # download any temp data from CHELSA
    getChelsa_cruts(myvar = vars[1], 
                    period = periods[1], 
                    dest.folder = scratch_dir, 
                    over.write = TRUE)
    
    my_temp <- here::here(scratch_dir,vars[1],paste0(vars[1],"_",periods[1],".tif"))
    my_temp <- terra::rast(my_temp)
    # add some range for the mask value representing the ocean value just for making sure we nothing weird happens
    my_mask_temp <- terra::ifel(my_temp < -10000, 0, 1)
    # plot(my_mask_temp)
    # Crop Antarctic cells
    ant_ext <- ext(my_mask_temp)
    ant_ext[4] <- -60
    my_mask_temp[ant_ext] <- 0
    
    # download any precip data from CHELSA
    getChelsa_cruts(myvar = vars[3], 
                    period = periods[3], 
                    dest.folder = scratch_dir, 
                    over.write = TRUE)
    
    my_precip <- here::here(scratch_dir,vars[3],paste0(vars[3],"_",periods[1],".tif"))
    my_precip <- terra::rast(my_precip)
    # for precip, there should not be any prec < 0
    my_mask_precip <- terra::ifel(my_precip < 0, 0, 1)
    # plot(my_mask_precip)
    # Crop Antarctic cells
    ant_ext <- ext(my_mask_precip)
    ant_ext[4] <- -60
    my_mask_precip[ant_ext] <- 0
    
    # Combine both masks
    my_mask <- my_mask_temp * my_mask_precip
    
    # save
    terra::writeRaster(my_mask, 
                       here::here(nc_path,"model_raster_ter_1km.tif"), 
                       overwrite = TRUE)
    
    # get land mask at 5 km
    my_mask <- terra::rast(here::here(nc_path,"model_raster_ter_1km.tif"))
    my_mask <- terra::aggregate(my_mask, fact = 5)
    my_mask <- terra::ifel(my_mask == 1, 1, 0)
    
    # save
    terra::writeRaster(my_mask, 
                       here::here(nc_path,"model_raster_ter_5km.tif"), 
                       overwrite = TRUE)
    
}

my_mask_1km <- rast(here::here(nc_path,"model_raster_ter_1km.tif"))
my_mask_5km <- rast(here::here(nc_path,"model_raster_ter_5km.tif"))

#############################
### DOWNLOAD raw cruts data

# DO NOT MASK RASTER FILES WITH 1KM BECAUSE:
# 1) IT IS TIME CONSUMING
# 2) NEW SAVED RASTER ARE MUCH LARGER THEN THE ONES DONWLOADED DIRECTLY FROM CHELSA
# SOLUTION:
# WHEN EXTRACTING ENVIRONMENTAL DATA FROM OCCURRENCE
# 1) REMOVE CELLS FALLING IN THE OCEAN BASED ON THE MASK RASTER
# WHEN CALCULATING BIOCLIMATICS
# 1) STACK LAYERS FOR THE YEAR i
# 2) MASK
# 3) CALCULATE BIOCLIMATICS

for (v in 1:length(vars)){
    
    var_going <- vars[v]
    
    cat("\rDownloading", var_going)
    
    ncores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
    cat("N of cores is", ncores)
    
    # loop trough periods
    parallel::mclapply(mc.cores = ncores,
                       1:length(periods), 
                       function(i){
                           
                           try({
                               myfile_1km <- here::here(scratch_dir_1km, var_going, paste0(var_going,"_",periods[i],".tif"))
                               myfile_5km <- here::here(scratch_dir_5km, var_going, paste0(var_going,"_",periods[i],".tif"))
                               
                               if(!file.exists(myfile_1km)){
                                   getChelsa_cruts(myvar = var_going, 
                                                   period = periods[i], 
                                                   dest.folder = scratch_dir_1km, 
                                                   over.write = TRUE)
                                   
                               }
                               
                               
                               if(!file.exists(myfile_5km)){
                                   # load in
                                   tmp <- terra::rast(myfile_1km)
                                   
                                   # convert to 5km
                                   tmp <- terra::aggregate(tmp, fact = 5)
                                   
                                   # save
                                   terra::writeRaster(
                                       tmp, 
                                       myfile_5km, 
                                       overwrite = TRUE)
                               }
                               
                           }, 
                           silent = TRUE)
                           
                       })
    
}


########################################################################
# test if everything has been downloaded
########################################################################

##########
# 1km ----

# First test: Check if all periods were downloaded
my.test <- list()
for(i in 1:length(vars)){
    # Files should look  like this 
    myfiles <- here::here(scratch_dir_1km, vars[i], paste0(vars[i],"_",periods,".tif"))
    # downloaded files
    mydowns <- list.files(here::here(scratch_dir_1km,vars[i]), full.names = TRUE)
    # test
    my.test[[i]] <- myfiles[which(!myfiles %in% mydowns)]
}
names(my.test) <- vars
# my.test[[3]]

for(v in 1:length(my.test)){
    
    tmp <- my.test[[v]]
    
    # check if there is any missing
    if(length(tmp)>0){ 
        
        # if yes, try downloading again...
        for(i in 1:length(tmp)){
            
            tmp_file <- tmp[i]
            
            var_going <- vars[v]
            
            period_i <- gsub(scratch_dir,"",tmp_file)
            period_i <- gsub(paste0("/",var_going,"/",var_going,"_"),"",period_i)
            period_i <- gsub(".tif","",period_i)
            
            cat("\rDownloading", var_going)
            
            getChelsa_cruts(myvar = var_going, 
                            period = period_i, 
                            dest.folder = scratch_dir, 
                            over.write = TRUE)
        }
    }
}


# Second test: Check if all files open correctly
for(v in 1:length(vars)){
    
    var_going <- vars[v]
    
    cat("\rChecking", var_going,"\n")
    
    for(i in 1:length(periods)){
        
        cat("\rperiod",i,"from",length(periods),"..... period", periods[i])
        
        filetosave <- here::here(scratch_dir_1km, var_going, paste0(var_going,"_",periods[i],".tif"))
        
        tmp <- try(rast(filetosave), silent = TRUE)
        
        if(class(tmp)=='try-error'){ # download again
            
            cat("Downloading again for period", periods[i],"\n")
            
            getChelsa_cruts(myvar = var_going, 
                            period = periods[i], 
                            dest.folder = scratch_dir_1km, 
                            over.write = TRUE)
        }
    }
}





##########
# 5km ----

# First test: Check if all periods were downloaded
my.test <- list()
for(i in 1:length(vars)){
    # Files should look  like this 
    myfiles <- here::here(scratch_dir_5km, vars[i], paste0(vars[i],"_",periods,".tif"))
    # downloaded files
    mydowns <- list.files(here::here(scratch_dir_5km,vars[i]), full.names = TRUE)
    # test
    my.test[[i]] <- myfiles[which(!myfiles %in% mydowns)]
}
names(my.test) <- vars
# my.test[[3]]

for(v in 1:length(my.test)){
    
    tmp <- my.test[[v]]
    
    # check if there is any missing
    if(length(tmp)>0){ 
        
        # if yes, try downloading again...
        for(i in 1:length(tmp)){
            
            tmp_file <- tmp[i]
            
            var_going <- vars[v]
            
            period_i <- gsub(scratch_dir,"",tmp_file)
            period_i <- gsub(paste0("/",var_going,"/",var_going,"_"),"",period_i)
            period_i <- gsub(".tif","",period_i)
            
            cat("\rDownloading", var_going)
            
            getChelsa_cruts(myvar = var_going, 
                            period = period_i, 
                            dest.folder = scratch_dir, 
                            over.write = TRUE)
        }
    }
}

########################################################################
# Second test: Check if all files open correctly
for(v in 1:length(vars)){
    
    var_going <- vars[v]
    
    cat("\rChecking", var_going,"\n")
    
    for(i in 1:length(periods)){
        
        cat("\rperiod",i,"from",length(periods),"..... period", periods[i])
        
        filetosave <- here::here(scratch_dir_5km, var_going, paste0(var_going,"_",periods[i],".tif"))
        
        tmp <- try(rast(filetosave), silent = TRUE)
        
        if(class(tmp)=='try-error'){ # download again
            
            cat("Downloading again for period", periods[i],"\n")
            
            getChelsa_cruts(myvar = var_going, 
                            period = periods[i], 
                            dest.folder = scratch_dir_5km, 
                            over.write = TRUE)
        }
    }
}

########################################################################
# Visual inspection
# plot all

png("test")
