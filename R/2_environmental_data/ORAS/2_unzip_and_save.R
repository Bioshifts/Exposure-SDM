#### Unzip ORAS data

library(ncdf4)
library(raster)
library(terra)
library(pbapply)
library(parallel)
library(sf)
library(ggplot2)


# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

# load functions
source("R/my_functions.R")
source("R/settings.R")

realm = "Mar"
nc_path <- vars_dir(realm = realm) 

# raw files will be saved to scratch dir
zipdir <- here::here(scratch_dir,"Data/SST_ORAS_downloads")
    
# create folder to store rasters
scratch_dir <- here::here(scratch_dir,"Data/SST")

if(!dir.exists(scratch_dir)){
    dir.create(scratch_dir)
}




# list of files
my_files <- list.files(zipdir,full.names = T,pattern = ".zip")

# temporary directory
# create temporary folder 
tmp_dir <- here::here(scratch_dir,"tmp")

# unzip and save data
for(i in 1:length(my_files)){ cat ("file", i, "from", length(my_files))
    
    # create temporary folder 
    if(!dir.exists(tmp_dir)){
        dir.create(tmp_dir)
    }
    
    # unzip year i
    unzip(my_files[i], exdir = tmp_dir)
    
    # list files year i
    my_files_i = list.files(tmp_dir,full.names = T,pattern = ".nc")
    
    file_names <- gsub(tmp_dir,"",my_files_i)
    file_names <- gsub("/sosstsst_control_monthly_highres_2D_","",file_names)
    file_names <- gsub("_CONS_v0.1.nc","",file_names)
    file_names <- sapply(file_names, function(x) {
        paste(c(substr(x,1,4),substr(x,5,6)),collapse = "_")
    })
    
    # load coord data
    nc_data <- nc_open(my_files_i[1])
    
    lon <- ncvar_get(nc_data, "nav_lon")
    lat <- ncvar_get(nc_data, "nav_lat", verbose = F)
    fillvalue <- ncatt_get(nc_data, "sosstsst", "_FillValue")
    
    nc_close(nc_data)
    
    # empty raster
    geo.r <- terra::rast(terra::ext(-180,180,-64,85))
    terra::res(geo.r) <- c(0.25,0.25)
    geo.r[] <- NA
    
    # get all sst data for year i
    tmp = mclapply(1:length(my_files_i), function(x) {
        
        nc_data <- ncdf4::nc_open(my_files_i[x])
        tmp <- ncdf4::ncvar_get(nc_data, "sosstsst")
        tmp[tmp==fillvalue$value] <- NA
        
        ncdf4::nc_close(nc_data)
        
        # vector from sst i
        remap.tbl <- data.frame(lon=as.vector(lon), lat=as.vector(lat), sst = as.vector(tmp))
        remap.tbl <- na.omit(remap.tbl)
        remap.tbl <- terra::vect(remap.tbl, geom=c("lon", "lat"))
        
        # vector from sst i based on geo.r
        new_raster <- terra::rasterize(remap.tbl, geo.r, field = "sst")
        # par(mfrow = c(1,2))
        plot(new_raster)
        
        # workaround for weird projection from ORAS
        w <- matrix(1, 3, 3)  
        new_raster <- terra::focal(new_raster, w, mean, na.rm=TRUE, NAonly=TRUE, pad=TRUE)
        plot(new_raster)
        
        terra::writeRaster(new_raster, 
                           here::here(scratch_dir,paste0("SST_",file_names[x],".tif")), 
                           overwrite = TRUE)
        
    }, mc.cores = 12)
    
    
    # clean temporary files
    unlink(tmp_dir, recursive = TRUE)
    
}

# renames SST data to match CHELSA: varname_month_year
my_vars <- list.files(here::here(scratch_dir))
rename_vars <- sapply(my_vars, function(x){
    tmp = strsplit(x,"_")[[1]]
    tmp = paste0(paste(tmp[1],gsub('.tif','',tmp[3]),tmp[2],sep = "_"),'.tif')
    return(tmp)
})

file.rename(here::here(scratch_dir,my_vars),here::here(scratch_dir,rename_vars))

# my_vars <- list.files(here::here(scratch_dir),full.names = T)
# todelete <- grep("1979",my_vars)
# todelete <- my_vars[todelete]
# unlink(todelete)

#######################
# Save model raster
# To be used for extracting cell ID

tmp <- list.files(here::here(nc_path))[1]
tmp <- rast(here::here(nc_path,tmp))

tmp <- is.na(tmp)
tmp <- ifel(tmp,0,1)
names(tmp) <- "model_raster_mar"
plot(tmp)

writeRaster(tmp, here::here(nc_path,"model_raster_mar.tif"), overwrite = TRUE)

