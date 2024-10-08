###########
# Calculate bioclimatics from monthly raster data using dismo

rm(list=ls())
gc()

# Run to year i
# load period to go
command_args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(command_args[1])

library(dismo)
library(raster)
library(terra)
library(parallel)
library(pbapply)
library(qs)
library(tictoc)

# detect slurm cores
N_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

# Check if file exists
check_if_file_exists <- TRUE

# set realm
realm = "Ter"

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

# load functions and settings
source("R/getChelsa.R")
source("R/my_functions.R")

# get vars
my_vars <- myvars(realm)

# get directory to save bioclimatics
vars_dir <- vars_dir(realm = realm)
bios_dir <- bios_dir(realm)

# create dir to store bioclimatics for each period
if(!dir.exists(bios_dir)){
    dir.create(bios_dir)
}

# get mask
my_mask <- rast(here::here(vars_dir,paste0("model_raster_ter_",my_res,".tif")))

# all layers
all_layers_names <- list.files(here(vars_dir,my_res), pattern = ".tif", recursive = TRUE)
all_layers_names <- gsub('.tif','',all_layers_names)
all_layers_names <- sapply(all_layers_names, function(x) strsplit(x,"/")[[1]][2])

all_layers <- list.files(here(vars_dir,my_res), pattern = ".tif",full.names = TRUE, recursive = TRUE)

# possible years
my_yrs <- sapply(all_layers_names, function(x) strsplit(x,"_")[[1]][3])
my_yrs <- sort(as.numeric(unique(my_yrs)))

# remove the last - n_yr_bioclimatic (cannot calculate bioclimate for the last year)
my_yrs <- my_yrs[-1:-n_yr_bioclimatic]

# This function cant be serialized as SpatRaster cannot be sent to computer nodes. A solution is submitting jobs for each year

cat("\rperiod",i,"from",length(my_yrs))

year_i <- my_yrs[i]

# file name to save
filetosave = here::here(bios_dir,paste0("bios_",realm,"_",year_i,".tif"))

# check if file exists
if(check_if_file_exists){
    RUN <- !file.exists(filetosave)
} else {
    RUN <- TRUE
}
if(RUN){
    
    # load vars for year_i - n_yr_bioclimatic
    periods_i <- as.Date(paste0("01_01_",year_i),"%d_%m_%Y")
    periods_i <- format(
        seq.Date(from = periods_i-364, 
                 to = periods_i, 
                 by = "month"), 
        "%m_%Y")
    
    layers_i_pos <- unique(grep(paste(periods_i,collapse = "|"),all_layers))
    layers_i <- all_layers[layers_i_pos]
    
    tmax <- layers_i[grep("tmax",layers_i)]
    tmin <- layers_i[grep("tmin",layers_i)]
    prec <- layers_i[grep("prec",layers_i)]
    
    tmax <- rast(tmax)
    tmin <- rast(tmin)
    prec <- rast(prec)
    
    # create a SpatRasterDataset for tmean
    tmean <- sds(tmin,tmax)
    tmean <- terra::app(tmean, mean, cores = N_cores)
    
    tsea <- terra::app(tmean, sd, cores = N_cores)
    tmean <- terra::app(tmean, mean, cores = N_cores)
    tmax <- terra::app(tmax, max, cores = N_cores)
    tmin <- terra::app(tmin, min, cores = N_cores)
    
    prmean <- terra::app(prec, mean, cores = N_cores)
    prmax <- terra::app(prec, max, cores = N_cores)
    prmin <- terra::app(prec, min, cores = N_cores)
    prsea <- terra::app(prec, sd, cores = N_cores)
    
    bios_year_i <- c(tmean, tmax, tmin, tsea,
                     prmean, prmax, prmin, prsea)
    names(bios_year_i) <- c("mat","maxt","mint","seat",
                            "map","maxp","minp","seap")
    
    # mask the final raster
    bios_year_i <- terra::mask(bios_year_i,my_mask,maskvalues=0)
    
    # save raster
    writeRaster(bios_year_i, 
                filetosave,
                overwrite = TRUE)
    
}