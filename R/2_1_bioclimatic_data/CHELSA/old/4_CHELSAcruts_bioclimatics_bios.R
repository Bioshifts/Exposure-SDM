###########
# Calculate bioclimatics from monthly raster


library(dismo)
library(raster)
library(terra)
library(parallel)
library(pbapply)
library(qs)

# set computer
computer = "muse"

if(computer == "muse"){
    wd <- "/storage/simple/projects/t_cesab/brunno/Exposure-SDM"
    setwd(wd)
    scratch_dir <- "/lustre/oliveirab"
}

# raw files are saved here
scratch_dir <- here::here(scratch_dir,"cruts")

# load functions
source("R/getChelsa.R")
source("R/my_functions.R")
source("R/settings.R")

# set n years to calculate bioclimatics
n_yr_bioclimatic

# set realm
realm = "Ter"

# get directory to save bioclimatics
vars_dir <- get_varsdir(computer = computer,
                        realm = realm)

# get vars
my_vars <- get_myvars(realm)


# create folder to store bios
dir_bios <- here::here(vars_dir,"bios")
if(!dir.exists(dir_bios)){
    dir.create(dir_bios, recursive = TRUE)
}


# all layers
all_layers_names <- list.files(scratch_dir, pattern = ".tif", recursive = TRUE)
all_layers_names <- gsub('.tif','',all_layers_names)
all_layers_names <- sapply(all_layers_names, function(x) strsplit(x,"/")[[1]][2])

all_layers <- list.files(scratch_dir, pattern = ".tif",full.names = TRUE, recursive = TRUE)

# possible periods
my_periods <- sapply(all_layers_names, function(x){
    tmp=strsplit(x,'_')[[1]]
    paste(tmp[2],tmp[3],sep = '_')
})
my_periods <- as.Date(paste0("01_",my_periods),"%d_%m_%Y")

# remove the first year + n_yr_bioclimatic (cannot calculate bioclimate for the last year)
my_periods <- my_periods[-grep(get_temporal_range_env_data(realm)[1], my_periods)]
my_periods <- format(my_periods,"%m_%Y")

# Run to period i
# load period to go
command_args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(command_args[1])


# This function cant be serialized as SpatRaster cannot be sent to computer nodes. A solution is submitting jobs for each year


cat("\rperiod",i,"from",length(my_periods))

# load vars for period_i - n_yr_bioclimatic
period_i <- my_periods[i]
cat("Running period", period_i)

# file name to save
filetosave = here::here(
    vars_dir,"bios",paste0("bios_",realm,"_",period_i,".tif"))

# check if files exists
if(!file.exists(filetosave)){
    
    # test if there is any issue loading the raster files
    tmp <- try(rast(filetosave),silent = TRUE)
    
    # if there was an issue, calculate bioclimatics again...
    if(class(tmp)=='try-error'){ 
        
        periods_i <- as.Date(paste0("01_",period_i),"%d_%m_%Y")
        periods_i <- format(
            seq.Date(from = periods_i-365, 
                     to = periods_i, 
                     by = "month"), 
            "%m_%Y")[1:12]
        
        layers_i_pos <- grep(paste(periods_i,collapse = "|"),all_layers)
        layers_i <- all_layers[layers_i_pos]
        
        tmax <- layers_i[grep("tmax",layers_i)]
        tmin <- layers_i[grep("tmin",layers_i)]
        prec <- layers_i[grep("prec",layers_i)]
        
        tmax <- rast(tmax)
        tmin <- rast(tmin)
        prec <- rast(prec)
        
        tmean <- lapply(1:nlyr(tmax), function(x){
            (tmin[[x]]+tmax[[x]])/2
        })
        tmean <- rast(tmean)
        
        tmean <- terra::app(tmean, mean)
        tmax <- terra::app(tmax, max)
        tmin <- terra::app(tmin, min)
        tsea <- terra::app(tmean, sd)
        
        prmean <- terra::app(prec, mean)
        prmax <- terra::app(prec, max)
        prmin <- terra::app(prec, min)
        prsea <- terra::app(prec, sd)
        
        bios_year_i <- c(tmean, tmax, tmin, tsea,
                         prmean, prmax, prmin, prsea)
        names(bios_year_i) <- c("mat","maxt","mint","seat",
                                "map","maxp","minp","seap")
        
        # save raster
        writeRaster(bios_year_i, 
                    filetosave,
                    overwrite = TRUE)
        
    }
    
}

