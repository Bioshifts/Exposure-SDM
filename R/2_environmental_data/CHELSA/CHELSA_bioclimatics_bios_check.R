
rm(list=ls())
gc()

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
vars_dir <- get_varsdir(realm = realm)

# get vars
my_vars <- get_myvars(realm)


# create folder to store bios
dir_bios <- here::here(vars_dir,"bios")
if(!dir.exists(dir_bios)){
    dir.create(dir_bios, recursive = TRUE)
}


# get mask
my_mask <- rast(here::here(vars_dir,"model_raster_ter_5km.tif"))
cellsNA <- which(my_mask[]==0)

# all layers
all_layers_names <- list.files(scratch_dir, pattern = ".tif", recursive = TRUE)
all_layers_names <- gsub('.tif','',all_layers_names)
all_layers_names <- sapply(all_layers_names, function(x) strsplit(x,"/")[[1]][2])

all_layers <- list.files(scratch_dir, pattern = ".tif",full.names = TRUE, recursive = TRUE)

# possible years
my_yrs <- sapply(all_layers_names, function(x) strsplit(x,"_")[[1]][3])
my_yrs <- sort(as.numeric(unique(my_yrs)))

# remove the last - n_yr_bioclimatic (cannot calculate bioclimate for the last year)
my_yrs <- my_yrs[-1:-n_yr_bioclimatic]



# check if all bioclimatics calculations went well
for(i in 1:length(my_yrs)) { 
    
    period_i <- my_yrs[i]   
    
    cat("\ryear",i,"from",length(my_yrs),"..... year", period_i)
    
    # file name to save
    filetosave = here::here(vars_dir,"bio_proj",paste0("bios_",realm,"_",my_yrs[i],".tif"))
    
    # test if there is any issue loading the raster files
    tmp <- try(rast(filetosave),silent = TRUE)
    
    # if there was an issue, calculate bioclimatics again...
    if(class(tmp)=='try-error'){ 
        
        cat('calculating again for period', period_i)
        
        # vector of dates
        periods_i <- as.Date(paste0("01_",period_i),"%d_%m_%Y")
        periods_i <- format(
            seq.Date(from = periods_i-365, 
                     to = periods_i, 
                     by = "month"), 
            "%m_%Y")
        
        layers_i_pos <- unique(grep(paste(periods_i,collapse = "|"),all_layers))
        layers_i <- all_layers[layers_i_pos]
        layers_i <- rast(layers_i)
        
        # layer names
        layers_names_i <- all_layers_names[layers_i_pos]
        names(layers_i) <- layers_names_i
        
        # get bioclimatics for period i
        mean_i <- terra::app(layers_i, mean)
        max_i <- terra::app(layers_i, max)
        min_i <- terra::app(layers_i, min)
        sd_i <- terra::app(layers_i, sd)
        bios_year_i <- c(mean_i, max_i, min_i, sd_i)
        
        # save raster
        writeRaster(bios_year_i, 
                    filetosave,
                    overwrite = TRUE)
        
        
    }
    
}