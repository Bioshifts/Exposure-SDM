###########
# Calculate bioclimatics from monthly raster data using dismo

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
}

# load functions
source("R/getChelsa.R")
source("R/my_functions.R")
source("R/settings.R")

# set n years to calculate bioclimatics
n_yr_bioclimatic

# set realm
realm = "Mar"

# raw files are saved here
scratch_dir <- here::here(vars_dir(realm),"SST")

# get directory to save bioclimatics
vars_dir <- vars_dir(realm)

# get vars
my_vars <- myvars(realm)

# create dir to store bioclimatics for each period
if(!dir.exists(here::here(vars_dir,"bio_proj"))){
    dir.create(here::here(vars_dir,"bio_proj"))
}

# load mask
my_mask <- terra::rast(here::here(vars_dir,"model_raster_mar.tif"))


# all layers
all_layers_names <- list.files(scratch_dir, pattern = ".tif")
all_layers_names <- gsub('.tif','',all_layers_names)
all_layers <- list.files(scratch_dir, pattern = ".tif",full.names = TRUE)

# possible years
my_yrs <- sapply(all_layers_names, function(x) strsplit(x,"_")[[1]][3])
my_yrs <- sort(as.numeric(unique(my_yrs)))

# remove the last - n_yr_bioclimatic (cannot calculate bioclimate for the last year)
my_yrs <- my_yrs[-1:-n_yr_bioclimatic]


# loop though possible years
# This function cant be serialized as SpatRaster cannot be sent to computer nodes. A solution is submitting jobs for each year
for(i in 1:length(my_yrs)) { 
    
    
    year_i <- my_yrs[i]
    
    cat("\ryear",i,"from",length(my_yrs), ".... nyear", year_i)
    
    # file name to save
    filetosave = here::here(vars_dir,"bio_proj",paste0("bios_",realm,"_",year_i,".tif"))
    
    if(!file.exists(filetosave)){
        # load vars for year_i - n_yr_bioclimatic
        periods_i <- as.Date(paste0("01_01_",year_i),"%d_%m_%Y")
        periods_i <- format(
            seq.Date(from = periods_i-364, 
                     to = periods_i, 
                     by = "month"), 
            "%m_%Y")
        
        layers_i_pos <- unique(grep(paste(periods_i,collapse = "|"),all_layers))
        layers_i <- all_layers[layers_i_pos]
        layers_i <- rast(layers_i)
        
        # layer names
        layers_names_i <- all_layers_names[layers_i_pos]
        names(layers_i) <- layers_names_i
        
        # get bioclimatics for year i
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