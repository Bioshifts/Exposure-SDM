
list.of.packages <- c("terra","here",
                      "reshape2","tidyr","tidyverse",
                      "foreach","doParallel","pbapply")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

# increase default timeout (=60s) prevents downloads to stop before finishing
options(timeout=3000) 

########################
# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}
work_dir <- getwd()

# source settings
source("R/settings.R")

# load functions and settings
source("R/getChelsa.R")
source("R/my_functions.R")

# N cpus to use
N_cpus <- parallelly::availableCores()

# settings
realm = "Ter"
my_res_down = "25km"

# path to save rasters
nc_path <- here::here(bios_dir(realm = realm))
if(!dir.exists(here::here(nc_path))){
    dir.create(here::here(nc_path), recursive = TRUE)
}

vars_dir_1km <- bios_dir(realm = realm)

model_raster_marine <- vars_dir(realm = "Mar")
model_raster_marine <- terra::rast(here(model_raster_marine,"model_raster_mar.tif"))

all_vars <- list.files(vars_dir_1km)

parallel::mclapply(all_vars, function(x){
    
    var_name = x
    var_address = here(vars_dir_1km,var_name)
    var_x <- terra::rast(here(vars_dir_1km,var_name))
    new_var_address = gsub("1km",my_res_down,var_address)
    var_new <- try(terra::rast(new_var_address), silent = TRUE)
    if(class(var_new)=="try-error"){
        terra::project(var_x, 
                       model_raster_marine, 
                       threads=TRUE, use_gdal=TRUE, gdal=TRUE,
                       filename = new_var_address,
                       overwrite=TRUE)
    }
})

