library(dismo)
library(raster)
library(terra)
library(parallel)
library(pbapply)
library(qs)
library(data.table)

# set computer
computer = "muse"

if(computer == "muse"){
    wd <- "/storage/simple/projects/t_cesab/brunno/Exposure-SDM"
    setwd(wd)
    scratch_dir <- "/lustre/oliveirab"
}

# raw CHELSE data are save here
scratch_dir <- here::here(scratch_dir,"cruts")

# load functions
source("R/getChelsa.R")
source("R/my_functions.R")
source("R/settings.R")

# set n years to calculate bioclimatics
n_yr_bioclimatic

# set realm
realm = "Ter"

# get directory for monthly raster
vars_dir <- get_varsdir(computer = computer,
                        realm = realm)

# get vars
my_vars <- get_myvars(realm)

# all layers
all_layers_names <- list.files(scratch_dir, 
                               pattern = ".tif", 
                               recursive = TRUE)
rem <- grep("model_raster",all_layers_names)
if(any(rem)){
    all_layers_names <- all_layers_names[-rem]
}
all_layers_names <- gsub(".tif","",all_layers_names)
all_layers_names <- sapply(all_layers_names, function(x) strsplit(x,"/")[[1]][2])

all_layers <- list.files(scratch_dir, 
                         pattern = ".tif",full.names = TRUE, 
                         recursive = TRUE)
rem <- grep("model_raster",all_layers_names)
if(any(rem)){
    all_layers_names <- all_layers_names[-rem]
}

# possible years
my_yrs <- sapply(all_layers_names, function(x) strsplit(x,"_")[[1]][3])
my_yrs <- sort(as.numeric(unique(my_yrs)))

# remove the last - n_yr_bioclimatic (cannot calculate bioclimate for the last year)
my_yrs <- my_yrs[-1:-n_yr_bioclimatic]

# possible periods
my_periods <- sapply(all_layers_names, function(x){
    tmp=strsplit(x,'_')[[1]]
    paste(tmp[2],tmp[3],sep = '_')
})
my_periods <- as.Date(paste0("01_",my_periods),"%d_%m_%Y")

# remove the first year + n_yr_bioclimatic (cannot calculate bioclimate for the last year)
my_periods <- format(
    seq.Date(from = min(my_periods)+365, 
             to = max(my_periods), 
             by = "month"), 
    "%m_%Y")

# create temporary folder to store chunks
tmp_dir_chunks <- here::here(scratch_dir,"tmp_chunks")
if(!dir.exists(tmp_dir_chunks)){
    dir.create(tmp_dir_chunks)
}

# get mask var
my_mask <- terra::rast(here::here(vars_dir,"cruts","model_raster_ter_5km.tif"))
# get useful cells
cells_i <- terra::cells(my_mask, 1)[[1]]

# get chunks
chunk_size <- 10^5
cells_chunks <- split(cells_i, ceiling(seq_along(cells_i)/chunk_size))
length(cells_chunks)
# save chunks for fast reading
mclapply(1:length(cells_chunks),function(x){
    saveRDS(cells_chunks[[x]],here::here(tmp_dir_chunks, paste0("cell_chunk_",x)))
},mc.cores = 20)



