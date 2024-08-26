###########
# Calculate bioclimatics from monthly raster

# set computer
computer = "muse"

if(computer == "muse"){
    Sys.getenv('R_LIBS_SITE')
    .libPaths(Sys.getenv('R_LIBS_SITE'))
}

library(dismo)
library(raster)
library(terra)
library(parallel)
library(pbapply)
library(qs)
library(data.table)

if(computer == "muse"){
    wd <- "/storage/simple/projects/t_cesab/brunno/Exposure-SDM"
    setwd(wd)
    scratch_dir <- "/lustre/oliveirab"
}


# raw CHELSE data are save here
scratch_dir <- here::here(scratch_dir,"cruts")


# source code
source(here::here(wd,"R/my_functions.R"))
source(here::here(wd,"R/settings.R"))

# load year to go
command_args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(command_args[1])


# load vars for year_i - n_yr_bioclimatic
year_i <- my_yrs[i]
cat("Running year", year_i)

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


# create temporary folder to store bios from chunks
tmp_dir_bios <- here::here(scratch_dir,"tmp_bios",year_i)
if(!dir.exists(tmp_dir_bios)){
    dir.create(tmp_dir_bios, recursive = TRUE)
}

# create folder to store bios proj
dir_bios_proj <- here::here(vars_dir,"bio_proj")
if(!dir.exists(dir_bios_proj)){
    dir.create(dir_bios_proj, recursive = TRUE)
}

my_mask <- terra::rast(here::here(vars_dir,"cruts","model_raster_ter_5km.tif"))

# check if file exists
my_file <- here::here(vars_dir,"bio_proj",paste0("bios_",realm,"_",year_i,".tif"))

# do not run if file exists...
if(!file.exists(my_file)){
    
    ncores_slurm = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
    ncores_parallel = parallel::detectCores()
    cat("N cores requested in the job is", ncores_slurm)
    cat("detectCores returns", ncores_parallel)
    
    # force maximum cores
    ncores = max(ncores_slurm,ncores_parallel)
    cat("using", ncores,"\nnot sure if this works...")
    
    
    # loop trough cell chunks
    parallel::mclapply(1:length(cells_chunks), function(x){
        # get cells to go
        cells_i <- readRDS(here::here(tmp_dir_chunks, paste0("cell_chunk_",x)))
        
        # get bioclimatics from cells i
        data_env <- data.frame(layers_i[cells_i])
        tasmax <- data_env[,grep("tmax",names(data_env))]
        tasmin <- data_env[,grep("tmin",names(data_env))]
        prec <- data_env[,grep("prec",names(data_env))]
        prec[prec<0] <- 0
        
        # fix layers
        data_env <- data.frame(cbind(tasmax, tasmin, prec))
        names_data <- names(data_env)
        names_data <- gsub("prec","pr",names_data)
        names_data <- gsub("tmax","tasmax",names_data)
        names_data <- gsub("tmin","tasmin",names_data)
        
        names(data_env) <- names_data
        
        # CHELSAcruts are in C/10, but CHELSA are in K/10. 
        # Transform to C because bioclimatic variables are usually in C.
        # years_i <- sapply(names(data_env), function(j){
        #     as.numeric(strsplit(j,"_")[[1]][3])
        # })
        # 
        # if(any(years_i >= 1980)){ # if year CHELSA
        #     
        #     t_data <- which(years_before_1980 >= 1980 & grepl("tas",names(data_env)))
        #     data_env[,t_data] <- ((data_env[,t_data]/10) - 273.15) * 10
        #     
        # }
        
        bios_i <- bioclimatics_land(env_data = data_env)
        
        # save 
        bios_year_i_xy <- terra::xyFromCell(my_mask, cells_i)
        bios_year_i <- cbind(bios_i, bios_year_i_xy)
        coordinates(bios_year_i) <- ~ x + y
        # coerce to SpatialPixelsDataFrame
        gridded(bios_year_i) <- TRUE
        # coerce to raster
        bios_year_i <- rast(bios_year_i)
        # save raster
        writeRaster(bios_year_i, 
                    here::here(tmp_dir_bios,paste0("bios_yr_",year_i,"_chunk_",x,".tif")),
                    overwrite = TRUE)
        
    }, mc.cores = ncores)
    
    # load all years
    bios_year_i <- list.files(here::here(tmp_dir_bios), full.names = TRUE, pattern = '.tif')
    # list of rasters
    x <- lapply(bios_year_i, terra::rast)
    # coerce spatial raster collection
    rsrc <- sprc(x)
    # merge rasters
    m <- merge(rsrc)
    # save raster
    writeRaster(m, 
                here::here(vars_dir,"bio_proj",paste0("bios_",realm,"_",year_i,".tif")),
                overwrite = TRUE)
    
    # delete tmp files
    bios_year_i <- list.files(here::here(tmp_dir_bios), full.names = TRUE)
}

# delete tmp folder
unlink(here::here(tmp_dir_bios), recursive = TRUE)



# test <- rast(here::here(vars_dir,"bio_proj",paste0("bios_",realm,"_",year_i,".tif")))
# plot(test);dev.off()

