
########################
# Setup
rm(list=ls())
gc()

list.of.packages <- c("raster","rgdal","rgeos","terra","stars",
                      "tidyverse","tictoc","qs",
                      "foreach","doParallel","pbapply")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
sptogo <- as.character(paste(command_args[1], collapse = " "))
realm <- as.character(paste(command_args[2], collapse = " "))
print(sptogo)
print(realm)

# run test with terrestrial
# sptogo <- "Caltha_palustris"
# sptogo="Formica_sanguinea"
# realm = "Ter"

# run test with marine
# sptogo <- "Alosa_pseudoharengus"
# realm = "Mar"

N_cores = parallelly::availableCores()
# N_cores = 3

cat("N cores = ", N_cores)

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

########################
# load functions
source("R/my_functions.R")
source("R/get_ecoregions.R")
# source settings
source("R/settings.R")

env_data_dir <- env_data_dir(realm)

########################
# Check if dir exists for saving files exists
if(!dir.exists(env_data_dir)){
    dir.create(env_data_dir, recursive = TRUE)
}



vars_dir <- vars_dir(realm)
mask.ras = 
    if(realm=="Ter") { terra::rast(here::here(vars_dir,paste0("model_raster_ter_",my_res,".tif"))) 
    } else {
        if(realm=="Mar") { terra::rast(here::here(vars_dir,"model_raster_mar.tif")) 
        }
    }

########################
# Start the clock
start_time <- Sys.time()


########################
# Load occurrences for species i

sp_occ <- qs::qread(here::here(occ_dir, paste0(sptogo, ".qs")))
sp_occ$pa <- 1

sp_occ <- sp_occ %>%
    mutate(x = decimalLongitude,
           y = decimalLatitude) %>%
    dplyr::select(-c('basisOfRecord','speciesKey','decimalLongitude','decimalLatitude')) 

nrow(sp_occ)

########################
# Load ecoregions
output_dir <- here::here(sdm_dir(realm), sptogo)
if(!dir.exists(output_dir)){
    dir.create(output_dir,recursive = TRUE)
}

BA <- get_ecoregions(realm = realm, 
                     sptogo = sptogo,
                     PresAbs = data.frame(sp_occ[,c("x","y")]), 
                     varsdir = vars_dir, 
                     mask.ras = mask.ras,
                     return.shp = TRUE,
                     return.raster = FALSE,
                     check_if_exists = TRUE,
                     output_dir = output_dir)

BA <- BA$shape_file
BA <- rasterize(BA, mask.ras)

########################
# Create random pseudo-absences
env_range <- temporal_range_env_data(realm) 
PA_occ <- create_temporal_pseudo_absences(BA, env_range, sp_occ)

# # Test plot
# back_coords <- vect(data.frame(PA_occ), geom = c('x','y'))
# sp_coords <- vect(data.frame(sp_occ), geom = c('x','y'))
# plot(BA)
# plot(back_coords,add=T,cex=.5)
# plot(sp_coords,add=T,col = "red",cex=.2)
# dev.off()

# Add PA to occ
sp_occ <- rbind(sp_occ, PA_occ)

table(sp_occ$pa)
# head(sp_occ)
# min(sp_occ$year)

########################
# Load bioclimatics from each date and cell

# layers
if(realm == "Ter"){
    all_layers <- list.files(here::here(vars_dir,my_res), pattern = ".tif", full.names = TRUE, recursive = TRUE)
    all_layers_names <- list.files(here::here(vars_dir,my_res), pattern = ".tif", recursive = TRUE)
    all_layers_names <- gsub('.tif','',all_layers_names)
    all_layers_names <- sapply(all_layers_names, function(x) strsplit(x,"/")[[1]][2])
}
if(realm == "Mar"){
    all_layers <- list.files(here::here(vars_dir,"SST"), pattern = ".tif", full.names = TRUE, recursive = TRUE)
    all_layers_names <- list.files(here::here(vars_dir,"SST"), pattern = ".tif", recursive = TRUE)
    all_layers_names <- gsub('.tif','',all_layers_names)
}


# 1) Get all possible dates
sp_occ$date <- paste("01",sp_occ$month, sp_occ$year,sep = "_")
sp_occ$date <- as.Date(sp_occ$date,"%d_%m_%Y")
sp_occ$date <- format(sp_occ$date,"%m_%Y")

possibledates <- unique(sp_occ$date)

dim(sp_occ)
length(possibledates)

# 3) Create ID for merge bioclimatics with occ
sp_occ$ID <- paste(sp_occ$cell,sp_occ$date,sep = "_")

# 4) create chunks for possible dates >> avoid overload memory
N_dates = N_cores
chunks <- split(possibledates, ceiling(seq_along(1:length(possibledates))/N_dates))

biosclim <- data.frame()

for(j in 1:length(chunks)){
    
    cat("\nchunk", j, "from", length(chunks))
    
    possibledates_j <- chunks[[j]]
    
    tmp <- NA
    trying = 0
    
    while(trying == 0 | any(sapply(tmp,class)=="try-error")){
        
        trying = trying+1
        
        if(trying > 1){ 
            
            cat("\nretrying", trying) 
            
        }
        
        if(trying > 10){ 
            
            stop("\nToo many retries\nStop and report!") 
            
        }
        
        tmp <- parallel::mclapply(1:length(possibledates_j), function(i) {
            
            
            # date i
            date_i <- possibledates_j[i]
            bios_i <- bioclimatics_from_date(date_i,realm,sp_occ,n_yr_bioclimatic)
            
        }, 
        mc.cores = N_cores)
        
    }
    
    tmp <- data.table::rbindlist(tmp)
    biosclim <- rbind(biosclim, tmp)
    
}

# Compile
biosclim <- filter(biosclim, !duplicated(ID))
sp_occ <- merge(sp_occ, biosclim, by = "ID", all.x = TRUE)
sp_occ <- sp_occ[,-"ID"]

table(sp_occ$pa)

# save
qs::qsave(sp_occ, here::here(env_data_dir,paste0(sptogo,"_",realm,".qs")))

########################
# Stop the clock
end_time <- Sys.time()
end_time - start_time
