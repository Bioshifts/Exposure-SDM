




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
sptogo <- gsub("_"," ",sptogo)
realm <- as.character(paste(command_args[2], collapse = " "))
print(sptogo)
print(realm)

# run test with terrestrial
# sptogo <- "Caltha_palustris"
# sptogo <- gsub("_"," ",sptogo)
# realm = "Ter"

# run test with marine
# sptogo <- "Zaprora_silenus"
# sptogo <- gsub("_"," ",sptogo)
# realm = "Mar"

N_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
# N_cores = 3

cat("N cores = ", N_cores)

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

########################
# load functions
source("R/my_functions.R")
source("R/Get_Bioclimatics.R")
source("R/get_ecoregions.R")
source("R/get_env_data_for_modeling.R")
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

sp_occ <- qs::qread(here::here(occ_dir, paste0(gsub(" ","_",sptogo), ".qs")))
sp_occ$pa <- 1
sp_occ <- sp_occ %>%
    mutate(x = decimalLongitude,
           y = decimalLatitude) %>%
    dplyr::select(-c('basisOfRecord','speciesKey','decimalLongitude','decimalLatitude')) 

nrow(sp_occ)

########################
# Load ecoregions
BA <- get_ecoregions(realm = realm, 
                     PresAbs = sp_occ, 
                     varsdir = vars_dir, 
                     mask.ras = mask.ras,
                     return.shp = FALSE,
                     return.raster = TRUE)
BA <- BA$raster_file
# plot(BA);dev.off()

########################
# Create random pseudo-absences

# 1) get cells at the BA
cells_BA <- terra::cells(BA, 1)[[1]]
# 2) get dates
env_range <- temporal_range_env_data(realm) 
env_range[1] <- env_range[1] + n_yr_bioclimatic # + n_yr_bioclimatic to be able to calculate bioclimatics for the previous n_yr_bioclimatic
all_dates <- format(
    seq.Date(from = as.Date(paste0(env_range[1],"/01/01")), 
             to = as.Date(paste0(env_range[2],"/12/01")), 
             by = "month"), 
    "%m_%Y")
# 3) sample random cells and dates
random_cells <- sample(cells_BA, size = nrow(sp_occ), replace = TRUE)
random_dates <- sample(all_dates, size = nrow(sp_occ), replace = TRUE)
PA_cell_dates <- paste(random_cells,random_dates)
# 4) check for the existence of any combination of cell dates from PA in sp_occ
occ_cell_dates <- paste(sp_occ$cell, paste(sp_occ$month,sp_occ$year,sep="_"))
while(any(PA_cell_dates %in% occ_cell_dates)){
    # resample
    random_cells <- sample(cells_BA, size = nrow(sp_occ), replace = TRUE)
    random_dates <- sample(all_dates, size = nrow(sp_occ), replace = TRUE)
    PA_cell_dates <- paste(random_cells,random_dates)
}
# Create a PA dataset
tmp <- lapply(random_dates, function(x) strsplit(x, "_")[[1]])
PA_years <- sapply(tmp, function(x) x[2])
PA_months <- sapply(tmp, function(x) x[1])

PA_occ <- data.frame(year = PA_years, 
                     month = PA_months,
                     species = sp_occ$species[1],
                     cell = random_cells,
                     pa = 0)
PA_xy <- terra::xyFromCell(BA,random_cells)    
PA_occ <- cbind(PA_occ, PA_xy) 
# # Test plot
# back_coords <- vect(data.frame(PA_occ), geom = c('x','y'))
# sp_coords <- vect(data.frame(sp_occ), geom = c('x','y'))
# plot(BA)
# plot(back_coords,add=T,cex=.5)
# plot(sp_coords,add=T,col = "red",cex=.5)
# dev.off()

# Add PA to occ
sp_occ <- rbind(sp_occ, PA_occ)

table(sp_occ$pa)
# head(sp_occ)
# min(sp_occ$year)

########################
# Load bioclimatics from each date and cell

# 1) Get all possible dates
sp_occ$date <- paste("01",sp_occ$month, sp_occ$year,sep = "_")
sp_occ$date <- as.Date(sp_occ$date,"%d_%m_%Y")
sp_occ$date <- format(sp_occ$date,"%m_%Y")

possibledates <- unique(sp_occ$date)

dim(sp_occ)
length(possibledates)

# 2) calculate bioclimatics for occurrences
if(realm == "Ter"){
    all_layers <- list.files(here::here(vars_dir,paste0("cruts_",my_res)), pattern = ".tif", full.names = TRUE, recursive = TRUE)
    all_layers_names <- list.files(here::here(vars_dir,paste0("cruts_",my_res)), pattern = ".tif", recursive = TRUE)
}
if(realm == "Mar"){
    all_layers <- list.files(here::here(vars_dir,"SST"), pattern = ".tif", full.names = TRUE, recursive = TRUE)
    all_layers_names <- list.files(here::here(vars_dir,"SST"), pattern = ".tif", recursive = TRUE)
    
}
all_layers_names <- gsub('.tif','',all_layers_names)
all_layers_names <- sapply(all_layers_names, function(x) strsplit(x,"/")[[1]][2])

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
            
            # subset occ data for date i
            sub_occ <- sp_occ %>% dplyr::filter(date == date_i)
            
            # get range of dates to get the bioclimatics
            periods_i <- paste0("01_",date_i)
            periods_i <- as.Date(periods_i,"%d_%m_%Y")
            periods_i <- format(
                seq.Date(from = periods_i-365, 
                         to = periods_i, 
                         by = "month"), 
                "%m_%Y")[1:12]
            
            # subset env layers for this period
            layers_i_pos <- unique(grep(paste(periods_i,collapse = "|"),all_layers))
            layers_i <- all_layers[layers_i_pos]
            layers_i <- rast(layers_i)
            # fix layers names
            layers_i_names <- all_layers_names[layers_i_pos]
            names(layers_i) <- layers_i_names
            
            # extract data
            layers_i <- layers_i[sub_occ$cell]
            
            # calculate bioclimatics
            if(realm == "Ter"){
                bios_i <- bioclimatics_land_simple(layers_i)
            }
            if(realm == "Mar"){
                bios_i <- bioclimatics_ocean_simple(layers_i)
            }
            
            bios_i$ID <- paste(sub_occ$cell,sub_occ$date,sep = "_")
            
            return(bios_i)
            
        }, 
        mc.cores = N_cores)
        
    }
    
    tmp <- data.table::rbindlist(tmp)
    
    biosclim <- rbind(biosclim, tmp)
    
}



# Compile
sp_occ <- merge(sp_occ, biosclim, by = "ID")
sp_occ <- sp_occ[,-"ID"]

table(sp_occ$pa)

# save
qs::qsave(sp_occ, here::here(env_data_dir,paste0(gsub(" ","_",sptogo),"_",realm,".qs")))

########################
# Stop the clock
end_time <- Sys.time()
end_time - start_time
