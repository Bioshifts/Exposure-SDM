# Setup
rm(list=ls())
gc()

list.of.packages <- c("raster","rgdal","rgeos","terra","rgis","stars",
                      "SDMtune","dismo","spThin",
                      "tidyverse","tictoc",
                      "RSQLite","DBI","odbc",
                      "foreach","doParallel","pbapply","doFuture","future","doRNG")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

########################
# load functions
source(here::here("R/my_functions.R"))
source(here::here("R/Get_Bioclimatics.R"))
source(here::here("R/get_ecoregions.R"))
source(here::here("R/get_env_data_for_modeling.R"))

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
sptogo <- as.character(paste(command_args[1], collapse = " "))
sptogo <- gsub("_"," ",sptogo)
realm <- as.character(paste(command_args[2], collapse = " "))
print(sptogo)
print(realm)

# run test with terrestrial
# sptogo <- "Abax_parallelepipedus"
# sptogo <- gsub("_"," ",sptogo)
# realm = "Ter"

# run test with marine
# sptogo <- "Acropora_cervicornis"
# sptogo <- gsub("_"," ",sptogo)
# realm = "Mar"

########################
# source settings
source(here::here("R/settings.R"))

output_dir <- here::here("/media/seagate/boliveira/SDMs/Env_data",realm)

########################
# Check if file exists
file.test <- here::here(output_dir,paste(gsub(" ","_",sptogo),realm,"bio.RDS",sep="_"))
if(file.exists(file.test)){
    stop("No need to run this. It seems here have data for this species already!")
} 
# Check if dir for saving files exists
if(!file.exists(output_dir)){
    dir.create(output_dir)
}

########################
# Load occurrences

# get occurrences for species i
drv <- dbDriver("SQLite") 
my_db_file <- here::here("Data/BioshiftsExposure.sqlite")
db <- dbConnect(drv, my_db_file) 

sp_occ <- dbGetQuery(
    db, 
    paste0("SELECT year, month, cell, decimalLongitude, decimalLatitude    
           FROM occur 
           WHERE species = '", sptogo, "'"))
sp_occ$pa <- 1
names(sp_occ)[which(names(sp_occ) %in% c('decimalLongitude','decimalLatitude'))] <- c('x','y')

# get occurrences for Class i
classgoind <- readRDS(here::here("Data/SDMsSpList.RDS"))
classgoind <- unique(classgoind$Class[which(classgoind$species == gsub("_"," ",sptogo))])

if(realm == "Ter"){
    back_occ <- dbGetQuery(
        db, 
        paste0("SELECT year, month, cell, decimalLongitude, decimalLatitude    
           FROM ClassOccTer 
           WHERE species = '", classgoind, "'"))
}
if(realm == "Mar"){
    back_occ <- dbGetQuery(
        db, 
        paste0("SELECT year, month, cell, decimalLongitude, decimalLatitude    
           FROM ClassOccMar 
           WHERE species = '", classgoind, "'"))
}
back_occ$pa <- 0
names(back_occ)[which(names(back_occ) %in% c('decimalLongitude','decimalLatitude'))] <- c('x','y')

DBI::dbDisconnect(db)

########################
# Load ecoregions and more

BA <- get_ecoregions(realm,sp_occ)

back_coords <- vect(data.frame(back_occ), geom = c('x','y'))
sp_coords <- vect(data.frame(sp_occ), geom = c('x','y'))

# plot(BA)
# plot(back_coords,add=T)
# plot(sp_coords,add=T,col = "red")

########################
# Crop background data (Class level ocurrences) to the background area (ecoregions)
keep <- terra::extract(BA,data.frame(back_occ[,c('x','y')]))
keep <- which(keep$BA==1)
back_occ <- back_coords[keep]
back_occ <- as.data.frame(back_occ, geom = "XY")

# test
# plot(BA)
# plot(back_coords[keep],add=T)
# plot(vect(data.frame(sp_occ), geom = c('x','y')),add=T,col = "red")

########################
# combine species and background records
PresAbs <- rbind(sp_occ, back_occ)

########################
# Load environmental variables

# load climate variables
files <- sapply(myvars, function(i) {
    list.files(here::here(varsdir,i,"Tiles"), 
               full.names = T, recursive = T, pattern = "tile")
})
files <- unlist(files)

tiles <- paste0("tile",1:n_tiles,".tif")
tiles <- lapply(tiles, function(x) files[grep(x, files)])
# fix names
for(j in 1:length(tiles)){
    names(tiles[[j]]) <- sapply(tiles[[j]], function(i){
        t1 <- strsplit(i,"/")[[1]]
        t1 <- t1[length(t1)]
        cleanNamesFiles(t1, tile = T)
    })
}

########################
# clean occurrences
# remove dates outside the time range of the environmental data
PresAbs <- PresAbs |> filter(year > temporal_range_env_data[1] & year <= temporal_range_env_data[2])

########################
# Get env data from species occurrences
# We set a limit of 10000 background records as this has been show as enough to describe variation in environmental conditions. More records do not improve model performance: 
# https://onlinelibrary.wiley.com/doi/10.1111/j.0906-7590.2008.5203.x 
# https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1486

# plot(BA)
# plot(vect(data.frame(back_occ), geom = c('x','y')), col = "black", add=T)
# plot(vect(data.frame(sp_occ), geom = c('x','y')), col = "red", add=T)

# filter presences
sp_occ <- PresAbs %>% filter(pa == 1)
sp_occ <- if(nrow(sp_occ) > limit_recs) {sp_occ[sample(1:nrow(sp_occ), limit_recs),]} else {sp_occ}

# filter pseudo-absences
back_occ <- PresAbs %>% filter(pa == 0)
back_occ <- PresAbs %>% filter(!cell %in% sp_occ$cell) # Remove pseudo-absences at the same cells with presences
back_occ <- if(nrow(back_occ) > limit_recs) {back_occ[sample(1:nrow(back_occ), limit_recs),]} else {back_occ}

PresAbs <- rbind(sp_occ, back_occ)

table(PresAbs$pa)

# plot(BA)
# plot(vect(data.frame(back_occ), geom = c('x','y')), col = "black", add=T)
# plot(vect(data.frame(sp_occ), geom = c('x','y')), col = "red", add=T)
bioclimatics <- Get_Bioclimatics(sp_name=sptogo,
                                 occu=PresAbs,
                                 layers=tiles,
                                 min_year_layers = min_year_layers,
                                 realm = realm,
                                 dir_env_vars = varsdir,
                                 myvars = myvars,
                                 n_yr_bioclimatic = n_yr_bioclimatic,
                                 limit_occ_recs = NULL,
                                 output_dir = output_dir, 
                                 map_sps_occ = "/media/seagate/boliveira/SDMs/Maps_occs", 
                                 check_if_file_exists = TRUE,
                                 return_data = FALSE)

# table(bioclimatics$pa)

########################
