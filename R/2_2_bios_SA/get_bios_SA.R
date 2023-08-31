# --------------------------------------------------------
# title: "Extract bioclimatics at each study area in the Bioshifts database"
# author: "Brunno F Oliveira & Bioshifts group"
# --------------------------------------------------------

# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","raster","sf","rgdal","Hmisc","dplyr", "data.table","here")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
print(command_args)
polygontogo <- command_args
# polygontogo <- "A112_P1" # Mar # South
# polygontogo <- "A43_P1" # Mar # North
# polygontogo <- "A102_P2" # Ter

print(polygontogo)

cat("\rrunning polygon", polygontogo)

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

# source settings
source("R/settings.R")


# load SA polygon i
SA_poly_i <- terra::vect(here::here(SA_shps_dir,paste0(polygontogo,".shp")))

# Is it Terrestrial or Marine?
my_test <- if(any(is.na(SA_poly_i$EleExtentm))){ 
    realm = "Mar"
} else {
    realm = "Ter"
}

# create dir to store results
if(!dir.exists(here::here(bios_SA_dir(realm),polygontogo))){
    dir.create(here::here(bios_SA_dir(realm),polygontogo))
}

# get time period
if(grepl("A",polygontogo)){ # if v1
    biov <- read.csv(here::here(Bioshifts_dir,"biov1_fixednames.csv"))   
}
if(grepl("B",polygontogo)){ # if v1
    biov <- read.csv(here::here(Bioshifts_dir,"biov2_fixednames.csv"))
    biov$ID <- paste0("B",biov$Paper.ID,"_",biov$Study.Period)
    biov$START <- biov$Start.Year
    biov$END <- biov$End.Year   
}
period <- biov %>% 
    filter(ID == polygontogo) %>%
    select(START,END) %>%
    round(digits = 0) %>%
    unique 
period = period$START:period$END

# select bioclimatics
bioclima <- list.files(bios_dir(realm))
bioclima <- bioclima[grep(paste0(period,collapse = "|"),bioclima)]

# load in bioclimatics
bioclima <- lapply(bioclima, function(x){
    terra::rast(here::here(bios_dir(realm), x))
})

# mask to the SA
bioclima <- lapply(bioclima, function(x){
    terra::mask(terra::crop(x, SA_poly_i), SA_poly_i)
})

# save
sapply(1:length(bioclima), function(x){
    terra::writeRaster(bioclima[[x]], 
                       here::here(bios_SA_dir(realm),polygontogo,paste0(polygontogo," bios ",period[x],".tif")), 
                       overwrite = TRUE)
})

