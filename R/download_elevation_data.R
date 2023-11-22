### Get elevation data
# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","MultiscaleDTM")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

########################
# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    work_dir <- getwd()
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
    work_dir <- getwd()
}

# source settings
source("R/settings.R")
source("R/velocity_functions.R")

# get elevation layer
if(!file.exists(here::here(work_dir,"Data/elevation_1km_Eckt.tif"))){
    download.file("https://data.earthenv.org/topography/elevation_1KMmn_GMTEDmn.tif",
                  here::here(work_dir,"Data/elevation_1km.tif"))
    elevation <- terra::rast(here::here(work_dir,"Data/elevation_1km.tif"))
    elevation <- terra::project(elevation,Eckt,filename=here::here(work_dir,"Data/elevation_1km_Eckt.tif"))
} else {
    elevation <- terra::rast(here::here(work_dir,"Data/elevation_1km_Eckt.tif"))
}

# calculate slope angle
slp_asp<- SlpAsp(r = elevation, w = c(5,5), unit = "degrees",
                 method = "queen", metrics = "slope",
                 filename=here::here(work_dir,"Data/elevation_slp_1km_Eckt.tif"))

elevation2 <- terra::rast(here::here(work_dir,"Data/elevation_1km.tif"))

slp_asp2<- terra::project(slp_asp,elevation2,
                          filename=here::here(work_dir,"Data/elevation_slp_1km.tif"))

plot(slp_asp)
dev.off()
