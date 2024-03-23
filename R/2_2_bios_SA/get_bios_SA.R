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
# polygontogo <- "A187_P2" # Ter

print(polygontogo)

cat("\rrunning polygon", polygontogo)

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
    dir.create(here::here(bios_SA_dir(realm),polygontogo),recursive = TRUE)
}

# get time period
if(grepl("A",polygontogo)){ # if v1
    biov <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v1))   
}
if(grepl("B",polygontogo)){ # if v1
    biov <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v2))
    biov$ID <- paste0("B",biov$Paper.ID,"_",biov$Study.Period)
    biov$START <- biov$Start.Year
    biov$END <- biov$End.Year   
}
period <- biov %>% 
    filter(ID == polygontogo) %>%
    dplyr::select(START,END) %>%
    round(digits = 0) %>%
    unique 
period = min(period$START):max(period$END)

# select bioclimatics for period
bioclima <- list.files(bios_dir(realm))
bioclima <- bioclima[grep(paste0(period,collapse = "|"),bioclima)]
bioclima <- data.frame(ID = polygontogo, bio = bioclima, year = period)

lapply(1:nrow(bioclima), function(i){
    # load in bioclimatics
    tmp <- terra::rast(here::here(bios_dir(realm), bioclima$bio[i]))
    tmp_name <- paste(bioclima$ID[i], "bios", bioclima$year[i])
    # check if file exists
    my_file <- here::here(bios_SA_dir(realm),polygontogo,paste0(tmp_name,".tif"))
    tmp_try <- try(rast(my_file), silent = TRUE)
    if(class(tmp_try)=='try-error'){ 
        # mask to the SA
        terra::window(tmp) <- terra::ext(SA_poly_i)
        terra::mask(tmp, mask = SA_poly_i, filename = my_file)
        terra::window(tmp) <- NULL
    }
})
