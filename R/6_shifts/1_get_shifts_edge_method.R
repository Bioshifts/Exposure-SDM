# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","Hmisc","dplyr", "tidyterra","data.table", "biomod2", "geosphere")
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

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
sptogo <- as.character(paste(command_args[1], collapse = " "))
realm <- as.character(paste(command_args[2], collapse = " "))
IDtogo <- as.character(paste(command_args[3], collapse = " "))

# sptogo="Acanthaluteres_spilomelanurus"
# realm <- "Mar"
# IDtogo = "A112_P1"

# sptogo="Veronica_montana"
# realm <- "Ter"
# IDtogo = "A19_P1"

# sptogo="Ardeotis_australis"
# realm <- "Ter"
# IDtogo = "A36_P1"

print(sptogo)
print(realm)
print(IDtogo)

########################
# source functions
source("R/my_functions.R")
source("R/range_shift_functions.R")
source("R/velocity_functions.R")
# source settings
source("R/settings.R")

# create dirs
if(!dir.exists(shift_dir(realm))){
    dir.create(shift_dir(realm),recursive = TRUE)
}

########################
# Load files

# shift info
shift_info <- read.csv(here::here(sdm_dir(realm),sptogo,paste(sptogo,"shift_info.csv")))
# select the shift of interest
shift_info <- shift_info %>% filter(ID == IDtogo)

# range shift SA
sdms <- list.files(here::here(sdm_dir(realm),sptogo,gsub("_",".",sptogo)), 
                   pattern = "SA ens", 
                   full.names = TRUE)
sdms_names <- list.files(here::here(sdm_dir(realm),sptogo,gsub("_",".",sptogo)), 
                         pattern = "SA ens")
# fix names
sdms_names <- gsub("proj_","",sdms_names)

# check if have projection for all years
test_sdms_names <- sapply(1:nrow(shift_info), function(i){
    tmp <- sdms_names[grep(shift_info$ID[i],sdms_names)]
    all(sapply(round(shift_info$Start[i],0):round(shift_info$End[i],0), function(x){
        any(grepl(x,tmp))
    }))
})

# filter only shifts with projections for all years
shift_info <- shift_info[test_sdms_names,]
time_perid <- round(min(shift_info$Start),0):round(max(shift_info$End),0)

sdms <- unlist(lapply(IDtogo, function(x) {
    sdms[grep(x,sdms)]
}))
sdms_names <- unlist(lapply(IDtogo, function(x) {
    sdms_names[grep(x,sdms_names)]
}))




#####################################
# load SDM ensemble projections for the shift of interest

tmp <- sdms[grep(IDtogo,sdms)]
tmp <- tmp[grep(paste(time_perid,collapse = "|"),tmp)]
tmp_tif = list.files(tmp, full.names = TRUE, pattern = "\\.tif")
tmp <- lapply(tmp_tif, function(x) {
    tmp <- rast(x)
    if(nlyr(tmp)>1){
        # get cross validation file
        CV_file <- get_CV_file(realm,sptogo,type="RUN")
        tmp_runs <- grep("RUN", unlist(strsplit(names(tmp), split = "_")), value = TRUE, ignore.case = FALSE)
        names(tmp) <- tmp_runs
        CV_file_runs <- CV_file$merged.by.run
        tmp <- tmp[CV_file_runs]
        CV_file <- filter(CV_file, merged.by.run == names(tmp))
        
        round(weighted.mean(tmp,CV_file$validation),0)
    } else {
        tmp
    }
})
names(tmp) <- time_perid
sdms_SA <- rast(tmp)

# plot(sdms_SA[[1]]);dev.off()
# plot(ifel(sdms_SA>cutoff,1,NA));dev.off()

#####################################
# Calc SA range shifts ensemble

# shift ID for saving objects
shift_ID <- paste(shift_info$Species, shift_info$ID, shift_info$time_period, realm, "SA", sep = "_")

# sdms
sdms_i <- sdms_SA

for(i in 1:nrow(shift_info)){
    
    shift_ID_togo <- shift_ID[i]
    periods_i <- as.character(round(min(shift_info$Start[i]),0):round(max(shift_info$End[i]),0))
    sdms_i_togo <- sdms_i[[periods_i]]
    
    # plot(sdms_i_togo[[1]]);dev.off()
    
    SA_i <- terra::vect(here(SA_shps_dir,paste0(IDtogo,".shp")))
    
    # calculate range shifts
    range_shifts <- range_shift(x = sdms_i_togo,
                                SA = SA_i,
                                periods = periods_i,
                                raster_size_tolerance = 10^5,
                                proj_equal = TRUE)
    
    ###
    # save results
    write.csv(range_shifts[[1]], 
              here::here(shift_dir(realm),
                         paste0(shift_ID_togo,".csv")),
              row.names = FALSE)
    
    write.csv(range_shifts[[2]], 
              here::here(shift_dir(realm),
                         paste0(shift_ID_togo,"_edges.csv")),
              row.names = FALSE)
    
}
