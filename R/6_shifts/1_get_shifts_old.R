# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","elevatr","raster","sf","rgdal","Hmisc","dplyr", "data.table","enmSdmX","biomod2","shape")
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
# source functions
source("R/my_functions.R")
source("R/bioshiftsFunction.R")
# source settings
source("R/settings.R")

check_if_file_exists <- FALSE

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
sptogo <- as.character(paste(command_args[1], collapse = " "))
realm <- as.character(paste(command_args[2], collapse = " "))

print(sptogo)
print(realm)

# sptogo="Anoplopoma_fimbria" # example of species with multiple shifts
# sptogo="Aplodactylus_arctidens"
# realm <- "Mar"

# sptogo="Amphipoea_oculea"
# sptogo="Myrmus_miriformis"
# sptogo="Carex_extensa"
# realm <- "Ter"

sptogo <- gsub(" ","_",sptogo)

shift_dir <- shift_dir(realm)
shiftplot_dir <- shiftplot_dir(realm)
sdm_dir <- sdm_dir(realm)

# create dirs
if(!dir.exists(shift_dir)){
    dir.create(shift_dir,recursive = TRUE)
}
if(!dir.exists(shiftplot_dir)){
    dir.create(shiftplot_dir,recursive = TRUE)
}

########################
# check if file exists
if(check_if_file_exists){
    test <- file.exists(
        here::here(shift_dir,
                   paste(sptogo,"shift SA.csv")))
    if(test){
        stop("Shift already calculated....Dont need to run for this species!") 
    }
}

########################
# Load files

# bioshift data
bioshift_info <- read.csv(here::here(sdm_dir,sptogo,paste(sptogo,"bioshift.csv")))
# only LAT
bioshift_info <- bioshift_info[which(bioshift_info$Type=="LAT"),]

# shift info
shift_info <- read.csv(here::here(sdm_dir,sptogo,paste(sptogo,"shift_info.csv")))
Time_periods = shift_info$time_period

# range shift SA
sdms <- list.files(here::here(sdm_dir,sptogo,gsub("_",".",sptogo)), 
                   pattern = "SA", 
                   full.names = TRUE)
sdms_names <- list.files(here::here(sdm_dir,sptogo,gsub("_",".",sptogo)), 
                         pattern = "SA")
# fix names
sdms_names <- gsub("proj_","",sdms_names)

# ensemble models
pos <- grep(" ens",sdms)
sdms_ens <- sdms[pos]
sdms_ens_names <- sdms_names[pos]

# single models
pos <- !grepl(" ens",sdms)
sdms <- sdms[pos]
sdms_names <- sdms_names[pos]

# check if have projection for all years
test_sdms_names <- sapply(1:nrow(shift_info), function(i){
    tmp <- sdms_ens_names[grep(shift_info$ID[i],sdms_ens_names)]
    all(sapply(round(shift_info$START[i],0):round(shift_info$END[i],0), function(x){
        any(grepl(x,tmp))
    }))
})

# filter only shifts with projections for all years
shift_info <- shift_info[test_sdms_names,]
Time_periods <- Time_periods[test_sdms_names]

sdms <- unlist(lapply(shift_info$ID, function(x) {
    sdms[grep(x,sdms)]
}))
sdms_names <- unlist(lapply(shift_info$ID, function(x) {
    sdms_names[grep(x,sdms_names)]
}))
sdms_ens <- unlist(lapply(shift_info$ID, function(x) {
    sdms_ens[grep(x,sdms_ens)]
}))
sdms_ens_names <- unlist(lapply(shift_info$ID, function(x) {
    sdms_ens_names[grep(x,sdms_ens_names)]
}))

# shift type
shift_info <- unique(merge(shift_info,bioshift_info[,c("Type","ID")],by="ID"))

#####################################
# load SDM ensemble projections for each shift in bioshifts database

sdms_ens_SA <- lapply(1:nrow(shift_info), function(i) {
    
    tmp <- sdms_ens[grep(shift_info$ID[i],sdms_ens)]
    
    tmp_tif = list.files(tmp,full.names = T,pattern = ".tif")
    
    tmp <- lapply(tmp_tif, function(x) mean(rast(x)))
    names(tmp) <- round(shift_info$START[i],0):round(shift_info$END[i],0)
    return(rast(tmp))
})
names(sdms_ens_SA) <- shift_info$ID


#####################################
# Calc SA range shifts ensemble

shifts_SA_ens <- list()

for(i in 1:length(sdms_ens_SA)){ # loop across shifts
    
    test <- try({
        # shift info
        info_i <- shift_info[i,]
        
        # sdms
        sdms_i <- sdms_ens_SA[[i]]
        
        # project to equal-area
        sdms_i <- terra::project(sdms_i,Eckt)
        
        # range 0-1
        sdms_i <- rast(lapply(sdms_i, range01raster))
        
        # quantiles
        quants <- c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99)
        
        names(sdms_i) <- paste0("t",1:nlyr(sdms_i))
        times <- round(shift_info$START[i],0):round(shift_info$END[i],0)
        
        # remove if all cells are unsuitable
        rem <- data.frame(minmax(sdms_i))
        rem <- apply(rem,2,function(x) all(is.na(x)))
        if(any(rem)){
            sdms_i <- sdms_i[[!rem]]
            times <- times[!rem]
        }
        
        tmp <- bioshifts(
            x = sdms_i,
            times = times,
            quants = quants,
            cores = 2,
            metrics = c("centroid","nsCentroid","nsQuants")
        )
        tmp <- tmp[,-2]
        names(tmp)[1:3] <- c("START","END","DUR")
        # add info
        tmp$ID <- info_i$ID
        tmp$Type <- info_i$Type
        tmp$Species <- info_i$Species
        tmp$time_period = info_i$time_period
        tmp$model = "ENSEMBLE"
        
        shifts_SA_ens[[i]] <- tmp
        
        pdf(here::here(shiftplot_dir,
                       paste(sptogo, info_i$ID, info_i$Type, info_i$time_period, tmp$model[1], "shiftplot SA.pdf")))
        shift_plot(r1 = sdms_i[[1]], 
                   r2 = sdms_i[[nlyr(sdms_i)]], 
                   quants = quants,
                   times = c(tmp$START[1], tmp$END[nrow(tmp)]))
        dev.off()
    },silent = TRUE)
    if(class(test)=="try-error"){
        shifts_SA_ens[[i]] <- NULL
    }
}


shifts_SA_ens <- rbindlist(shifts_SA_ens)

if(nrow(shifts_SA_ens)>0){
    write.csv(shifts_SA_ens, 
              here::here(shift_dir,
                         paste(sptogo,"shift ens SA.csv")),
              row.names = FALSE)
    
}

#####################################
# Calc SA range shifts single algorithms
sdm_algo <- c("GLM","GAM","GBM","MAXNET")

# load SDM single algorithms projections for each shift in bioshifts database
sdms_SA <- lapply(1:nrow(shift_info), function(i) {
    
    tmp <- sdms[grep(shift_info$ID[i],sdms)]
    
    tmp_tif = list.files(tmp,full.names = T,pattern = ".tif")
    tmp_tif <- tmp_tif[-grep("ClampingMask",tmp_tif)]
    
    tmp <- lapply(tmp_tif, function(x) {
        tmp_years <- strsplit(x, "/")[[1]]
        tmp_years <- tmp_years[length(tmp_years)]
        tmp_years <- strsplit(tmp_years, " ")[[1]]
        tmp_years <- tmp_years[3]
        
        tmp_tif_i <- rast(x)
        tmp_tif_i <- lapply(sdm_algo, function(i){
            mean(tmp_tif_i[[grep(i,names(tmp_tif_i))]])
        })
        tmp_tif_i <- rast(tmp_tif_i)
        names(tmp_tif_i) <- paste(sdm_algo, tmp_years)
        return(tmp_tif_i)
    })
    
    return(rast(tmp))
})
names(sdms_SA) <- shift_info$ID

shifts_SA <- list()

for(i in 1:length(sdms_SA)){ # loop across shifts
    
    test <- try({
        # shift info
        info_i <- shift_info[i,]
        
        # sdms
        sdms_i <- sdms_SA[[i]]
        
        shifts_SA_algo <- list()
        # run for each algorithm
        for(j in 1:length(sdm_algo)){ # loop across shifts
            
            algo_j <- sdm_algo[j]
            
            sdms_i_j <- sdms_i[[grep(algo_j,names(sdms_i))]]
            
            # project to equal-area
            sdms_i_j <- terra::project(sdms_i_j,Eckt)
            
            # range 0-1
            sdms_i_j <- rast(lapply(sdms_i_j, range01raster))
            
            # quantiles
            quants <- c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99)
            
            names(sdms_i_j) <- paste0("t",1:nlyr(sdms_i_j))
            times <- round(shift_info$START[1],0):round(shift_info$END[1],0)
            
            # remove if all cells are unsuitable
            rem <- data.frame(minmax(sdms_i_j))
            rem <- apply(rem,2,function(x) all(is.na(x)))
            if(any(rem)){
                sdms_i_j <- sdms_i_j[[!rem]]
                times <- times[!rem]
            }
            
            tmp <- bioshifts(
                x = sdms_i_j,
                times = times,
                quants = quants,
                cores = 2,
                metrics = c("centroid","nsCentroid","nsQuants")
            )
            tmp <- tmp[,-2]
            names(tmp)[1:3] <- c("START","END","DUR")
            # add info
            tmp$ID <- info_i$ID
            tmp$Type <- info_i$Type
            tmp$Species <- info_i$Species
            tmp$time_period = info_i$time_period
            tmp$model = "ENSEMBLE"
            tmp$algo = algo_j
            
            shifts_SA_algo[[j]] <- tmp
        }
        
    },silent = TRUE)
    if(class(test)=="try-error"){
        shifts_SA_algo[[j]] <- NULL
    }
    
    
    shifts_SA[[i]] <- rbindlist(shifts_SA_algo)
    
}

shifts_SA <- rbindlist(shifts_SA)

if(nrow(shifts_SA)>0){
    write.csv(shifts_SA, 
              here::here(shift_dir,
                         paste(sptogo,"shift SA.csv")),
              row.names = FALSE)
    
}