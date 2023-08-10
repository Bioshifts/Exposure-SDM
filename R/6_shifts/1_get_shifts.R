# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","elevatr","raster","sf","rgdal","Hmisc","dplyr", "data.table","enmSdmX")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)


computer = "muse"

# Eckert 4 equal-area projection
Eckt <- enmSdmX::getCRS('Eckert 4')

check_if_file_exists <- FALSE

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
sptogo <- as.character(paste(command_args[1], collapse = " "))
realm <- as.character(paste(command_args[2], collapse = " "))

print(sptogo)
print(realm)

# sptogo="Anoplopoma_fimbria" # example of species with multiple shifts
# sptogo="Abra_alba"
# sptogo="Deania_profundorum"
# realm <- "Mar"

# sptogo="Amphipoea_oculea"
# sptogo="Myrmus_miriformis"
# sptogo="Allophyes_oxyacanthae"
# realm <- "Ter"

sptogo <- gsub(" ","_",sptogo)

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

work_dir <- getwd()

shift_dir <- here::here(work_dir,"Data/SHIFT",realm)
shiftplot_dir <- here::here(work_dir,"Data/SHIFTplot",realm,sptogo)
scratch_dir <- "/lustre/oliveirab"
sdm_dir <- here::here(scratch_dir,"SDMs",realm)
script.dir <- here::here(work_dir,"R/6_shifts")

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
# load functions
source("R/my_functions.R")
source("R/bioshiftsFunction.R")
# source settings
source("R/settings.R")

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
pos <- grep("ens",sdms)
sdms_ens <- sdms[pos]
sdms_ens_names <- sdms_names[pos]
# single models
pos <- !grepl("ens",sdms)
sdms <- sdms[pos]
sdms_names <- sdms_names[pos]

# load SDM ensemble projections for each shift in bioshifts database
sdms_ens_SA <- lapply(1:nrow(shift_info), function(i) {
    
    time_tmp <- shift_info$time_period[i]
    time_tmp <- strsplit(time_tmp,"-")[[1]][1]
    
    tmp_start <- here::here(sdm_dir,sptogo,gsub("_",".",sptogo),
                            paste0("proj_",sptogo," ",shift_info$ID[i], " ",time_tmp," SA start ens"))
    tmp_start = list.files(tmp_start,full.names = T,pattern = ".out")
    tmp_start <- get(load(tmp_start))
    tmp_start = tmp_start@proj.out@val
    tmp_start = rast(tmp_start)
    
    tmp_end <- here::here(sdm_dir,sptogo,gsub("_",".",sptogo),
                          paste0("proj_",sptogo," ",shift_info$ID[i], " ",time_tmp," SA end ens"))
    tmp_end = list.files(tmp_end,full.names = T,pattern = ".out")
    tmp_end <- get(load(tmp_end))
    tmp_end = tmp_end@proj.out@val
    tmp_end = rast(tmp_end)
    
    # model names
    model_names <- names(tmp_start)
    model_names <- strsplit(model_names,"_")
    model_names <- sapply(model_names, function(x){
        paste(x[4])
    })
    
    names(tmp_start) <- model_names
    names(tmp_end) <- model_names
    
    list(start = tmp_start, end = tmp_end)
})
names(sdms_ens_SA) <- shift_info$ID
names(sdms_ens_SA[[1]]$start)

# load SDM projections (single models) for each shift in bioshifts database
make.unique.2 = function(x, sep='.'){
    ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}

sdms_SA <- lapply(1:nrow(shift_info), function(i) {
    
    time_tmp <- shift_info$time_period[i]
    time_tmp <- strsplit(time_tmp,"-")[[1]][1]
    
    tmp_start <- here::here(sdm_dir,sptogo,gsub("_",".",sptogo),
                            paste0("proj_",sptogo," ",shift_info$ID[i], " ",time_tmp," SA start"))
    tmp_start = list.files(tmp_start,full.names = T,pattern = ".out")
    tmp_start <- get(load(tmp_start))
    tmp_start = tmp_start@proj.out@val
    tmp_start = rast(tmp_start)
    
    tmp_end <- here::here(sdm_dir,sptogo,gsub("_",".",sptogo),
                          paste0("proj_",sptogo," ",shift_info$ID[i], " ",time_tmp," SA end"))
    tmp_end = list.files(tmp_end,full.names = T,pattern = ".out")
    tmp_end <- get(load(tmp_end))
    tmp_end = tmp_end@proj.out@val
    tmp_end = rast(tmp_end)
    
    # model names
    model_names <- names(tmp_start)
    model_names <- strsplit(model_names,"_")
    model_names <- sapply(model_names, function(x){
        paste(x[4])
    })
    model_names <- make.unique.2(model_names,sep = "_")
    
    names(tmp_start) <- model_names
    names(tmp_end) <- model_names
    
    list(start = tmp_start, end = tmp_end)
})
names(sdms_SA) <- shift_info$ID
names(sdms_SA)
names(sdms_SA[[1]]$start)

# # range shift BG
# sdms <- list.files(here::here(sdm_dir,sptogo,gsub("_",".",sptogo)), 
#                    pattern = "BG", 
#                    full.names = TRUE)
# sdms <- sdms[grep("ens",sdms)]
# 
# sdms <- lapply(sdms, function(x) {
#     tmp = list.files(x,full.names = T,pattern = ".out")
#     get(load(tmp))
# })
# sdms_BG <- lapply(sdms, function(x) {
#     tmp = x@proj.out@val
#     tmp = rast(tmp)
#     names(tmp) = paste(names(tmp), shift_info$ID)
#     return(tmp)
# })
# names(sdms_BG) <- c("end","start")

# shift type
shift_info <- unique(merge(shift_info,bioshift_info[,c("Type","ID")],by="ID"))
sel <- which(names(sdms_SA) %in% unique(shift_info$ID))
sdms_SA <- sdms_SA[sel]
sdms_ens_SA <- sdms_ens_SA[sel]

#####################################
# Calc SA range shifts ensemble

nmodels <- terra::nlyr(sdms_ens_SA[[1]]$start)

shifts_SA_ens <- list()

for(i in 1:length(sdms_ens_SA)){ # loop across shifts
    
    sh <- lapply(1:nmodels, function(j){ # loop across models
        
        test <- try({
            # shift info
            info_j <- shift_info[i,]
            
            # load in start and end periods
            r1 <- sdms_ens_SA[[i]]$start[[j]]
            r2 <- sdms_ens_SA[[i]]$end[[j]]
            
            # range 0-1
            r1 <- range01raster(r1)
            r2 <- range01raster(r2)
            
            # project to equal-area
            r1 <- terra::project(r1,Eckt)
            r2 <- terra::project(r2,Eckt)
            
            # quantiles
            quants <- c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)
            
            series <- c(
                r1,
                r2
            )
            names(series) <- c('t1', 't2')
            times <- c(round(shift_info$START[1],0), round(shift_info$END[1],0))
            
            tmp <- bioshifts(
                x = series,
                times = times,
                quants = quants,
                cores = 2,
                metrics = c("centroid","nsCentroid","nsQuants")
            )
            tmp <- tmp[,-2]
            names(tmp)[1:3] <- c("START","END","DUR")
            # add info
            tmp$ID <- info_j$ID
            tmp$Type <- info_j$Type
            tmp$Species <- info_j$Species
            tmp$time_period = info_j$time_period
            tmp$model = "ENSEMBLE"
            tmp$RUN = paste(names(sdms_ens_SA[[i]]$start[[j]]))
            
            pdf(here::here(shiftplot_dir,
                           paste(sptogo, info_j$ID, info_j$Type, info_j$time_period, tmp$model, tmp$RUN, "shiftplot SA.pdf")))
            shift_plot(r1, r2, times = c(tmp$START, tmp$END))
            dev.off()
            
            return(tmp)
            
        }, silent = TRUE)
        
        if(class(test)=="try-error"){
            return(NULL)
        } else {
            return(test)
        }
        
    })
    
    sh = rbindlist(sh)
    
    shifts_SA_ens[[i]] <- sh
}
shifts_SA_ens <- rbindlist(shifts_SA_ens)

if(nrow(shifts_SA_ens)>0){
    write.csv(shifts_SA_ens, 
              here::here(shift_dir,
                         paste(sptogo,"ens shift SA.csv")),
              row.names = FALSE)
    
}

#####################################
# Calc SA range shifts single algorithms

nmodels <- terra::nlyr(sdms_SA[[1]]$start)

shifts_SA <- list()

for(i in 1:length(sdms_SA)){ # across shifts
    
    sh <- lapply(1:nmodels, function(j){ # loop across models
        
        test <- try({
            
            # shift info
            info_j <- shift_info[i,]
            
            # load in start and end periods
            r1 <- sdms_ens_SA[[i]]$start[[j]]
            r2 <- sdms_ens_SA[[i]]$end[[j]]
            
            # range 0-1
            r1 <- range01raster(r1)
            r2 <- range01raster(r2)
            
            # project to equal-area
            r1 <- terra::project(r1,Eckt)
            r2 <- terra::project(r2,Eckt)
            
            # quantiles
            quants <- c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)
            
            series <- c(
                r1,
                r2
            )
            names(series) <- c('t1', 't2')
            times <- c(round(shift_info$START[1],0), round(shift_info$END[1],0))
            
            tmp <- bioshifts(
                x = series,
                times = times,
                quants = quants,
                cores = 2,
                metrics = c("centroid","nsCentroid","nsQuants")
            )
            tmp <- tmp[,-2]
            names(tmp)[1:3] <- c("START","END","DUR")
            # add info
            tmp$ID <- info_j$ID
            tmp$Type <- info_j$Type
            tmp$Species <- info_j$Species
            tmp$time_period = info_j$time_period
            
            tmp2 = paste(names(sdms_SA[[i]]$start[[j]]))
            tmp2 <- strsplit(tmp2,"_")[[1]]
            tmp$model <- tmp2[1]
            
            tmp$RUN = paste0("RUN_",tmp2[2])
            
            pdf(here::here(shiftplot_dir,
                           paste(sptogo, info_j$ID, info_j$Type, info_j$time_period, tmp$model, tmp$RUN, "shiftplot SA.pdf")))
            shift_plot(r1, r2, times = c(tmp$START, tmp$END))
            dev.off()
            
            return(tmp)
        }, silent = TRUE)
        
        if(class(test)=="try-error"){
            return(NULL)
        } else {
            return(test)
        }
        
    })
    
    sh = rbindlist(sh)
    
    shifts_SA[[i]] <- sh
}
shifts_SA <- rbindlist(shifts_SA)

if(nrow(shifts_SA)>0){
    
    write.csv(shifts_SA, 
              here::here(shift_dir,
                         paste(sptogo,"shift SA.csv")),
              row.names = FALSE)
}