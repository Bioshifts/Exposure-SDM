# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","Hmisc","dplyr", "tidyterra","data.table")
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

# N cores
ncores <- parallelly::availableCores()


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
# sptogo="Abies alba"
# realm <- "Ter"

sptogo <- gsub(" ","_",sptogo)

########################
# source functions
source("R/my_functions.R")
source("R/bioshiftsFunction.R")
source("R/velocity_functions.R")
# source settings
source("R/settings.R")

check_if_file_exists <- FALSE

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
    
    test <- rast(tmp_tif[1])
    
    tmp <- lapply(tmp_tif, function(x) mean(rast(x)))
    names(tmp) <- round(shift_info$START[i],0):round(shift_info$END[i],0)
    return(rast(tmp))
})
names(sdms_ens_SA) <- shift_info$ID

# plot(sdms_ens_SA[[1]][[1]]);dev.off()

#####################################
# Calc SA range shifts ensemble

shifts_SA_ens <- list()

for(i in 1:length(sdms_ens_SA)){ # loop across shifts
    
    test <- try({
        
        # shift info
        info_i <- shift_info[i,]
        
        # shift ID for saving objects
        shift_ID <- paste(c("SA",info_i[,-2:-3],realm),collapse = "_")
        
        # sdms
        sdms_i <- sdms_ens_SA[[i]]
        
        # project to equal-area
        sdms_i <- terra::project(sdms_i,Eckt)
        
        #######
        # calculate the trend (C/year)
        cat("calculate the trend\n")
        
        ttrend_file <- here::here(tmp_dir, paste0("trend_",shift_ID,".tif"))
        ttrend <- try(terra::rast(ttrend_file),silent = TRUE)
        if(class(ttrend)=="try-error"){
            ttrend = temp_grad(
                sdms_i,
                th = 0.25*nlyr(sdms_i), ## set minimum N obs. to 1/4 time series length
                tempfile = ttrend_file,
                overwrite = TRUE,
                ncores = ncores)
        }
        
        #######
        # Get averaged climate layers
        avg_sdm_file <- here::here(tmp_dir, paste0("AvgSDM_",shift_ID,".tif"))
        avg_sdms <- try(terra::rast(avg_sdm_file),silent = TRUE)
        if(class(avg_sdms)=="try-error"){
            avg_sdms <- terra::app(
                sdms_i, 
                fun = mean, na.rm = TRUE, 
                filename=avg_sdm_file,
                overwrite=TRUE,
                cores = ncores)
        }
        
        
        #######
        # Get the spatial gradient (C/km)
        cat("calculate the spatial gradient\n")
        
        spgrad_file <- here::here(tmp_dir, paste0("spgrad_",shift_ID,".qs"))
        spgrad <- try(qs::qread(spgrad_file, nthreads = ncores),silent = TRUE)
        if(any(class(spgrad)=="try-error")){
            
            if(realm=="Ter"){
                
                # load SA polygon
                SA_i <- terra::vect(here::here(SA_shps_dir,paste0(shift_info$ID[i],".shp")))
                
                # Big area test
                area_i <- SA_i$Areakm2[1]
                big <- area_i > 10^4
                
                if(big){
                    avg_sdms_tiles <- avg_sdms
                    # define tile resolution 
                    nrow(avg_sdms_tiles) <- round(nrow(avg_sdms)/100,0)
                    ncol(avg_sdms_tiles) <- round(ncol(avg_sdms)/100,0)
                    avg_sdms_tiles[] <- 1:ncell(avg_sdms_tiles)
                    # mask raster
                    SA_i_Eck <- terra::project(SA_i,avg_sdms_tiles)
                    avg_sdms_tiles <- terra::mask(avg_sdms_tiles, SA_i_Eck)
                    
                    # plot(avg_sdms_tiles);plot(SA_i_Eck,add=T);dev.off()
                    
                    parallel::mclapply(terra::cells(avg_sdms_tiles), function(x) {
                        tmp_file <- here::here(tmp_dir, paste0("spgrad_tile",x,"_",shift_ID,".qs"))
                        test <- try(qs::qread(tmp_file),silent = TRUE)
                        if(any(class(test)=="try-error")){
                            tmp_ext <- ext(avg_sdms_tiles, x)
                            terra::window(avg_sdms) <- tmp_ext
                            tmp_data <- spatial_grad(avg_sdms)
                            terra::window(avg_sdms) <- NULL
                            # plug in "real" icells
                            real_cells <- terra::cells(avg_sdms, tmp_ext)
                            tmp_data$icell <- real_cells
                            qs::qsave(tmp_data, tmp_file)
                        }
                    }, mc.cores = ncores)
                    
                    spgrad <- lapply(terra::cells(avg_sdms_tiles), function(x) {
                        tmp_file <- here::here(tmp_dir, paste0("spgrad_tile",x,"_",shift_ID,".qs"))
                        qs::qread(tmp_file)
                    })
                    
                    spgrad <- data.frame(data.table::rbindlist(spgrad))
                    
                    # delete temporary files
                    del_files <- lapply(terra::cells(avg_sdms_tiles), function(x) {
                        tmp_file <- here::here(tmp_dir, paste0("spgrad_tile",x,"_",shift_ID,".qs"))
                        unlink(tmp_file)
                    })
                    
                    qs::qsave(spgrad, spgrad_file)
                    
                } else {
                    spgrad <- spatial_grad(avg_sdms)
                    qs::qsave(spgrad, spgrad_file)
                }
                
            } else {
                spgrad <- spatial_grad(avg_sdms)
                qs::qsave(spgrad, spgrad_file)
            }
            
        }
        
        ## calculate gradient-based climate velocity:
        cat("Calculate Bioclimatic velocities\n")
        
        bVel <- gVelocity(grad = spgrad, slope = ttrend, truncate = TRUE)
        
        ## Across latitude
        cat("Velocity across latitude\n")
        bVelLat <- bVel$Vel * cos(deg_to_rad(bVel$Ang))
        
        # change sign of bVelLat if in the south hemisphere to reflect a velocity away of the tropics
        Cells <- na.omit(terra::as.data.frame(bVelLat$Vel, xy = TRUE, cell = TRUE) )
        SouthCells <- Cells %>% filter(y<0)
        SouthCells <- SouthCells$cell
        if(length(SouthCells)>0){
            bVelLat[SouthCells] <- bVelLat[SouthCells] * -1
        }
        
        #######
        # summary shift velocity over study area
        bv.trend.mean <- as.numeric(global(ttrend,mean,na.rm = TRUE)[,1])
        bv.trend.sd <- as.numeric(global(ttrend,sd,na.rm = TRUE)[,1])
        
        bv.mean <- terra::global(bVel$Vel, mean, na.rm=TRUE)[,1]
        bv.median <- terra::global(bVel$Vel, median, na.rm=TRUE)[,1]
        bv.sd <- terra::global(bVel$Vel, sd, na.rm=TRUE)[,1]
        
        bv.lat.mean <- terra::global(bVelLat$Vel, mean, na.rm=TRUE)[,1]
        bv.lat.median <- terra::global(bVelLat$Vel, median, na.rm=TRUE)[,1]
        bv.lat.sd <- terra::global(bVelLat$Vel, sd, na.rm=TRUE)[,1]
        
        bVelSA_i <- data.frame(bv.trend.mean, bv.trend.sd, 
                               bv.mean, bv.median, bv.sd, 
                               bv.lat.mean, bv.lat.median, bv.lat.sd)
        
        shifts_SA_ens[[i]] <- bVelSA_i
        
        
        #######
        cat("Save bioclimatic velocity maps\n")
        # save shift velocity maps
        terra::writeRaster(bVel, 
                           here::here(shift_dir, paste0("bVel_",shift_ID,".tif")),
                           overwrite=TRUE)
        
        terra::writeRaster(bVelLat, 
                           here::here(shift_dir, paste0("bVelLat_",shift_ID,".tif")),
                           overwrite=TRUE)
        
        #######
        # delete temporary files
        unlink(ttrend_file)
        unlink(spgrad_file)
        
        
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
# # Calc SA range shifts single algorithms
# sdm_algo <- c("GLM","GAM","GBM","MAXNET")
# 
# # load SDM single algorithms projections for each shift in bioshifts database
# sdms_SA <- lapply(1:nrow(shift_info), function(i) {
#     
#     tmp <- sdms[grep(shift_info$ID[i],sdms)]
#     
#     tmp_tif = list.files(tmp,full.names = T,pattern = ".tif")
#     tmp_tif <- tmp_tif[-grep("ClampingMask",tmp_tif)]
#     
#     tmp <- lapply(tmp_tif, function(x) {
#         tmp_years <- strsplit(x, "/")[[1]]
#         tmp_years <- tmp_years[length(tmp_years)]
#         tmp_years <- strsplit(tmp_years, " ")[[1]]
#         tmp_years <- tmp_years[3]
#         
#         tmp_tif_i <- rast(x)
#         tmp_tif_i <- lapply(sdm_algo, function(i){
#             mean(tmp_tif_i[[grep(i,names(tmp_tif_i))]])
#         })
#         tmp_tif_i <- rast(tmp_tif_i)
#         names(tmp_tif_i) <- paste(sdm_algo, tmp_years)
#         return(tmp_tif_i)
#     })
#     
#     return(rast(tmp))
# })
# names(sdms_SA) <- shift_info$ID
# 
# shifts_SA <- list()
# 
# for(i in 1:length(sdms_SA)){ # loop across shifts
#     
#     test <- try({
#         # shift info
#         info_i <- shift_info[i,]
#         
#         # sdms
#         sdms_i <- sdms_SA[[i]]
#         
#         shifts_SA_algo <- list()
#         # run for each algorithm
#         for(j in 1:length(sdm_algo)){ # loop across shifts
#             
#             algo_j <- sdm_algo[j]
#             
#             sdms_i_j <- sdms_i[[grep(algo_j,names(sdms_i))]]
#             
#             # project to equal-area
#             sdms_i_j <- terra::project(sdms_i_j,Eckt)
#             
#             # range 0-1
#             sdms_i_j <- rast(lapply(sdms_i_j, range01raster))
#             
#             # quantiles
#             quants <- c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99)
#             
#             names(sdms_i_j) <- paste0("t",1:nlyr(sdms_i_j))
#             times <- round(shift_info$START[1],0):round(shift_info$END[1],0)
#             
#             # remove if all cells are unsuitable
#             rem <- data.frame(minmax(sdms_i_j))
#             rem <- apply(rem,2,function(x) all(is.na(x)))
#             if(any(rem)){
#                 sdms_i_j <- sdms_i_j[[!rem]]
#                 times <- times[!rem]
#             }
#             
#             tmp <- bioshifts(
#                 x = sdms_i_j,
#                 times = times,
#                 quants = quants,
#                 cores = 2,
#                 metrics = c("centroid","nsCentroid","nsQuants")
#             )
#             tmp <- tmp[,-2]
#             names(tmp)[1:3] <- c("START","END","DUR")
#             # add info
#             tmp$ID <- info_i$ID
#             tmp$Type <- info_i$Type
#             tmp$Species <- info_i$Species
#             tmp$time_period = info_i$time_period
#             tmp$model = "ENSEMBLE"
#             tmp$algo = algo_j
#             
#             shifts_SA_algo[[j]] <- tmp
#         }
#         
#     },silent = TRUE)
#     if(class(test)=="try-error"){
#         shifts_SA_algo[[j]] <- NULL
#     }
#     
#     
#     shifts_SA[[i]] <- rbindlist(shifts_SA_algo)
#     
# }
# 
# shifts_SA <- rbindlist(shifts_SA)
# 
# if(nrow(shifts_SA)>0){
#     write.csv(shifts_SA, 
#               here::here(shift_dir,
#                          paste(sptogo,"shift algo SA.csv")),
#               row.names = FALSE)
#     
# }