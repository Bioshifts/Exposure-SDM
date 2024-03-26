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
IDtogo <- as.character(paste(command_args[3], collapse = " "))

# sptogo="Fucus_vesiculosus"
# realm <- "Mar"
# IDtogo = "A176_P1"

print(sptogo)
print(realm)
print(IDtogo)


########################
# source functions
source("R/my_functions.R")
source("R/bioshiftsFunction.R")
source("R/velocity_functions.R")
# source settings
source("R/settings.R")

check_if_file_exists <- TRUE

shift_dir <- shift_dir(realm)
shiftplot_dir <- shiftplot_dir(realm)
sdm_dir <- sdm_dir(realm)

# create dirs
if(!dir.exists(shift_dir)){
    dir.create(shift_dir,recursive = TRUE)
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

# select the shift of interest
shift_info <- shift_info %>% filter(ID == IDtogo)

#####################################
# load SDM ensemble projections for the shift of interest

tmp <- sdms_ens[grep(shift_info$ID,sdms_ens)]
tmp_tif = list.files(tmp,full.names = T,pattern = ".tif")
tmp <- lapply(tmp_tif, function(x) mean(rast(x)))
names(tmp) <- round(shift_info$START,0):round(shift_info$END,0)
sdms_ens_SA <- rast(tmp)

# plot(sdms_ens_SA[[1]]);dev.off()

#####################################
# Calc SA range shifts ensemble

# shift ID for saving objects
shift_ID <- paste(shift_info$Species, shift_info$ID, shift_info$time_period, shift_info$Type, realm, "SA",sep = "_")

# sdms
sdms_i <- sdms_ens_SA

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
        SA_i <- terra::vect(here::here(SA_shps_dir,paste0(shift_info$ID,".shp")))
        
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

############
# get edge shift
edges <- terra::as.data.frame(sdms_i[[1]],xy=TRUE,na.rm=TRUE)
edges <- wtd.quantile(edges[,2], weights = edges[,3],
                      probs = c(0.01, 0.05, 0.1, 0.25, 0.75, 0.90, 0.95, 0.99))

edges_ext_01 <- edges_ext_05 <- edges_ext_1 <- edges_ext_25 <- edges_ext_75 <- edges_ext_9 <- edges_ext_95 <- edges_ext_99 <- ext(sdms_i)

edges_ext_01[4] <- edges[1]
edges_ext_05[4] <- edges[2]
edges_ext_1[4] <- edges[3]
edges_ext_25[4] <- edges[4]
edges_ext_75[3] <- edges[5]
edges_ext_9[3] <- edges[6]
edges_ext_95[3] <- edges[7]
edges_ext_99[3] <- edges[8]

edges_ext <- list(edges_ext_01,edges_ext_05,edges_ext_1,edges_ext_25,edges_ext_75,edges_ext_9,edges_ext_95,edges_ext_99)

edge_shift <- sapply(edges_ext, function(x){
    tmp <- bVel$Vel
    test <- try({
        window(tmp) <- x
        shift_i <- global(tmp,mean,na.rm=TRUE)
        window(tmp) <- NULL
    },silent = TRUE)
    if(class(test)=="try-error"){
        tmp <- crop(tmp, x)
        shift_i <- global(tmp,mean,na.rm=TRUE)
    }
    return(shift_i)
})
edge_shift <- data.frame(t(unlist(edge_shift)))
names(edge_shift) <- paste0("bv.",c("01","05","1","25","75","9","95","99"))

edge_shift_lat <- sapply(edges_ext, function(x){
    tmp <- bVelLat$Vel
    test <- try({
        window(tmp) <- x
        shift_i <- global(tmp,mean,na.rm=TRUE)
        window(tmp) <- NULL
    },silent = TRUE)
    if(class(test)=="try-error"){
        tmp <- crop(tmp, x)
        shift_i <- global(tmp,mean,na.rm=TRUE)
    }
    return(shift_i)
})
edge_shift_lat <- data.frame(t(unlist(edge_shift_lat)))
names(edge_shift_lat) <- paste0("bv.lat.",c("01","05","1","25","75","9","95","99"))

edge_lat <- data.frame(t(edges))
names(edge_lat) <- paste("lat",c("01","05","1","25","75","9","95","99"),sep=".")

edge_lat <- cbind(edge_lat,edge_shift,edge_shift_lat)

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

shifts_SA_ens <- data.frame(bv.trend.mean, bv.trend.sd, 
                            bv.mean, bv.median, bv.sd, 
                            bv.lat.mean, bv.lat.median, bv.lat.sd,
                            edge_lat)


#######
cat("Save bioclimatic velocity maps\n")
# save shift velocity maps
terra::writeRaster(bVel, 
                   here::here(shift_dir, "raster_files", paste0(shift_ID,"_bVel.tif")),
                   overwrite=TRUE)

terra::writeRaster(bVelLat, 
                   here::here(shift_dir, "raster_files", paste0(shift_ID,"_bVelLat.tif")),
                   overwrite=TRUE)

#######
# delete temporary files
unlink(ttrend_file)
unlink(spgrad_file)
unlink(avg_sdm_file)


if(nrow(shifts_SA_ens)>0){
    write.csv(shifts_SA_ens, 
              here::here(shift_dir,
                         paste0(shift_ID,".csv")),
              row.names = FALSE)
    
}

