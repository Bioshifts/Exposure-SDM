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

# sptogo="Citharichthys_arctifrons"
# realm <- "Mar"
# IDtogo = "A196_P4"

# sptogo="Veronica_montana"
# realm <- "Ter"
# IDtogo = "A19_P1"

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
# create dirs
tmp_dir_sps <- here::here(tmp_dir,paste(sptogo,IDtogo,realm))
if(!dir.exists(tmp_dir_sps)){
    dir.create(tmp_dir_sps,recursive = TRUE)
}

########################
# Load files

# bioshift data
bioshift_info <- read.csv(here::here(sdm_dir,sptogo,paste(sptogo,"bioshift.csv")))
# only LAT
bioshift_info <- bioshift_info[which(bioshift_info$Type=="LAT"),]

# shift info
shift_info <- read.csv(here::here(sdm_dir,sptogo,paste(sptogo,"shift_info.csv")))
# select the shift of interest
shift_info <- shift_info %>% filter(ID == IDtogo)

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
    all(sapply(round(shift_info$Start[i],0):round(shift_info$End[i],0), function(x){
        any(grepl(x,tmp))
    }))
})

# filter only shifts with projections for all years
shift_info <- shift_info[test_sdms_names,]

# shift type
shift_info <- unique(merge(shift_info,bioshift_info[,c("Type","ID")],by="ID"))
time_perid <- round(min(shift_info$Start),0):round(max(shift_info$End),0)


sdms <- unlist(lapply(IDtogo, function(x) {
    sdms[grep(x,sdms)]
}))
sdms_names <- unlist(lapply(IDtogo, function(x) {
    sdms_names[grep(x,sdms_names)]
}))
sdms_ens <- unlist(lapply(IDtogo, function(x) {
    sdms_ens[grep(x,sdms_ens)]
}))
sdms_ens_names <- unlist(lapply(IDtogo, function(x) {
    sdms_ens_names[grep(x,sdms_ens_names)]
}))


#####################################
# load SDM ensemble projections for the shift of interest

tmp <- sdms_ens[grep(IDtogo,sdms_ens)]
tmp <- tmp[grep(paste(time_perid,collapse = "|"),tmp)]
tmp_tif = list.files(tmp, full.names = T, pattern = "ensemble.tif")
tmp <- lapply(tmp_tif, function(x) mean(rast(x)))
names(tmp) <- time_perid
sdms_ens_SA <- rast(tmp)

# plot(sdms_ens_SA[[1]]);dev.off()

#####################################
# Calc SA range shifts ensemble

# shift ID for saving objects
shift_ID <- paste(shift_info$Species, shift_info$ID, shift_info$time_period, shift_info$Type, realm, "SA", sep = "_")

# sdms
sdms_i <- sdms_ens_SA


for(i in 1:nrow(shift_info)){
    
    shift_ID_togo <- shift_ID[i]
    periods_i <- as.character(round(min(shift_info$Start[i]),0):round(max(shift_info$End[i]),0))
    sdms_i_togo <- sdms_i[[periods_i]]
    
    #######
    # calculate the trend (C/year)
    cat("calculate the trend\n")
    
    ttrend_file <- here::here(tmp_dir_sps, paste0("trend_",shift_ID_togo,".tif"))
    ttrend <- try(terra::rast(ttrend_file),silent = TRUE)
    if(class(ttrend)=="try-error"){
        ttrend = temp_grad(
            sdms_i_togo,
            th = 0.25*nlyr(sdms_i_togo), ## set minimum N obs. to 1/4 time series length
            tempfile = ttrend_file,
            overwrite = TRUE,
            ncores = NULL) # no improvement in speed for trend function
    }
    
    #######
    # Get averaged sdms
    avg_sdm_file <- here::here(tmp_dir_sps, paste0("AvgSDM_",shift_ID_togo,".tif"))
    avg_sdms <- try(terra::rast(avg_sdm_file),silent = TRUE)
    if(class(avg_sdms)=="try-error"){
        avg_sdms <- terra::app(
            sdms_i_togo, 
            fun = mean, na.rm = TRUE, 
            filename=avg_sdm_file,
            overwrite=TRUE,
            cores = ncores)
    }
    
    
    #######
    # Get the spatial gradient (C/km)
    cat("calculate the spatial gradient\n")
    
    spgrad_file <- here::here(tmp_dir_sps, paste0("spgrad_",shift_ID_togo,".tif"))
    spgrad <- try(terra::rast(spgrad_file),silent = TRUE)
    
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
                # avg_sdms_tiles[] <- 1:ncell(avg_sdms_tiles)
                # plot(avg_sdms_tiles);dev.off()
                avg_sdms_tiles <- terra::makeTiles(
                    avg_sdms, 
                    avg_sdms_tiles, 
                    filename = here::here(tmp_dir_sps,paste(shift_ID_togo,"tile_.tif",sep="_")),
                    na.rm = TRUE, overwrite=TRUE)
                
                # x=avg_sdms_tiles[1]
                test <- parallel::mclapply(avg_sdms_tiles, function(x) {
                    tmp_file <- gsub("tile","spgrad",x)
                    test <- try(rast(tmp_file),silent = TRUE)
                    if(any(class(test)=="try-error")){
                        tmp_tile <- terra::rast(x)
                        tmp_data <- spatial_grad(tmp_tile)
                        terra::writeRaster(tmp_data,tmp_file,overwrite = TRUE)
                    }
                }, mc.cores = ncores)
                
                test <- sapply(test,class)
                remove_these <- which(test=="try-error")
                if(length(remove_these)>0){
                    tiles_good <- avg_sdms_tiles[-remove_these]
                } else {
                    tiles_good <- avg_sdms_tiles
                }
                
                spgrad <- terra::vrt(gsub("tile","spgrad",tiles_good),
                                     set_names = TRUE,
                                     overwrite = TRUE)
                
                terra::writeRaster(spgrad,
                                   spgrad_file, 
                                   overwrite = TRUE)
                
                spgrad <- terra::rast(spgrad_file)
                
                # delete temporary files
                unlink(avg_sdms_tiles)
                unlink(gsub("tile","spgrad",avg_sdms_tiles))
                
            } else {
                spgrad = spatial_grad(avg_sdms)
                
                terra::writeRaster(spgrad, spgrad_file)
                
                spgrad <- terra::rast(spgrad_file)
            } 
        } else {
            spgrad = spatial_grad(avg_sdms)
            
            terra::writeRaster(spgrad, spgrad_file)
            
            spgrad <- terra::rast(spgrad_file)
        }
        
    }
    
    if(!terra::ext(spgrad) == terra::ext(ttrend)){
        spgrad <- terra::project(spgrad, ttrend)
    }
    
    ## calculate gradient-based climate velocity:
    cat("Calculate Bioclimatic velocities\n")
    
    spgrad <- terra::project(spgrad, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    ttrend <- terra::project(ttrend, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    
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
    ## Project to equal area for more accurate statistics
    bVel <- terra::project(bVel, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    bVelLat <- terra::project(bVelLat, bVel, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    ttrend <- terra::project(ttrend, bVel, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    sdms_year_1 <- terra::project(sdms_i[[1]], bVel, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
    
    ############
    # get edge shift
    edges <- terra::as.data.frame(sdms_year_1,xy=TRUE,na.rm=TRUE)
    edges <- wtd.quantile(edges[,2], weights = edges[,3],
                          probs = c(0.01, 0.05, 0.1, 0.25, 0.75, 0.90, 0.95, 0.99))
    
    edges_ext_01 <- edges_ext_05 <- edges_ext_1 <- edges_ext_25 <- edges_ext_75 <- edges_ext_9 <- edges_ext_95 <- edges_ext_99 <- ext(sdms_year_1)
    
    edges_ext_01[4] <- edges[1]
    edges_ext_05[4] <- edges[2]
    edges_ext_1[4] <- edges[3]
    edges_ext_25[4] <- edges[4]
    edges_ext_75[3] <- edges[5]
    edges_ext_9[3] <- edges[6]
    edges_ext_95[3] <- edges[7]
    edges_ext_99[3] <- edges[8]
    
    edges_ext <- list(edges_ext_01,edges_ext_05,edges_ext_1,edges_ext_25,edges_ext_75,edges_ext_9,edges_ext_95,edges_ext_99)
    x=edges_ext[[1]]
    edge_shift <- sapply(edges_ext, function(x){
        tmp <- bVel$Vel
        tmp <- crop(tmp, x)
        tmp2 <- sdms_year_1
        tmp2 <- crop(tmp2, x)
        shift_i <- as.numeric(terra::global(tmp,"mean",weights=tmp2,na.rm=TRUE))
        return(shift_i)
    })
    edge_shift <- data.frame(t(unlist(edge_shift)))
    names(edge_shift) <- paste0("bv.",c("01","05","1","25","75","9","95","99"))
    
    edge_shift_lat <- sapply(edges_ext, function(x){
        tmp <- bVelLat$Vel
        tmp <- crop(tmp, x)
        tmp2 <- sdms_year_1
        tmp2 <- crop(tmp2, x)
        shift_i <- as.numeric(terra::global(tmp,"mean",weights=tmp2,na.rm=TRUE))
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
    
    shifts_SA_ens <- data.frame(shift_info[i,],
                                bv.trend.mean, bv.trend.sd, 
                                bv.mean, bv.median, bv.sd, 
                                bv.lat.mean, bv.lat.median, bv.lat.sd,
                                edge_lat)
    
    if(nrow(shifts_SA_ens)>0){
        write.csv(shifts_SA_ens, 
                  here::here(shift_dir,
                             paste0(shift_ID_togo,".csv")),
                  row.names = FALSE)
        
    }
    
    #######
    cat("Save bioclimatic velocity maps\n")
    # save shift velocity maps
    terra::writeRaster(bVel, 
                       here::here(shift_dir, "raster_files", paste0(shift_ID_togo,"_bVel.tif")),
                       overwrite=TRUE)
    
    terra::writeRaster(bVelLat, 
                       here::here(shift_dir, "raster_files", paste0(shift_ID_togo,"_bVelLat.tif")),
                       overwrite=TRUE)
    
}

#######
# delete temporary files
unlink(tmp_dir_sps,recursive = TRUE)
