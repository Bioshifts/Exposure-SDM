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
# source functions
source("R/my_functions.R")
source("R/bioshiftsFunction.R")
source("R/velocity_functions.R")
# source settings
source("R/settings.R")

# sps I have shift
sp_shift <- rbind(data.frame(sps = list.files(shift_dir("Ter"), pattern = ".csv"),
                             dir = shift_dir("Ter"),
                             realm = "Ter"),
                  data.frame(sps = list.files(shift_dir("Mar"), pattern = ".csv"),
                             dir = shift_dir("Mar"),
                             realm = "Mar"))
length(unique(sp_shift$ID))

sp_shift_data <- lapply(1:nrow(sp_shift), function(i){
    tmp <- read.csv(here::here(sp_shift$dir[i],sp_shift$sps[i]))
    tmp$realm <- sp_shift$realm[i]
    return(tmp)
})
sp_shift_data <- rbindlist(sp_shift_data,fill = TRUE)
length(unique(sp_shift_data$ID))

# filter the ones without edge shift
sp_shift_data <- sp_shift_data[is.na(sp_shift_data$lat.01),]
length(unique(sp_shift_data$ID))

all_sps <- unique(sp_shift_data$Species)
length(all_sps)

for(j in 1:length(all_sps)){ 
    
    sptogo <- all_sps[j]
    # sptogo <- "Acrocephalus_arundinaceus"
    
    sp_shift_data_j <- sp_shift_data %>% filter(Species==sptogo)
    realtogo <- sp_shift_data_j$realm[1]
    
    tmp <- try({
        
        shifts_SA_ens <- data.frame()
        
        for(i in 1:nrow(sp_shift_data_j)){ cat("\r sp",j,"from",length(all_sps),"- shift", i, "from", nrow(sp_shift_data_j))
            
            sp_shift_data_i <- sp_shift_data_j[i, ]
            
            # load shifts
            bVel <- paste(sp_shift_data_i$Species, sp_shift_data_i$ID, sp_shift_data_i$time_period, sp_shift_data_i$Type, realtogo, "SA_bVel.tif", sep = "_")
            bVel <- rast(here::here(shift_dir(realtogo), "raster_files", bVel))
            
            bVelLat <- paste(sp_shift_data_i$Species, sp_shift_data_i$ID, sp_shift_data_i$time_period, sp_shift_data_i$Type, realtogo, "SA_bVelLat.tif", sep = "_")
            bVelLat <- rast(here::here(shift_dir(realtogo), "raster_files", bVelLat))
            
            # load sdm T1
            sdms_i <- here::here(sdm_dir(realtogo),sptogo,gsub("_",".",sptogo),
                                 paste(paste0("proj_",sptogo), sp_shift_data_i$ID, round(sp_shift_data_i$START,0), "SA ens"),
                                 paste(
                                     paste(paste0("proj_",sptogo), sp_shift_data_i$ID, round(sp_shift_data_i$START,0), "SA ens"),
                                     gsub("_",".",sptogo), "ensemble.tif",sep = "_"))
            sdms_i <- rast(sdms_i)
            sdms_i <- mean(sdms_i)
            
            sdms_i <- project(sdms_i, bVel)
            # plot(bVel);dev.off()
            # plot(sdms_i);dev.off()
            
            ########################
            # get edge shift
            edges <- terra::as.data.frame(sdms_i,xy=TRUE,na.rm=TRUE)
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
                window(tmp) <- x
                shift_i <- global(tmp,mean,na.rm=TRUE)
                window(tmp) <- NULL
                return(shift_i)
            })
            edge_shift <- data.frame(t(unlist(edge_shift)))
            names(edge_shift) <- paste0("bv.",c("01","05","1","25","75","9","95","99"))
            
            edge_shift_lat <- sapply(edges_ext, function(x){
                tmp <- bVelLat$Vel
                window(tmp) <- x
                shift_i <- global(tmp,mean,na.rm=TRUE)
                window(tmp) <- NULL
                return(shift_i)
            })
            edge_shift_lat <- data.frame(t(unlist(edge_shift_lat)))
            names(edge_shift_lat) <- paste0("bv.lat.",c("01","05","1","25","75","9","95","99"))
            
            edge_lat <- data.frame(t(edges))
            names(edge_lat) <- paste("lat",c("01","05","1","25","75","9","95","99"),sep=".")
            
            edge_lat <- cbind(edge_lat,edge_shift,edge_shift_lat)
            
            # remove the edge columns
            sp_shift_data_i <- sp_shift_data_i[,1:14]
            
            # add edge columns    
            sp_shift_data_i <- data.frame(sp_shift_data_i,edge_lat)
            
            shifts_SA_ens <- rbind(shifts_SA_ens,sp_shift_data_i)
            
        }
        
    }, silent = TRUE)
    
    if(!class(tmp)=="try-error"){
        
        write.csv(shifts_SA_ens, 
              here::here(shift_dir(realtogo),
                         paste(sptogo,"shift ens SA.csv")),
              row.names = FALSE)
        
    }
    
}



