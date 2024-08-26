
rm(list=ls())
gc()

list.of.packages <- c("dplyr","terra", "tidyr","Hmisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)


###
# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}


source("R/settings.R")

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
print(command_args)

spstogo <- as.character(paste(command_args[1], collapse = " "))
SAtogo <- as.character(paste(command_args[2], collapse = " "))
ECO <- as.character(paste(command_args[3], collapse = " "))
Start <- as.character(paste(command_args[4], collapse = " "))
End <- as.character(paste(command_args[5], collapse = " "))

# spstogo <- "Abies_alba"
# SAtogo <- "A156_P1"
# ECO <- "Ter"
# Start <- "1950"
# End <- "2005"

cat("\rrunning species", spstogo)

# load connectivity raster
norm_curr <- try(rast(here::here(connectivity_data_dir(ECO),"norm_curr.tif")), silent = TRUE)
norm_curr_rc <- try(rast(here::here(connectivity_data_dir(ECO),"norm_curr_rc.tif")), silent = TRUE)
# plot(norm_curr);dev.off()
# plot(norm_curr_rc);dev.off()

# create dir to store connectivity SA
conn_SA_dir <- here(scratch_dir,"Data/Connectivity/Connectivity_SA")
if(!dir.exists(conn_SA_dir)){
    dir.create(conn_SA_dir)
}

# load shift
time_period <- paste(Start,End,sep="-")
shift_ID <- paste(spstogo, SAtogo, time_period, ECO, "SA_edges.csv", sep = "_")
tmp <- here(shift_dir(ECO),shift_ID)
shifts_ens <- data.frame(read.csv(tmp))[1,]

shift_ID <- paste(spstogo, SAtogo, time_period, ECO, "Conn", sep = "_")

# create dir to store results
output_dir <- here(work_dir,paste0("Data/Connectivity_SA_edges/",ECO))
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}

# load SA polygon
SA <- terra::vect(here::here(SA_shps_dir,paste0(SAtogo,".shp")))

# Get connectivity at the study area
norm_curr_SA_file <- here(conn_SA_dir,paste0(SAtogo,"_norm_curr.tif"))
if(!file.exists(norm_curr_SA_file)){
    terra::window(norm_curr) <- ext(SA)
    norm_curr_SA <- terra::mask(norm_curr, SA, 
                                filename = norm_curr_SA_file, 
                                overwrite = TRUE)
    terra::window(norm_curr) <- NULL
} else {
    norm_curr_SA <- terra::rast(norm_curr_SA_file)
}
# plot(norm_curr);dev.off()

norm_curr_rc_SA_file <- here(conn_SA_dir,paste0(SAtogo,"_norm_curr_rc.tif"))
if(!file.exists(norm_curr_rc_SA_file)){
    terra::window(norm_curr_rc) <- ext(SA)
    norm_curr_rc_SA <- terra::mask(norm_curr_rc, SA, 
                                   filename = norm_curr_rc_SA_file, overwrite = TRUE)
    terra::window(norm_curr_rc) <- NULL
} else {
    norm_curr_rc_SA <- terra::rast(norm_curr_rc_SA_file)
}
# plot(norm_curr_rc);dev.off()

# load sdm i
proj_name <- paste(paste0("proj_",spstogo),SAtogo, Start, "SA ens")
sdm <- try(rast(here::here(sdm_dir(ECO),spstogo,gsub("_",".",spstogo),proj_name,
                             paste(proj_name,gsub("_",".",spstogo),"ensemble.tif",sep = "_"))),
             silent = TRUE)

if(!class(SA)=="try-error" & !class(sdm)=="try-error"){
    
    filetosave <- here(output_dir, paste0(shift_ID,".csv"))
    
    if(!file.exists(filetosave)){
        
        # avg sdm
        sdm <- mean(sdm)
        # proj sdm
        sdm <- terra::project(sdm, norm_curr_SA)
        # plot(sdm);dev.off()
        # plot(SA);dev.off()
        
        # weighted connectivity
        w_norm_curr_SA <- norm_curr_SA*sdm
        # plot(w_norm_curr_SA);dev.off()
        
        ## Project to equal area for more accurate statistics
        w_norm_curr_SA <- terra::project(w_norm_curr_SA, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
        norm_curr_rc_SA <- terra::project(norm_curr_rc_SA, Eckt, threads=TRUE, use_gdal=TRUE, gdal=TRUE)
        SA <- terra::project(SA, Eckt)
        
        #######
        # Connectivity over the species edges
        edges <- data.frame(lower = range(c(shifts_ens$edge_y_0.25, shifts_ens$edge_y_0.05)),
                            CE = range(c(shifts_ens$edge_y_0.25, shifts_ens$edge_y_0.75)),
                            upper = range(c(shifts_ens$edge_y_0.75, shifts_ens$edge_y_0.95)))
        
        conn_j <- list()
        
        for(j in 1:ncol(edges)){
            
            ext_j <- ext(SA)
            ext_j[3] <- edges[1,j]
            ext_j[4] <- edges[2,j]
            
            # average connectivity
            conn.mean <- as.numeric(terra::global(w_norm_curr_SA, mean, na.rm=TRUE)[,1])
            conn.median <- as.numeric(terra::global(w_norm_curr_SA, median, na.rm=TRUE)[,1])
            conn.sd <- as.numeric(terra::global(w_norm_curr_SA, sd, na.rm=TRUE)[,1])
            
            # classified connectivity
            conn.class <- table(terra::values(round(norm_curr_rc_SA,0)))
            # 1) impeded--no movement, completely blocked by barriers or resistance
            # 2) diffuse--movement is unimpeded (good thing)
            # 3) intensified--movement is partly bottlenecked
            # 4) channelized--movement is completely bottlenecked.
            impeded_Ncell <- conn.class[1]
            diffuse_Ncell <- conn.class[2]
            intensified_Ncell <- conn.class[3]
            channelized_Ncell <- conn.class[4]
            Ncell <- sum(impeded_Ncell,diffuse_Ncell,intensified_Ncell,channelized_Ncell)
            impeded_prop <- impeded_Ncell/Ncell
            diffuse_prop <- diffuse_Ncell/Ncell
            intensified_prop <- intensified_Ncell/Ncell
            channelized_prop <- channelized_Ncell/Ncell
            
            conn <- data.frame(conn.mean,
                               conn.median,
                               conn.sd,
                               impeded_Ncell,
                               diffuse_Ncell,
                               intensified_Ncell,
                               channelized_Ncell,
                               Ncell,
                               impeded_prop,
                               diffuse_prop,
                               intensified_prop,
                               channelized_prop)
            names(conn) <- paste(names(conn),names(edges)[j],sep="_")
            
            conn_j[[j]] <- conn
            
        }
        conn_j <- do.call(cbind,conn_j)
        conn_j <- data.frame(ID = SAtogo,
                             Species = spstogo,
                             time_period = time_period,
                             conn_j)
        
        #######
        # save connectivity results
        write.csv(conn_j, filetosave, row.names = FALSE)
    }
}
