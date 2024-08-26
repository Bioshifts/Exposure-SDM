
# run in singularity container using:
# srun -N 1 -n 1 -c 7 --time=05:00:00 --mem=64G -J vel_edge -p bigmem-amd singularity run brunnospatial.sif 

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

work_dir <- getwd()

# create dir to store results
output_dir <- here(work_dir,"Data/Velocity_SA_edges/Ter")
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}
output_dir <- here(work_dir,"Data/Velocity_SA_edges/Mar")
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}

# get all shifts
shifts_ens_files <- list.files(shift_dir("Mar"), pattern = ".csv", full.names = TRUE)
length(shifts_ens_files)

shifts_ens <- lapply(shifts_ens_files, function(x) {
    tmp <- data.frame(read.csv(x))
    if(any(names(tmp)=="START")){
        pos <- which(names(tmp)=="START")
        names(tmp)[pos] <- "Start"
        pos <- which(names(tmp)=="END")
        names(tmp)[pos] <- "End"
        write.csv(tmp,x,row.names = FALSE)
    }
    if(ncol(tmp)>38){
        write.csv(tmp[,1:38],x,row.names = FALSE)
    }
    return(tmp)
})

shifts_ens <- data.table::rbindlist(shifts_ens)
shifts_ens_mar <- data.frame(shifts_ens)
shifts_ens_mar$ECO = "Marine"

length(unique(shifts_ens_mar$Species))
nrow(shifts_ens_mar)

shifts_ens_files <- list.files(shift_dir("Ter"), pattern = ".csv", full.names = TRUE)
length(shifts_ens_files)

shifts_ens <- lapply(shifts_ens_files, function(x){
    tmp <- data.frame(read.csv(x))
    if(any(names(tmp)=="START")){
        pos <- which(names(tmp)=="START")
        names(tmp)[pos] <- "Start"
        pos <- which(names(tmp)=="END")
        names(tmp)[pos] <- "End"
        write.csv(tmp,x,row.names = FALSE)
    }
    if(ncol(tmp)>38){
        write.csv(tmp[,1:38],x,row.names = FALSE)
    }
    return(tmp)
})

shifts_ens <- data.table::rbindlist(shifts_ens)
shifts_ens_ter <- data.frame(shifts_ens)
shifts_ens_ter$ECO = "Terrestrial"

length(unique(shifts_ens_ter$Species))
nrow(shifts_ens_ter)

shifts_ens <- rbind(shifts_ens_mar, shifts_ens_ter)

length(unique(shifts_ens$Species))
length(unique(shifts_ens$Species[which(shifts_ens$ECO=="Terrestrial")]))
length(unique(shifts_ens$Species[which(shifts_ens$ECO=="Marine")]))


for(i in 1:nrow(shifts_ens)){ cat("\r",i,"from",nrow(shifts_ens))
    
    shifts_ens_i <- shifts_ens[i,]
    
    ECO_i <- shifts_ens_i$ECO
    if(ECO_i=="Marine"){
        var_i <- "sst"
        eco_sdm <- "Mar"
    } else {
        var_i <- "mat"
        eco_sdm <- "Ter"
    }
    name_SA_i <- shifts_ens_i$ID
    sp_i <- shifts_ens_i$Species
    Start_i <- shifts_ens_i$Start
    # load SA i
    SA_i <- try(rast(here::here(velocity_SA_dir,paste(name_SA_i,var_i,"gVelLat.tif",sep="_"))),
                silent = TRUE)
    
    # load sdm i
    proj_name <- paste(paste0("proj_",sp_i),name_SA_i, Start_i, "SA ens")
    sdm_i <- try(rast(here::here(sdm_dir(eco_sdm),sp_i,gsub("_",".",sp_i),proj_name,
                                 paste(proj_name,gsub("_",".",sp_i),"ensemble.tif",sep = "_"))),
                 silent = TRUE)
    
    if(!class(SA_i)=="try-error" & !class(sdm_i)=="try-error"){
        
        filetosave <- here(work_dir,"Data/Velocity_SA_edges", eco_sdm,
                           paste(sp_i,name_SA_i,"vel_SA_edges.csv"))
        
        if(!file.exists(filetosave)){
            
            # avg sdm
            sdm_i <- mean(sdm_i)
            # proj sdm
            sdm_i <- terra::project(sdm_i,SA_i)
            # plot(sdm_i);dev.off()
            # plot(SA_i);dev.off()
            # remove very low suitability
            thrs <- as.numeric(terra::global(sdm_i, quantile, probs = 0.05, na.rm = TRUE))
            sdm_i <- ifel(sdm_i < thrs, NA, sdm_i)
            
            edges <- terra::as.data.frame(sdm_i,xy=TRUE,na.rm=TRUE)
            edges_pos <- wtd.quantile(edges[,2], weights = edges[,3],
                                      probs = c(0.01, 0.05, 0.1, 0.25, 0.75, 0.90, 0.95, 0.99))
            
            # get extents
            ext01 <- ext05 <- ext1 <- ext25 <- ext75 <- ext9 <- ext95 <- ext99 <- terra::ext(SA_i)
            # TE
            ext01[4] <- edges_pos[1]
            ext05[4] <- edges_pos[2]
            ext1[4] <- edges_pos[3]
            ext25[4] <- edges_pos[4]
            # LE
            ext75[3] <- edges_pos[5]
            ext9[3] <- edges_pos[6]
            ext95[3] <- edges_pos[7]
            ext99[3] <- edges_pos[8]
            
            edges <- list(edge01 = ext01,
                          edge05 = ext05,
                          edge1 = ext1,
                          edge25 = ext25,
                          edge5 = ext(SA_i),
                          edge75 = ext75,
                          edge9 = ext9,
                          edge95 = ext95,
                          edge99 = ext99)
            
            vel_edge_j <- parallel::mclapply(edges, function(x){
                sdm_crop <- crop(sdm_i, x)
                SA_crop <- crop(SA_i, x)
                as.numeric(terra::global(SA_crop,"mean",weights=sdm_crop,na.rm=TRUE))
            }, mc.cores = 7)
            
            vel_edge_j <- unlist(vel_edge_j)
            vel_edge_j <- data.frame(t(vel_edge_j))
            names(vel_edge_j) <- paste("vel",names(vel_edge_j),sep = "_")
            vel_edge_j <- cbind(shifts_ens[i,1:6],vel_edge_j)
            
            if(any(is.na(shifts_ens_i[,15:22]))){
                shifts_ens_i[,15:22] <- edges_pos
                
                write.csv(shifts_ens_i[,1:38], 
                          here(shift_dir("Mar"),
                               paste(shifts_ens_i$Species,shifts_ens_i$ID,shifts_ens_i$time_period,"LAT",eco_sdm,"SA.csv",
                                     sep="_")),
                          row.names = FALSE)
            }
            
            write.csv(vel_edge_j, 
                      filetosave,
                      row.names = FALSE)
        }
    }
}

