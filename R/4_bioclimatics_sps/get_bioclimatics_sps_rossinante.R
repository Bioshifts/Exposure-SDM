# Rossinante
# docker run -it --rm --cpus="10.0" --memory="100000m"  -v /home/boliveira/Exposure-SDM:/home/boliveira/project -v /media/seagate/boliveira:/home/boliveira/outputs brunnospatial 

# matrics
# singularity container
# srun -N 1 -n 1 -c 64 --time=9-24:00:00 --mem=64G -J bios -p normal-amd singularity run brunnospatial.sif 

########################
# Setup
rm(list=ls())
gc()

list.of.packages <- c("raster","terra","here",
                      "tidyverse","tictoc","qs","parallelly",
                      "foreach","doParallel","pbapply","sp")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)


########################
N_cores = parallelly::availableCores()

cat("N cores = ", N_cores)

########################

computer="matrics"

if(computer=="rossinante"){
    setwd("/home/boliveira/Exposure-SDM")
    work_dir <- getwd()
    scratch_dir <- "/home/boliveira/outputs"
    output_dir <- "/home/boliveira/outputs"
    occ_dir <- here::here(work_dir,"Data/GBIF_data")
    tmp_dir <- here(scratch_dir,"tmp_dir")
}
if(computer=="matrics"){
    setwd("/users/boliveira/Exposure-SDM")
    work_dir <- getwd()
    scratch_dir <- "/scratch/boliveira/Data"
    output_dir <- work_dir
    occ_dir <- here::here(work_dir,"Data/GBIF_data")
    tmp_dir <- here(scratch_dir,"tmp_dir")
}

# load functions
source("R/my_functions.R")
source("R/Get_Bioclimatics.R")
source("R/get_ecoregions.R")
source("R/get_env_data_for_modeling.R")
# source settings
source("R/settings.R")


N_cores = parallelly::availableCores()

if(computer=="rossinante"){
    ########################
    # get list of species missing env data
    all_sps_togo <- read.csv(here(work_dir,"missing_env_data_sps.csv"))
    nrow(all_sps_togo)
    
    sps_I_have_GBIF <- list.files(occ_dir)
    sps_I_have_GBIF <- gsub(".qs","",sps_I_have_GBIF)
    length(sps_I_have_GBIF)
    
    all_sps_togo <- all_sps_togo[gsub(" ","_",all_sps_togo$sps) %in% sps_I_have_GBIF,]
    nrow(all_sps_togo)
    
    # remove species I already have bioclimatics
    sps_I_have_env <- list.files(here(work_dir,"/Data/Env_data"),recursive = TRUE)
    sps_I_have_env <- gsub(".qs","",sps_I_have_env)
    if(length(sps_I_have_env)>0){
        sps_I_have_env <- sapply(sps_I_have_env, function(x){
            tmp <- strsplit(x,"/")[[1]][2]
            tmp <- strsplit(tmp,"_")[[1]]
            tmp <- paste(tmp[1],tmp[2],sep="_")
            return(tmp)
        })
    }
    
    all_sps <- all_sps_togo[!gsub(" ","_",all_sps_togo$sps) %in% sps_I_have_env,]
    nrow(all_sps)
}
if(computer=="matrics"){
    
    sps_I_have_GBIF <- list.files(occ_dir)
    sps_I_have_GBIF <- gsub(".qs","",sps_I_have_GBIF)
    length(sps_I_have_GBIF)
    
    sps_I_have_env <- list.files(here(work_dir,"/Data/Env_data"),recursive = TRUE)
    sps_I_have_env <- gsub(".qs","",sps_I_have_env)
    if(length(sps_I_have_env)>0){
        sps_I_have_env <- sapply(sps_I_have_env, function(x){
            tmp <- strsplit(x,"/")[[1]][2]
            tmp <- strsplit(tmp,"_")[[1]]
            tmp <- paste(tmp[1],tmp[2],sep="_")
            return(tmp)
        })
    }
    length(sps_I_have_env)
    
    all_sps <- sps_I_have_GBIF[which(!sps_I_have_GBIF %in% sps_I_have_env)]
    length(all_sps)
    
    # merge class
    bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
    bioshifts <- bioshifts_sdms_selection(bioshifts)
    
    # bioshifts %>% group_by(class) %>% tally %>% data.frame
    # bioshifts %>% group_by(class) %>% summarise(n = length(unique(sp_name_std))) %>% data.frame
    
    all_sps <- unique(bioshifts[which(bioshifts$sp_name_std %in% all_sps),c("sp_name_std","class","Eco")])
    nrow(all_sps)
    
    names(all_sps) <- c("sps","class","realm")
    table(all_sps$class)
}

# sort occurrences by size
sizes <- sapply(1:nrow(all_sps), function(i){
    sptogo <- all_sps$sps[i]
    sptogo <- gsub(" ","_",sptogo)
    
    sp_occ <- qs::qread(here::here(occ_dir, paste0(sptogo, ".qs")))
    return(nrow(sp_occ))
})
all_sps$N <- sizes
all_sps <- all_sps[order(all_sps$N),]

# priorities
priorities <- rev(c("Mammalia", "Squamata", "Amphibia", "Liliopsida", "Magnoliopsida", "Pinopsida", "Aves"))
for(i in 1:length(priorities)){
    p1 <- priorities[i]
    rowns_p1 <- which(all_sps$class == p1)
    if(length(rowns_p1)>0){
        sub_p1 <- all_sps[rowns_p1,]
        all_sps <- all_sps[-rowns_p1,]
        all_sps <- rbind(sub_p1,all_sps)
    }
}
# pos_in <- which(all_sps$class %in% priorities)
# pos_out <- which(!1:nrow(all_sps) %in% pos_in)
# pos <- c(sample(pos_in,length(pos_in)), pos_out)
# 
# all_sps <- all_sps[pos,]
# head(all_sps)
# tail(all_sps)

for(i in 1:nrow(all_sps)){ 
    
    ########################
    # Start the clock
    start_time <- Sys.time()
    
    sptogo <- all_sps$sps[i]
    sptogo <- gsub(" ","_",sptogo)
    realm <- all_sps$realm[i]
    
    cat("\n-----------------\nGoing sps", i, "from", nrow(all_sps), "-", sptogo, "\n")
    
    if(computer=="matrics"){
        env_data_dir_togo <- here(work_dir,"Data/Env_data",realm)
    }
    if(computer=="rossinante"){
        env_data_dir_togo <- here(scratch_dir,"Data/Env_data",realm)
    }
    
    vars_dir_togo <- vars_dir(realm)
    
    ########################
    # Check if dir exists for saving files exists
    if(!dir.exists(env_data_dir_togo)){
        dir.create(env_data_dir_togo, recursive = TRUE)
    }  
    
    mask.ras = 
        if(realm=="Ter") { terra::rast(here::here(vars_dir_togo,paste0("model_raster_ter_",my_res,".tif"))) 
        } else {
            if(realm=="Mar") { terra::rast(here::here(vars_dir_togo,"model_raster_mar.tif")) 
            }
        }
    
    ########################
    # Load occurrences for species i
    cat("\nLoading occurrence records\n")
    
    sp_occ <- qs::qread(here::here(occ_dir, paste0(sptogo, ".qs")))
    sp_occ$pa <- 1
    
    sp_occ <- sp_occ %>%
        mutate(x = decimalLongitude,
               y = decimalLatitude) %>%
        dplyr::select(-c('basisOfRecord','speciesKey','decimalLongitude','decimalLatitude')) 
    
    nrow(sp_occ)
    
    if(nrow(sp_occ) > limit_recs){
        sp_occ <- sp_occ[sample(1:nrow(sp_occ), limit_recs, replace = FALSE),]
    }
    
    ########################
    # Load ecoregions
    output_dir <- here::here(tmp_dir, sptogo)
    if(!dir.exists(output_dir)){
        dir.create(output_dir,recursive = TRUE)
    }
    
    cat("\nGetting ecoregions\n")
    
    BA <- get_ecoregions(realm = realm, 
                         sptogo = sptogo,
                         PresAbs = data.frame(sp_occ[,c("x","y")]), 
                         varsdir = vars_dir_togo, 
                         mask.ras = mask.ras,
                         return.shp = TRUE,
                         return.raster = FALSE,
                         check_if_exists = TRUE, 
                         output_dir = output_dir)
    
    BA <- BA$shape_file
    BA <- rasterize(BA, mask.ras)
    
    ########################
    # Create random pseudo-absences
    cat("\nCreating random pseudo-absences\n")
    
    env_range <- temporal_range_env_data(realm) 
    PA_occ <- create_temporal_pseudo_absences(BA, env_range, sp_occ)
    
    # # Test plot
    # back_coords <- vect(data.frame(PA_occ), geom = c('x','y'))
    # sp_coords <- vect(data.frame(sp_occ), geom = c('x','y'))
    # plot(BA)
    # plot(back_coords,add=T,cex=.5)
    # plot(sp_coords,add=T,col = "red",cex=.2)
    # dev.off()
    
    # Add PA to occ
    sp_occ <- rbind(sp_occ, PA_occ)
    
    cat("\nN Presences and absences\n")
    cat(table(sp_occ$pa))
    
    # head(sp_occ)
    # min(sp_occ$year)
    
    ########################
    # Load bioclimatics from each date and cell
    cat("\nCalculating bioclimatics from each date and cell")
    
    # 1) Get layers
    if(realm == "Ter"){
        all_layers <- list.files(here::here(vars_dir_togo,my_res), pattern = ".tif", full.names = TRUE, recursive = TRUE)
        all_layers_names <- list.files(here::here(vars_dir_togo,my_res), pattern = ".tif", recursive = TRUE)
        all_layers_names <- gsub('.tif','',all_layers_names)
        all_layers_names <- sapply(all_layers_names, function(x) strsplit(x,"/")[[1]][2])
    }
    if(realm == "Mar"){
        all_layers <- list.files(here::here(vars_dir_togo,"SST"), pattern = ".tif", full.names = TRUE, recursive = TRUE)
        all_layers_names <- list.files(here::here(vars_dir_togo,"SST"), pattern = ".tif", recursive = TRUE)
        all_layers_names <- gsub('.tif','',all_layers_names)
    }
    
    # 1) Get all possible dates
    sp_occ$date <- paste("01",sp_occ$month, sp_occ$year,sep = "_")
    sp_occ$date <- as.Date(sp_occ$date,"%d_%m_%Y")
    sp_occ$date <- format(sp_occ$date,"%m_%Y")
    
    possibledates <- unique(sp_occ$date)
    
    dim(sp_occ)
    length(possibledates)
    
    # 3) Create ID for merge bioclimatics with occ
    sp_occ$ID <- paste(sp_occ$cell,sp_occ$date,sep = "_")
    
    # 4) Calculate bioclimatics
    # biosclim <- list()
    # for(j in 1:length(possibledates)){ print(cat("\rDate",j,"from",length(possibledates)))
    #     date_j <- possibledates[j]
    #     bios_j <- bioclimatics_from_date(all_layers,
    #                                      date_i = date_j,
    #                                      realm,
    #                                      sp_occ,
    #                                      n_yr_bioclimatic,
    #                                      n_cores=ifelse(realm=="Ter", N_cores, NULL))
    #     biosclim[[j]] <- bios_j
    # }
    biosclim <- parallel::mclapply(1:length(possibledates), function(j){
        date_j <- possibledates[j]
        bios_j <- bioclimatics_from_date(all_layers,
                                         date_i = date_j,
                                         realm,
                                         sp_occ,
                                         n_yr_bioclimatic)
        return(bios_j)
    }, mc.cores = ifelse(realm=="Ter", N_cores, 1))
    
    classes <- sapply(biosclim, class)
    if(any(!classes == "data.frame")){
        keep <- which(classes == "data.frame")
        biosclim <- biosclim[keep]
    }
    biosclim <- data.table::rbindlist(biosclim)
    
    ########################
    # Stop the clock
    end_time <- Sys.time()
    end_time - start_time
    
    cat("\nTook:",round(end_time - start_time,1),"minutes")
    
    # Compile
    biosclim <- filter(biosclim, !duplicated(ID))
    sp_occ <- merge(sp_occ, biosclim, by = "ID", all.x = TRUE)
    sp_occ <- sp_occ[,-"ID"]
    
    table(sp_occ$pa)
    
    # save
    qs::qsave(sp_occ, here::here(env_data_dir_togo,paste0(sptogo,"_",realm,".qs")))
    
    gc()
    
}
