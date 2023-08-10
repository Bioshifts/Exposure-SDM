# Setup
rm(list=ls())
gc()

# devtools::install_github("sjevelazco/flexsdm@HEAD")

list.of.packages <- c("terra","rnaturalearthdata","biomod2",
                      "tidyverse","tictoc","tidyterra","ggplot2","data.table")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)



########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
sptogo <- as.character(paste(command_args[1], collapse = " "))
realm <- as.character(paste(command_args[2], collapse = " "))

print(sptogo)
print(realm)

# sptogo="Deania_calcea"
# sptogo="Abra_alba"
# realm <- "Mar"
# 
# sptogo="Platanthera_cooperi"
# realm <- "Ter"

check_if_exists = FALSE

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    
    scratch_dir <- "/lustre/oliveirab"
    work_dir <- getwd()
    
    setwd(work_dir)
    
    output_dir <- here::here(scratch_dir,"SDMs",realm,sptogo)
    occ_dir <- here::here(work_dir,"Data/Env_data",realm)
    
    bios_dir <- here::here(work_dir,"Data",realm,"bio_proj")
    if(realm == "Ter") { bios_dir <- gsub(realm,"Land",bios_dir)}
    if(realm == "Mar") { bios_dir <- gsub(realm,"Marine",bios_dir)}
}
# create dirs
if(!dir.exists(output_dir)){
    dir.create(output_dir,recursive = TRUE)
}

########################
# Check if file exists

if(check_if_exists){
    file.test <- here::here(output_dir,gsub("_",".",sptogo))
    
    # check if there a file with ensemble model outputs
    tmp = list.files(here::here(output_dir,
                                realm,
                                sptogo))
    tmp = any(grepl("ens_SDM",tmp))
    
    if(tmp){
        stop("No need to run this. SDMs were fitted for this species!")
    } 
}

########################
# load functions
source("R/my_functions.R")
source("R/get_ecoregions.R")
source("R/correct_colinvar_pca.R")
source("R/plot_SA_location.R")
# source settings
source("R/settings.R")

vars_dir <- get_varsdir(computer,realm)
mask.ras = 
    if(realm=="Ter") { terra::rast(here::here(vars_dir,"model_raster_ter_5km.tif")) 
    } else {
        if(realm=="Mar") { terra::rast(here::here(vars_dir,"model_raster_mar.tif")) 
        }
    }

if(realm=="Ter") { scratch_dir <- here::here(scratch_dir,"cruts") } 
if(realm=="Mar") { scratch_dir <- here::here(scratch_dir,"SST") }

########################
# Load env data for species i
PresAbs <- qs::qread(here::here(occ_dir,paste0(sptogo,"_",realm,".qs")))
PresAbs <- na.omit(PresAbs)

table(PresAbs$pa)

# filter presences
sp_occ <- PresAbs %>% filter(pa == 1)
# filter pseudo-absences
back_occ <- PresAbs %>% filter(pa == 0)


########################
# Load ecoregions
BA <- get_ecoregions(realm = realm, 
                     PresAbs = sp_occ, 
                     varsdir = vars_dir, 
                     mask.ras = mask.ras,
                     return.shp = TRUE)
BA_shp <- BA[[2]]
BA <- BA[[1]]
# test
# plot(BA)
# plot(vect(data.frame(back_occ), geom = c("x","y")),add=T,col = "red")
# plot(vect(data.frame(sp_occ), geom = c("x","y")),add=T)
# dev.off()

########################
# Study info

biov1 <- read.csv(here::here(work_dir,"Data/Bioshifts/biov1_fixednames.csv"),header = T)
biov1 <- biov1 %>% filter(sp_name_std_v1 == sptogo)

# Use LAT ELE shifts
biov1$Type[which(!is.na(biov1$Azimuth))] <- "LAT" # All obs type == HOR contain Azimuth value
biov1 <- biov1[which(biov1$Type=="ELE" | biov1$Type=="LAT"),]

# Use temporal period from the environmental data
biov1 <- biov1 %>% filter(START >= temporal_range_env_data[[realm]][1] + n_yr_bioclimatic)

# save
write.csv(biov1, here::here(output_dir,paste(sptogo,"bioshift.csv")),row.names = FALSE)

# get info relevant for SDM  
shift_info<- biov1[,c("ID","START","END")]
if(any(duplicated(shift_info))){
    shift_info<- shift_info[-which(duplicated(shift_info)),]
}
shift_info$time_period <- paste(round(shift_info$START), round(shift_info$END), sep = "-")
shift_info$Species <- sptogo
#save
write.csv(shift_info, here::here(output_dir,paste(sptogo,"shift_info.csv")),row.names = FALSE)

# get study ID
StudyID <- shift_info$ID

# get start and end periods
Start_period <- round(shift_info$START) # transform 2020.5 in 2020 << Think if makes sense to use 2020.5
End_period <- round(shift_info$END)
Time_periods = unique(sapply(1:nrow(shift_info), function(i) paste(Start_period[i], End_period[i], sep = "-")))

########################
# Load study area
fgdb <- here::here(work_dir,"Data/Bioshifts/Study_Areas_v1/Study_Areas.gdb")

StudyArea <- lapply(StudyID, function(x) {
    tmp <- terra::vect(fgdb,layer=x)
})
StudyArea <- vect(StudyArea)
StudyArea$NAME = StudyID

########################
# Plot for study area location
# create dir if doesn't exist
plot_SA_location(sptogo, 
                 BA,
                 realm, 
                 PresAbs, 
                 StudyArea, 
                 maps.dir = output_dir, 
                 show.plot = FALSE)

########################
# Get bioclimatics for the BG at each time period
bioclimatics_periods <- list()

for(i in 1:length(Time_periods)){
    
    if(file.exists(here::here(output_dir,paste(sptogo,Time_periods[i],"BG start bios.tif")))){
        t1 <- terra::rast(here::here(output_dir,paste(sptogo,Time_periods[i],"BG start bios.tif")))
        t2 <- terra::rast(here::here(output_dir,paste(sptogo,Time_periods[i],"BG end bios.tif")))
        
        tmp <- list(back_start = t1,
                    back_end = t2)
        
        bioclimatics_periods[[i]] <- tmp
    } else {
        p_i <- as.character(Time_periods[i])
        p_i <- strsplit(p_i,"-")[[1]]
        p_i <- as.numeric(p_i)
        
        p_1 = p_i[1]
        p_2 = p_i[2]
        
        # load in
        back_start <- terra::rast(here::here(bios_dir,paste0("bios_",realm,"_",p_1,".tif")))
        back_end <- terra::rast(here::here(bios_dir,paste0("bios_",realm,"_",p_2,".tif")))
        
        # mask
        back_start <- terra::mask(terra::crop(back_start,BA_shp),BA_shp)
        back_end <- terra::mask(terra::crop(back_end,BA_shp),BA_shp)
        
        tmp <- list(back_start = back_start, back_end = back_end)
        
        bioclimatics_periods[[i]] <- tmp
        
        # save
        writeRaster(bioclimatics_periods[[i]]$back_start,
                    here::here(output_dir,paste(sptogo,Time_periods[i],"BG start bios.tif")), overwrite=TRUE)
        writeRaster(bioclimatics_periods[[i]]$back_end,
                    here::here(output_dir,paste(sptogo,Time_periods[i],"BG end bios.tif")), overwrite=TRUE)
        
    }
}
names(bioclimatics_periods) <- Time_periods


# ########################
# Reduce collinearity among the predictors
# Run SDMs of PCA axes

# create an assemble of environmental data coming from species occurrences + background data at each time period

# 1) add background data ()
data_env <- lapply(bioclimatics_periods, function(x){
    t1 = as.data.frame(x$back_start, na.rm=TRUE)
    t2 = as.data.frame(x$back_end, na.rm=TRUE)
    
    if(nrow(t1) > 20000){
        t1 <- t1[sample(1:nrow(t1),20000),]
        t2 <- t1[sample(1:nrow(t1),20000),]
    }
    
    return(rbind(t1,t2))
})

data_env <- data.table::rbindlist(data_env)

# 2) add pres/abs data
data_env <- rbind(data_env, 
                  data.frame(PresAbs)[,names(data_env)])

dim(data_env)

# set a cap for the size of dataset
mycap <- 100000
if(nrow(data_env) > mycap){
    data_env <- data_env[sample(1:nrow(data_env),mycap), ]
}

# get uncorrelated PCs
new_data <- correct_colinvar_pca(env_layer = data_env, 
                                 proj = bioclimatics_periods, 
                                 PresAbs = PresAbs)

# update data
bioclimatics_periods <- new_data$proj

PresAbsFull <- data.table(data.frame(PresAbs[,c("x","y","pa")], new_data$PresAbs))
PresAbsFull <- PresAbsFull[order(PresAbsFull$pa, decreasing = TRUE),]

########################
# Crop to the SA 

bioclimatics_periods_SA <- list()

for(i in 1:nrow(shift_info)){
    if(file.exists(here::here(output_dir,paste(sptogo,paste(shift_info$ID, shift_info$time_period)[i],"SA start bios.tif")))){
        t1 <- terra::rast(here::here(output_dir,paste(sptogo,paste(shift_info$ID, shift_info$time_period)[i],"SA start bios.tif")))
        t2 <- terra::rast(here::here(output_dir,paste(sptogo,paste(shift_info$ID, shift_info$time_period)[i],"SA end bios.tif")))
        
        tmp <- list(back_start = t1,
                    back_end = t2)
        
        bioclimatics_periods_SA[[i]] <- tmp
        
    } else {
        
        tmp1 <- bioclimatics_periods[[shift_info$time_period[i]]]$back_start
        tmp1 <- terra::crop(
            terra::mask(tmp1,
                        StudyArea[StudyArea$NAME==shift_info$ID[i]]),
            StudyArea[StudyArea$NAME==shift_info$ID[i]])
        
        tmp2 <- bioclimatics_periods[[shift_info$time_period[i]]]$back_end
        tmp2 <- terra::crop(
            terra::mask(tmp2,
                        StudyArea[StudyArea$NAME==shift_info$ID[i]]),
            StudyArea[StudyArea$NAME==shift_info$ID[i]])
        
        tmp <- list(back_start = tmp1,
                    back_end = tmp2)
        
        bioclimatics_periods_SA[[i]] <- tmp
        
        # save
        writeRaster(bioclimatics_periods_SA[[i]]$back_start,
                    here::here(output_dir,paste(sptogo,
                                                paste(shift_info$ID, shift_info$time_period)[i],"SA start bios.tif")), 
                    overwrite=TRUE)
        writeRaster(bioclimatics_periods_SA[[i]]$back_end,
                    here::here(output_dir,paste(sptogo,
                                                paste(shift_info$ID, shift_info$time_period)[i],"SA end bios.tif")), 
                    overwrite=TRUE)
    }
}
names(bioclimatics_periods_SA) <- paste(shift_info$ID, shift_info$time_period)


########################

# P.S.: The code bellow works using the version 4.2-3 of biomod2
# installation :
# devtools::install_github("biomodhub/biomod2", dependencies = TRUE)

# Documentation for functions can be found here: 
# https://biomodhub.github.io/biomod2/index.html

# Fit SDMs

# N cpus to use
N_cpus = 12

table(PresAbsFull$pa)

resp <- as.vector(PresAbsFull$pa)
resp[resp==0] <- NA

data_sp <- BIOMOD_FormatingData(
    resp.name = sptogo, 
    resp.var = resp, 
    expl.var = data.frame(PresAbsFull[,-1:-3]), 
    resp.xy = data.frame(PresAbsFull[,1:2]),
    dir.name = here::here(output_dir),
    PA.strategy = "random",
    PA.nb.rep = 1,
    PA.nb.absences = length(which(is.na(resp))))

model_sp <- BIOMOD_Modeling(
    data_sp, 
    models = c("GLM","GBM","GAM","MAXENT.Phillips.2"),
    CV.perc = 70,
    CV.strategy = "kfold",
    CV.k = 3,
    CV.do.full.models = FALSE,
    prevalence = 0.5,
    metric.eval = c("TSS", "ROC"),
    scale.models = TRUE,
    nb.cpu = N_cpus)


ens_model_sp <- BIOMOD_EnsembleModeling(
    model_sp, 
    metric.select = "TSS",
    metric.eval = "TSS")



########################
# Evaluate
myevals = get_evaluations(model_sp)
myevals[which(myevals$metric.eval=="TSS"),]

write.csv(myevals, 
          here::here(output_dir,paste0(sptogo,"_CV.csv")),
          row.names = FALSE)

p1 <- ggplot(myevals, aes(x=validation, y = algo, color = algo))+
    geom_boxplot()+
    facet_wrap(.~metric.eval, scales = "free")

{
    pdf(width = 10, height = 5, file = here::here(output_dir,paste0(sptogo,"_eval.pdf")))
    print(p1)
    dev.off()
    
    }

########################
# save basic data
basicD <- data.frame(Species = sptogo,
                     N_occ = length(which(PresAbsFull$pa==1)),
                     N_background = length(which(PresAbsFull$pa==0)),
                     N_study_areas = length(StudyArea),
                     N_time_periods = length(Time_periods),
                     Time_periods = paste(Time_periods, collapse = ", "))

write.csv(basicD,
          here::here(output_dir,paste0(sptogo,"_SDM_info.csv")),
          row.names = FALSE)

########################
# Project models to the Study area and time periods

### project for the BG at each time period

for(i in 1:nrow(shift_info)){
    
    projname_start <- paste(sptogo, 
                            shift_info$ID[i], 
                            strsplit(split = "[-]",as.character(shift_info$time_period[i]))[[1]][1],
                            "BG", "start")
    
    m_start <- BIOMOD_Projection(
        bm.mod = model_sp, 
        proj.name = projname_start, 
        new.env = bioclimatics_periods[[i]]$back_start,
        nb.cpu = N_cpus)
    
    pdf(here::here(output_dir,paste(projname_start,"SDM.pdf")),
        width = 10, height = 10)
    plot(m_start)
    dev.off()
    
    ens_m_start <- BIOMOD_EnsembleForecasting(
        bm.em = ens_model_sp, 
        bm.proj = m_start,
        proj.name = paste(projname_start,"ens"),
        nb.cpu = N_cpus)
    
    pdf(here::here(output_dir,paste(projname_start,"ens_SDM.pdf")),
        width = 10, height = 10)
    plot(ens_m_start)
    dev.off()
    
    projname_end <- paste(sptogo, 
                          shift_info$ID[i], 
                          strsplit(split = "[-]",as.character(shift_info$time_period[i]))[[1]][1],
                          "BG", "end")
    
    m_end <- BIOMOD_Projection(
        bm.mod = model_sp, 
        proj.name = projname_end, 
        new.env = bioclimatics_periods[[i]]$back_end, 
        new.env.xy = PresAbs[,c("x","y")],
        nb.cpu = N_cpus)
    
    pdf(here::here(output_dir,paste(projname_end,"SDM.pdf")),
        width = 10, height = 10)
    plot(m_end)
    dev.off()
    
    ens_m_end <- BIOMOD_EnsembleForecasting(
        bm.em = ens_model_sp, 
        bm.proj = m_end,
        proj.name = paste(projname_end,"ens"),
        nb.cpu = N_cpus)
    
    pdf(here::here(output_dir,paste(projname_end,"ens_SDM.pdf")),
        width = 10, height = 10)
    plot(ens_m_end)
    dev.off()
    
}

### project for each SA and time period

for(i in 1:nrow(shift_info)){
    
    projname_start <- paste(sptogo, 
                            shift_info$ID[i], 
                            strsplit(split = "[-]",as.character(shift_info$time_period[i]))[[1]][1],
                            "SA", "start")
    
    m_start <- BIOMOD_Projection(
        bm.mod = model_sp, 
        proj.name = projname_start, 
        new.env = bioclimatics_periods_SA[[i]]$back_start, 
        new.env.xy = PresAbs[,c("x","y")],
        nb.cpu = N_cpus)
    
    pdf(here::here(output_dir,paste(projname_start,"SDM.pdf")),
        width = 10, height = 10)
    plot(m_start)
    dev.off()
    
    ens_m_start <- BIOMOD_EnsembleForecasting(
        bm.em = ens_model_sp, 
        bm.proj = m_start,
        proj.name = paste(projname_start,"ens"),
        nb.cpu = N_cpus)
    
    pdf(here::here(output_dir,paste(projname_start,"ens_SDM.pdf")),
        width = 10, height = 10)
    plot(ens_m_start)
    dev.off()
    
    projname_end <- paste(sptogo, 
                          shift_info$ID[i], 
                          strsplit(split = "[-]",as.character(shift_info$time_period[i]))[[1]][1],
                          "SA", "end")
    
    m_end <- BIOMOD_Projection(
        bm.mod = model_sp, 
        proj.name = projname_end, 
        new.env = bioclimatics_periods_SA[[i]]$back_end, 
        new.env.xy = PresAbs[,c("x","y")],
        nb.cpu = N_cpus)
    
    pdf(here::here(output_dir,paste(projname_end,"SDM.pdf")),
        width = 10, height = 10)
    plot(m_end)
    dev.off()
    
    ens_m_end <- BIOMOD_EnsembleForecasting(
        bm.em = ens_model_sp, 
        bm.proj = m_end,
        proj.name = paste(projname_end,"ens"),
        nb.cpu = N_cpus)
    
    pdf(here::here(output_dir,paste(projname_end,"ens_SDM.pdf")),
        width = 10, height = 10)
    plot(ens_m_end)
    dev.off()
    
    
}
