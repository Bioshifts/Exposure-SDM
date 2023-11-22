# Setup
rm(list=ls())
gc()

# devtools::install_github("sjevelazco/flexsdm@HEAD")
# devtools::install_github("azizka/speciesgeocodeR")

list.of.packages <- c("terra","rnaturalearthdata","biomod2",
                      "tidyverse","tictoc","tidyterra","ggplot2","data.table",
                      "parallelly","speciesgeocodeR")

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

# sptogo="Abra_alba"
# realm <- "Mar"
# sptogo="Arnoglossus_laterna"
# realm <- "Mar"
# 
# sptogo="Abies_alba"
# realm <- "Ter"
# sptogo="Abies_concolor"
# realm <- "Ter"
# sptogo="Formica_sanguinea"
# realm <- "Ter"


########################
# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}
work_dir <- getwd()

########################
# source functions
source("R/my_functions.R")
source("R/get_ecoregions.R")
source("R/correct_colinvar_pca.R")
source("R/plot_SA_location.R")
source("R/mask_bios_BG.R")

########################
# settings
source("R/settings.R")

# N cpus to use
N_cpus <- parallelly::availableCores()

########################
# get directories
bios_SA_dir <- bios_SA_dir(realm)
bios_dir <- bios_dir(realm)
output_dir <- here::here(sdm_dir(realm), sptogo)
occ_dir <- env_data_dir
vars_dir <- vars_dir(realm)

# create dirs
if(!dir.exists(output_dir)){
    dir.create(output_dir,recursive = TRUE)
}

########################
# Load env data for species i
PresAbs <- qs::qread(here::here(occ_dir(realm),paste0(sptogo,"_",realm,".qs")))
PresAbs <- na.omit(PresAbs)

table(PresAbs$pa)

PresAbs <- PresAbs %>%
    filter(mat < 400 & mat > -400)

# filter presences
sp_occ <- PresAbs %>% filter(pa == 1)
# filter pseudo-absences
back_occ <- PresAbs %>% filter(pa == 0)

########################
# Study info

biov1 <- read.csv(here::here(Bioshifts_dir,"biov1_fixednames.csv"))
# v1 from which we have shape files
shp_files <- list.files(SA_shps_dir,pattern = ".shp")
shp_files <- gsub(".shp","",shp_files)
biov1 <- biov1[biov1$ID %in% shp_files,]
# only LAT
biov1$Type[which(biov1$Type=="HOR")] <- "LAT"
biov1 <- biov1[which(biov1$Type=="LAT"),]

biov1 <- biov1 %>% filter(sp_name_std_v1 == sptogo)

# Use LAT ELE shifts
biov1$Type[which(!is.na(biov1$Azimuth))] <- "LAT" # All obs type == HOR contain Azimuth value
biov1 <- biov1[which(biov1$Type=="ELE" | biov1$Type=="LAT"),]

# Use temporal period from the environmental data
biov1 <- biov1 %>% filter(START >= temporal_range_env_data(realm)[1] + n_yr_bioclimatic)

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

# shift duration
S_time <- lapply(1:nrow(shift_info), function(i) {
    shift_info$START[i]:shift_info$END[i]
})
S_time <- unique(round(do.call(c,S_time),0))


########################
# Load ecoregions
cat("Load ecoregions\n")

# load mask raster
mask.ras <- if(realm=="Ter") { 
    terra::rast(here::here(vars_dir,paste0("model_raster_ter_",my_res,".tif")))
} else {
    if(realm=="Mar") { 
        terra::rast(here::here(vars_dir,"model_raster_mar.tif")) 
    }
}

BA <- get_ecoregions(realm = realm, 
                     sptogo = sptogo,
                     PresAbs = unique(sp_occ[,c('x','y')]), # to remove duplicated coordinates (in the same dates) 
                     varsdir = vars_dir, 
                     output_dir = output_dir,
                     mask.ras = mask.ras,
                     return.shp = TRUE,
                     return.raster = FALSE,
                     check_if_exists = TRUE)
# test
# {plot(BA$shape_file)
# plot(vect(data.frame(back_occ), geom = c("x","y")),add=T,col = "red")
# plot(vect(data.frame(sp_occ), geom = c("x","y")),add=T)
# dev.off()}


########################
# Load study area
fgdb <- here::here(SA_shps_dir)
# get study ID
StudyID <- shift_info$ID

StudyArea <- lapply(StudyID, function(x) {
    tmp <- terra::vect(fgdb,layer=x)
})
StudyArea <- vect(StudyArea)
StudyArea$NAME = StudyID

########################
# Plot for study area location
plot_SA_location(sptogo, 
                 BA_shp = BA$shape_file,
                 realm = realm, 
                 sp_occ = sp_occ %>% 
                     dplyr::filter(!duplicated('x','y')), # remove duplicated coordinates (in the same dates) 
                 StudyArea = StudyArea, 
                 maps.dir = output_dir, 
                 show.plot = FALSE)


########################
# Get bioclimatics for the BG at each time period

# select bioclimatics
bioclimatics_BG <- list.files(bios_dir)
bioclimatics_BG_pos <- grepl(paste(S_time,collapse = "|"),bioclimatics_BG)
bioclimatics_BG <- bioclimatics_BG[bioclimatics_BG_pos]

# load in bioclimatics
bioclimatics_BG <- lapply(bioclimatics_BG, function(x) terra::rast(here::here(bios_dir,x)))
names_bioclimatics_BG <- sapply(bioclimatics_BG, function(x){
    tmp <- strsplit(terra::sources(x),"/")[[1]]
    tmp <- tmp[length(tmp)]
    tmp <- strsplit(tmp,"_")[[1]]
    tmp <- tmp[length(tmp)]
    gsub(".tif","",tmp)
})
names(bioclimatics_BG) <- names_bioclimatics_BG


########################
# Mask bioclimatics to the BG
cat("Mask the BG\n")

# do not parallelize this step with all cores - it will overload memory
output_BG <- here::here(output_dir,"BG")
if(!dir.exists(output_BG)){
    dir.create(output_BG)
}

bioclimatics_BG <- mask_bios_BG(bioclimatics_BG, output_BG)


########################
# Load in bioclimatics at the SA 
bioclimatics_SA <- lapply(StudyID, function(x) {
    tmp_names <- list.files(here::here(bios_SA_dir,x))
    tmp_names <- gsub(paste0(c(x,"bios",".tif"," "),collapse = "|"),"",tmp_names)
    tmp <- list.files(here::here(bios_SA_dir,x), full.names = TRUE)
    tmp <- lapply(tmp, terra::rast)
    names(tmp) <- tmp_names
    return(tmp)
})
names(bioclimatics_SA) <- StudyID


#########################
# Reduce collinearity among the predictors
# Run SDMs of PCA axes
cat("Treating collinearity\n")

new_data <- PCA_env(sptogo = sptogo,
                    bioclimatics_BG = bioclimatics_BG,
                    bioclimatics_SA = bioclimatics_SA,
                    output_dir = output_dir,
                    shift_info = shift_info,
                    PresAbs = PresAbs)

# Update data
# BG
bioclimatics_BG <- new_data$bios_BG

# SA
bioclimatics_SA <- new_data$bios_SA

# PresAbs
PresAbsFull <- new_data$bios_PresAbs

rm(new_data)

########################
cat("Fitting SDMs\n")
# P.S.: The code bellow works using the version 4.2-3 of biomod2
# installation :
# devtools::install_github("biomodhub/biomod2", dependencies = TRUE)

# Documentation for functions can be found here: 
# https://biomodhub.github.io/biomod2/index.html

# Fit SDMs

table(PresAbsFull$pa)
PresAbsFull <- PresAbsFull[which(PresAbs$mat < 400 & PresAbs$mat > -400),]

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


# check if model exists
models_sp <- list.files(here::here(output_dir, gsub("_",".",sptogo)),
                        pattern = "ensemble.models.out",full.names = TRUE)

if(length(models_sp)>1){
    delete_duplicated_models(realm = realm,species = sptogo)
}

if(length(models_sp) == 0){
    
    # Fit sdms
    model_sp <- BIOMOD_Modeling(
        data_sp, 
        models = c("GLM","GBM","GAM","MAXNET"),
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
    
} else {
    
    # Load sdms
    files_sdms <- list.files(
        here::here(output_dir,gsub("_",".",sptogo)), 
        full.names = TRUE,
        pattern = "models.out")  
    
    pos <- grep("ensemble",files_sdms)
    files_sdms_ensemble <- files_sdms[pos]
    files_sdms <- files_sdms[-pos]
    
    model_sp <- get(load(files_sdms))
    
    ens_model_sp <- get(load(files_sdms_ensemble))
    
}


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
                     N_background = length(which(PresAbsFull$pa==0)))

write.csv(basicD,
          here::here(output_dir,paste0(sptogo,"_SDM_info.csv")),
          row.names = FALSE)

########################
### Project to the BG at each time period

# for(i in 1:length(bioclimatics_BG)){
#     
#     projname <- paste(sptogo, 
#                       names(bioclimatics_BG)[[i]],
#                       "BG")
#     
#     file_proj <- here::here(output_dir,gsub("_",".",sptogo),
#                             paste0("proj_",projname),
#                             paste0(gsub("_",".",sptogo),".",projname,".projection.out"))
# 
#     projname_ens <- paste(projname,"ens")
#     
#     file_ens <- gsub("BG","BG ens",file_proj)
#     file_ens <- gsub(".projection.out",".ensemble.projection.out",file_ens)
#     
#     
#     if(!file.exists(file_proj)){
#         m <- BIOMOD_Projection(
#             bm.mod = model_sp, 
#             proj.name = projname, 
#             new.env = bioclimatics_BG[[i]],
#             nb.cpu = N_cpus,
#             keep.in.memory = FALSE)
#     }
#     
#     if(!file.exists(file_ens)){
#         ens_m <- BIOMOD_EnsembleForecasting(
#             bm.em = ens_model_sp, 
#             bm.proj = m,
#             proj.name = projname_ens,
#             nb.cpu = N_cpus,
#             keep.in.memory = FALSE)
#     }
#     
#     gc()
#     
# }

### Project to each SA and each time period

for(j in 1:length(bioclimatics_SA)){
    
    bioclimatics_SA_j <- bioclimatics_SA[[j]]
    
    for(i in 1:length(bioclimatics_SA_j)){
        
        projname <- paste(sptogo, 
                          names(bioclimatics_SA)[[j]],
                          names(bioclimatics_SA_j)[[i]],
                          "SA")
        
        file_proj <- here::here(output_dir,gsub("_",".",sptogo),
                                paste0("proj_",projname),
                                paste0(gsub("_",".",sptogo),".",projname,".projection.out"))
        
        projname_ens <- paste(projname,"ens")
        
        file_ens <- gsub("SA","SA ens",file_proj)
        file_ens <- gsub(".projection.out",".ensemble.projection.out",file_ens)
        
        if(!file.exists(file_proj)){
            m <- BIOMOD_Projection(
                bm.mod = model_sp, 
                proj.name = projname, 
                new.env = bioclimatics_SA_j[[i]],
                nb.cpu = N_cpus,
                keep.in.memory = FALSE)
        }
        
        if(!file.exists(file_ens)){
            ens_m <- BIOMOD_EnsembleForecasting(
                bm.em = ens_model_sp, 
                bm.proj = m,
                proj.name = projname_ens,
                nb.cpu = N_cpus,
                keep.in.memory = FALSE)
        }
        
        gc()
        
    }
    
}

# delete temporary data
possible_sdms_SA <- here::here(output_dir,gsub("_",".",sptogo),paste(paste0("proj_",sptogo),shift_info$ID,S_time,"SA"))
possible_sdms_BG <- here::here(output_dir,gsub("_",".",sptogo),paste(paste0("proj_",sptogo),S_time,"BG"))

if(all(dir.exists(possible_sdms_BG))){
    unlink(here::here(output_dir,"BG"), recursive = TRUE)
    unlink(here::here(output_dir,"BG_PC"), recursive = TRUE)
} 

if(all(dir.exists(possible_sdms_SA))){
    unlink(here::here(output_dir,"SA_PC"), recursive = TRUE)
} 
