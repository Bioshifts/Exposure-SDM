# Setup
rm(list=ls())
gc()

# devtools::install_github("sjevelazco/flexsdm@HEAD")

list.of.packages <- c("terra","rnaturalearthdata","biomod2",
                      "tidyverse","tictoc","tidyterra","ggplot2","data.table",
                      "parallelly")

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
# sptogo="Acacia_melanoxylon"
# realm <- "Ter"
# sptogo="Abies_concolor"
# realm <- "Ter"

# N cpus to use
N_cpus <- parallelly::availableCores()

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

########################
# source functions
source("R/my_functions.R")
source("R/get_ecoregions.R")
source("R/correct_colinvar_pca.R")
source("R/plot_SA_location.R")

########################
# settings
# source settings
source("R/settings.R")

check_if_exists = FALSE

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

mask.ras <- if(realm=="Ter") { 
    terra::rast(here::here(vars_dir,paste0("model_raster_ter_",my_res,".tif")))
} else {
    if(realm=="Mar") { 
        terra::rast(here::here(vars_dir,"model_raster_mar.tif")) 
    }
}

########################
# Check if file exists

if(check_if_exists){
    
    # check if dir with sdm output exists
    if(dir.exists(here::here(output_dir,gsub("_",".",sptogo)))){
        stop("No need to run this. SDMs were fitted for this species!")
    } 
}

########################
# Load env data for species i
PresAbs <- qs::qread(here::here(occ_dir(realm),paste0(sptogo,"_",realm,".qs")))
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

biov1 <- read.csv(here::here(Bioshifts_dir,"biov1_fixednames.csv"))
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
# create dir if doesn't exist
plot_SA_location(sptogo, 
                 BA_shp,
                 realm, 
                 PresAbs, 
                 StudyArea, 
                 maps.dir = output_dir, 
                 show.plot = FALSE)

########################
# Get bioclimatics for the BG at each time period

S_time <- lapply(1:nrow(shift_info), function(i) {
    shift_info$START[i]:shift_info$END[i]
})
S_time <- unique(round(do.call(c,S_time),0))

bioclimatics_BG <- list.files(bios_dir)
bioclimatics_BG_pos <- grepl(paste(S_time,collapse = "|"),bioclimatics_BG)
bioclimatics_BG <- bioclimatics_BG[bioclimatics_BG_pos]
bioclimatics_BG <- lapply(bioclimatics_BG, function(x) terra::rast(here::here(bios_dir,x)))
names(bioclimatics_BG) <- S_time

# mask
if(!dir.exists(here::here(output_dir,"BG"))){
    dir.create(here::here(output_dir,"BG"))
}

bioclimatics_BG <- parallel::mclapply(1:length(bioclimatics_BG), function(x) {
    bios_year_i = bioclimatics_BG[[x]]
    year_i = names(bioclimatics_BG)[x]
    terra::window(bios_year_i) <- terra::ext(BA_shp)
    terra::mask(bios_year_i, BA_shp,
                filename = here::here(output_dir,"BG",paste0(year_i,".tif")),
                overwrite = TRUE)
}, mc.cores = N_cpus)
# load in
bioclimatics_BG <- lapply(list.files(here::here(output_dir,"BG"),full.names = TRUE), 
                          function(x) terra::rast(x))
names(bioclimatics_BG) <- S_time

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

# ########################
# Reduce collinearity among the predictors
# Run SDMs of PCA axes

# create an assemble of environmental data coming from species occurrences + background data at each time period

# 1) add background data 
data_env <- as.data.frame(bioclimatics_BG[[1]], na.rm=TRUE)
if(nrow(data_env) > 20000){
    data_env <- data_env[sample(1:nrow(data_env),20000),]
}

# 2) add pres/abs data
data_env <- rbind(data_env, 
                  data.frame(PresAbs)[,names(data_env)])

dim(data_env)

# set a cap for the size of dataset
mycap <- 100000
if(nrow(data_env) > mycap){
    data_env <- data_env[sample(1:nrow(data_env),mycap), ]
}

# 3) get uncorrelated PCs
new_data <- correct_colinvar_pca(env_layer = data_env, 
                                 proj = bioclimatics_BG, 
                                 proj2 = bioclimatics_SA, 
                                 PresAbs = PresAbs)

# 4) update data
bioclimatics_BG <- new_data$proj
bioclimatics_SA <- new_data$proj2
PresAbsFull <- data.table(data.frame(PresAbs[,c("x","y","pa")], new_data$PresAbs))
PresAbsFull <- PresAbsFull[order(PresAbsFull$pa, decreasing = TRUE),]


########################

# P.S.: The code bellow works using the version 4.2-3 of biomod2
# installation :
# devtools::install_github("biomodhub/biomod2", dependencies = TRUE)

# Documentation for functions can be found here: 
# https://biomodhub.github.io/biomod2/index.html

# Fit SDMs

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
                     N_background = length(which(PresAbsFull$pa==0)))

write.csv(basicD,
          here::here(output_dir,paste0(sptogo,"_SDM_info.csv")),
          row.names = FALSE)

########################
# Project models to the Study area and time periods

### project at the BG at each time period

for(i in 1:length(bioclimatics_BG)){
    
    projname <- paste(sptogo, 
                      names(bioclimatics_BG)[[i]],
                      "BG")
    
    m <- BIOMOD_Projection(
        bm.mod = model_sp, 
        proj.name = projname, 
        new.env = bioclimatics_BG[[i]],
        nb.cpu = N_cpus)
    
    ens_m <- BIOMOD_EnsembleForecasting(
        bm.em = ens_model_sp, 
        bm.proj = m,
        proj.name = paste(projname,"ens"),
        nb.cpu = N_cpus)
    
}

### project for each SA and time period

for(j in 1:length(bioclimatics_SA)){
    
    bioclimatics_SA_j <- bioclimatics_SA[[j]]
    
    for(i in 1:length(bioclimatics_SA_j)){
        
        projname <- paste(sptogo, 
                          shift_info$ID[j], 
                          names(bioclimatics_SA)[[j]],
                          names(bioclimatics_SA_j)[[i]],
                          "SA")
        
        m <- BIOMOD_Projection(
            bm.mod = model_sp, 
            proj.name = projname, 
            new.env = bioclimatics_SA_j[[i]],
            nb.cpu = N_cpus)
        
        ens_m <- BIOMOD_EnsembleForecasting(
            bm.em = ens_model_sp, 
            bm.proj = m,
            proj.name = paste(projname,"ens"),
            nb.cpu = N_cpus)
        
    }
    
}


# delete saved bioclimatics at the background area
unlink(here::here(output_dir,"BG"), recursive = TRUE)
