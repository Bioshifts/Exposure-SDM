# Setup
rm(list=ls())
gc()

# devtools::install_github("sjevelazco/flexsdm@HEAD")
# devtools::install_github("azizka/speciesgeocodeR")
# devtools::install_github("rpatin/gbm")
# devtools::install_github("mrmaxent/maxnet")
# devtools::install_github("biomodhub/biomod2", dependencies = FALSE)

list.of.packages <- c("terra","rnaturalearthdata","biomod2","gbm","maxnet",
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

# sptogo="Abies_alba"
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

em.by = "all"


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

# output dir
output_dir <- here::here(sdm_dir(realm), sptogo)
if(!dir.exists(output_dir)){
    dir.create(output_dir,recursive = TRUE)
}

########################
# Load env data for species i
PresAbs <- qs::qread(here::here(env_data_dir(realm),paste0(sptogo,"_",realm,".qs")))

table(PresAbs$pa)

# PresAbs <- PresAbs %>%
#     filter(mat < 400 & mat > -400)

# filter presences
sp_occ <- PresAbs %>% filter(pa == 1)
# filter pseudo-absences
back_occ <- PresAbs %>% filter(pa == 0)

########################
# Study info
bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_sdms_selection(bioshifts)
bioshifts <- filter(bioshifts, sp_name_std == sptogo)

# save
write.csv(bioshifts, here::here(output_dir,paste(sptogo,"bioshift.csv")),row.names = FALSE)

# get info relevant for SDM  
shift_info<- bioshifts[,c("ID","Start","End")]
if(any(duplicated(shift_info))){
    shift_info<- shift_info[-which(duplicated(shift_info)),]
}
shift_info$time_period <- paste(round(shift_info$Start), round(shift_info$End), sep = "-")
shift_info$Species <- sptogo
#save
write.csv(shift_info, here::here(output_dir,paste(sptogo,"shift_info.csv")),row.names = FALSE)

# shift duration
S_time <- lapply(1:nrow(shift_info), function(i) {
    shift_info$Start[i]:shift_info$End[i]
})
S_time <- unique(round(do.call(c,S_time),0))

########################
# Load ecoregions
cat("Load ecoregions\n")

# load mask raster
mask.ras <- if(realm=="Ter") { 
    terra::rast(here::here(vars_dir(realm),paste0("model_raster_ter_",my_res,".tif")))
} else {
    if(realm=="Mar") { 
        terra::rast(here::here(vars_dir(realm),"model_raster_mar.tif")) 
    }
}

BA <- get_ecoregions(realm = realm, 
                     sptogo = sptogo,
                     PresAbs = unique(sp_occ[,c('x','y')]), # to remove duplicated coordinates (in the same dates) 
                     varsdir = vars_dir(realm), 
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
bioclimatics_BG <- list.files(bios_dir(realm))
bioclimatics_BG_pos <- grepl(paste(S_time,collapse = "|"),bioclimatics_BG)
bioclimatics_BG <- sort(bioclimatics_BG[bioclimatics_BG_pos])

# load in bioclimatics
bioclimatics_BG <- lapply(bioclimatics_BG, function(x) terra::rast(here::here(bios_dir(realm),x)))
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
bioclimatics_BG <- mask_bios_BG(bioclimatics_BG, output_BG, BA$shape_file)


########################
# Load in bioclimatics at the SA 
bioclimatics_SA <- lapply(StudyID, function(x) {
    tmp_names <- list.files(here::here(bios_SA_dir(realm),x))
    tmp_names <- gsub(paste0(c(x,"bios",".tif"," "),collapse = "|"),"",tmp_names)
    tmp <- list.files(here::here(bios_SA_dir(realm),x), full.names = TRUE)
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
                    which_bioclimatics_BG = 1,
                    output_dir = output_dir,
                    shift_info = shift_info,
                    PresAbs = PresAbs,
                    check_if_PCA_model_exists = TRUE,
                    N_cpus = N_cpus)
names(new_data)

pc_names_bg <- names(new_data$bios_BG[[1]])
pc_names_sa <- names(new_data$bios_SA[[1]][[1]])
pc_names_pa <- names(new_data$bios_PresAbs[,-1:-3])

test1 <- identical(pc_names_bg,pc_names_sa)
test2 <- identical(pc_names_sa,pc_names_pa)
test3 <- identical(pc_names_bg,pc_names_pa)

if(!all(test1, test2, test3)){
    new_data <- PCA_env(sptogo = sptogo,
                        bioclimatics_BG = bioclimatics_BG,
                        bioclimatics_SA = bioclimatics_SA,
                        which_bioclimatics_BG = 1,
                        output_dir = output_dir,
                        shift_info = shift_info,
                        PresAbs = PresAbs,
                        check_if_PCA_model_exists = FALSE,
                        check_if_data_exists = FALSE,
                        N_cpus = N_cpus)
}

# save PCA results
saveRDS(new_data$PCA_model, here::here(output_dir,"PCA_model.RDS"))
write.csv(new_data$coefficients, here::here(output_dir,"PCA_coefficients.csv"), row.names = FALSE)
write.csv(new_data$cumulative_variance, here::here(output_dir,"PCA_cumulative_variance.csv"), row.names = FALSE)

# Update data
# BG
bioclimatics_BG <- new_data$bios_BG
# SA
bioclimatics_SA <- new_data$bios_SA
# PresAbs
PresAbsFull <- na.omit(new_data$bios_PresAbs)

dim(new_data$bios_PresAbs)
dim(PresAbs)
# clean
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

resp <- as.vector(PresAbsFull$pa)
resp[resp==0] <- NA

data_sp <- BIOMOD_FormatingData(
    resp.name = sptogo, 
    resp.var = resp, 
    expl.var = data.frame(PresAbsFull[,-1:-3]), 
    resp.xy = data.frame(PresAbsFull[,1:2]),
    dir.name = output_dir,
    na.rm = TRUE,
    PA.strategy = "random",
    PA.nb.rep = 1,
    PA.nb.absences = length(which(is.na(resp))))


# check if model exists
try({
    delete_duplicated_models(realm = realm, species = sptogo)
}, silent = TRUE)

models_sp <- list.files(here::here(output_dir, gsub("_",".",sptogo)),
                        pattern = "ensemble.models.out",full.names = TRUE)

# model exists?
if(length(models_sp)>0){
    
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
    
    model_sp <- apply_gsub_to_s4(model_sp, pattern = "/lustre/oliveirab", replacement = "/scratch/boliveira")
    ens_model_sp <- apply_gsub_to_s4(ens_model_sp, pattern = "/lustre/oliveirab", replacement = "/scratch/boliveira")
    
    if(!ens_model_sp@em.by == em.by){
        ens_model_sp <- BIOMOD_EnsembleModeling(
            model_sp, 
            em.by = em.by,
            metric.select = "TSS",
            metric.eval = "TSS")
    }
    
} else {
    
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
        em.by = em.by,
        metric.select = "TSS",
        metric.eval = "TSS")
    
}


########################
# Evaluate
myevals = get_evaluations(model_sp, metric.eval = "TSS")
myevals

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

myevals_ens = get_evaluations(ens_model_sp)
myevals_ens

if(em.by=="all"){
    write.csv(myevals_ens, 
              here::here(output_dir,paste0(sptogo,"_CV_ens_all.csv")),
              row.names = FALSE)
} else {
    write.csv(myevals_ens, 
              here::here(output_dir,paste0(sptogo,"_CV_ens_RUN.csv")),
              row.names = FALSE)
}


########################
# save basic data
basicD <- data.frame(Species = sptogo,
                     N_occ = length(which(PresAbsFull$pa==1)),
                     N_background = length(which(PresAbsFull$pa==0)))

write.csv(basicD,
          here::here(output_dir,paste0(sptogo,"_SDM_info.csv")),
          row.names = FALSE)

### Project to each SA and each time period

for(j in 1:length(bioclimatics_SA)){
    
    bioclimatics_SA_j <- bioclimatics_SA[[j]]
    
    if(length(bioclimatics_SA_j)>0){
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
                    build.clamping.mask = FALSE,
                    keep.in.memory = FALSE)
            }
            
            if(!file.exists(file_ens)){
                ens_m <- BIOMOD_EnsembleForecasting(
                    bm.em = ens_model_sp, 
                    bm.proj = m,
                    proj.name = projname_ens,
                    nb.cpu = N_cpus,
                    keep.in.memory = FALSE)
            } else{
                
                test <- get(load(file_ens))
                test <- test@models.projected == ens_model_sp@em.computed
                
                if(!test){
                    ens_m <- BIOMOD_EnsembleForecasting(
                        bm.em = ens_model_sp, 
                        bm.proj = m,
                        proj.name = projname_ens,
                        nb.cpu = N_cpus,
                        keep.in.memory = FALSE)
                }
            }
        }
    }
}

########################
### Project to the BG at t1

projname <- paste(sptogo,
                  names(bioclimatics_BG)[[1]],
                  "BG")

file_proj <- here::here(output_dir,gsub("_",".",sptogo),
                        paste0("proj_",projname),
                        paste0(gsub("_",".",sptogo),".",projname,".projection.out"))

projname_ens <- paste(projname,"ens")

file_ens <- gsub("BG","BG ens",file_proj)
file_ens <- gsub(".projection.out",".ensemble.projection.out",file_ens)


if(!file.exists(file_proj)){
    m <- BIOMOD_Projection(
        bm.mod = model_sp,
        proj.name = projname,
        new.env = bioclimatics_BG[[1]],
        nb.cpu = N_cpus,
        build.clamping.mask = FALSE,
        keep.in.memory = FALSE)
} else {
    m <- get(load(file_proj))
}

if(!file.exists(file_ens)){
    ens_m <- BIOMOD_EnsembleForecasting(
        bm.em = ens_model_sp,
        bm.proj = m,
        proj.name = projname_ens,
        nb.cpu = N_cpus,
        keep.in.memory = FALSE)
} else{
    test <- get(load(file_ens))
    test <- test@models.projected == ens_model_sp@em.computed
    
    if(!test){
        ens_m <- BIOMOD_EnsembleForecasting(
            bm.em = ens_model_sp, 
            bm.proj = m,
            proj.name = projname_ens,
            nb.cpu = N_cpus,
            keep.in.memory = FALSE)
    }
}

gc()

### Clean files ###
unlink(here::here(output_dir,"BG"), recursive = TRUE)
unlink(here::here(output_dir,"BG_PC"), recursive = TRUE)
unlink(here::here(output_dir,"SA_PC"), recursive = TRUE)

to_del <- list.files(here(output_dir,gsub("_",".",sptogo)), full.names = TRUE)
sel <- grep("models",to_del)
to_del <- to_del[-sel]
sel <- grep("SA ens",to_del)
to_del <- to_del[-sel]
sel <- grep("BG ens",to_del)
to_del <- to_del[-sel]

unlink(to_del, recursive = TRUE)

### END ###