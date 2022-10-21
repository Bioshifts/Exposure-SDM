# Setup
rm(list=ls())
gc()

# devtools::install_github("sjevelazco/flexsdm@HEAD")

list.of.packages <- c("rgdal","terra","rgis","rasterVis",
                      "PresenceAbsence","flexsdm","SDMtune",
                      "tidyverse","tictoc","tidyterra","ggplot2","data.table",
                      "RSQLite","DBI","odbc",
                      "foreach","doParallel","pbapply","doRNG")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

########################
# load functions
source(here::here("R/my_functions.R"))
source(here::here("R/MaxNet_functions.R"))
source(here::here("R/get_ecoregions.R"))
source(here::here("R/get_env_data_for_modeling.R"))
source(here::here("R/plot_SA_location.R"))
source(here::here("R/get_bioclimatics_period.R"))
source(here::here("R/correct_colinvar_pca.R"))

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
sptogo <- as.character(paste(command_args[1], collapse = " "))
realm <- as.character(paste(command_args[2], collapse = " "))

# sptogo="Amblyraja_radiata"
# sptogo="Hemilepidotus_jordani"

# sptogo="Asterias_amurensis"
# realm <- "Mar"

########################
# source settings
source(here::here("R/settings.R"))

n_jobs <- 20

output_dir = here::here("/media/seagate/boliveira/SDMs/MaxNet",realm,sptogo)
if(!dir.exists(output_dir)){
    dir.create(output_dir,recursive = T)
}

########################
# Check if file exists
file.test <- here::here(output_dir,paste0(sptogo,"_model.RDS"))
if(file.exists(file.test)){
    stop("No need to run this. SDMs were fitted for this species!")
} 

########################
# Load env data for species i
PresAbs <- readRDS(here::here("/media/seagate/boliveira/SDMs/Env_data",realm,paste0(sptogo,"_",realm,"_bio.RDS")))
PresAbs <- na.omit(PresAbs)
# filter presences
sp_occ <- PresAbs %>% filter(pa == 1)
sp_occ <- if(nrow(sp_occ) > limit_recs) {sp_occ[sample(1:nrow(sp_occ), limit_recs),]} else {sp_occ}
# filter pseudo-absences
back_occ <- PresAbs %>% filter(pa == 0)
back_occ <- PresAbs %>% filter(!cell %in% sp_occ$cell) # Remove pseudo-absences at the same cells with presences
back_occ <- if(nrow(back_occ) > limit_recs) {back_occ[sample(1:nrow(back_occ), limit_recs),]} else {back_occ}

PresAbs <- rbind(sp_occ, back_occ)

table(PresAbs$pa)

########################
# Load environmental variables
env_data <- get_env_data_for_modeling(realm,myvars,varsdir)

########################
# Load ecoregions and more
BA <- get_ecoregions(realm,PresAbs)

# test
# plot(BA)
# plot(vect(data.frame(PresAbs[which(PresAbs$pa==0),]), geom = c("x","y")),add=T)
# plot(vect(data.frame(PresAbs[which(PresAbs$pa==1),]), geom = c("x","y")),add=T,col = "red")

########################
# Study info

biov1 <- read.table(here::here("Data/Shifts2018_checkedtaxo.txt"),header = T)
biov1 <- biov1 %>% filter(New_name == sptogo)

# Use LAT ELE shifts
biov1$Type[which(!is.na(biov1$Azimuth))] <- "LAT" # All obs type == HOR contain Azimuth value
biov1 <- biov1[which(biov1$Type=="ELE" | biov1$Type=="LAT"),]

# Use temporal period from the environmental data
biov1 <- biov1 %>% filter(START >= temporal_range_env_data[1] + n_yr_bioclimatic)

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
fgdb <- here::here("Data/Study_Areas_v1/Study_Areas.gdb")

StudyArea <- lapply(StudyID, function(x) {
    tmp <- readOGR(dsn=fgdb,layer=x)
    vect(tmp)
})
StudyArea <- vect(StudyArea)
StudyArea$NAME = StudyID

########################
# Calculate bioclimatics for the BG at each time period
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
        
        tmp <- get_bioclimatics_period(p_1 = p_i[1],
                                       p_2 = p_i[2],
                                       n_yr_bioclimatic = n_yr_bioclimatic,
                                       BA = BA,
                                       env_data = env_data, 
                                       n_cores = n_jobs)
        
        # rasterize 
        t1 <- terra::rast(data.frame(tmp$xy,
                                     tmp$back_start))
        t2 <- terra::rast(data.frame(tmp$xy,
                                     tmp$back_end))
        
        tmp <- list(back_start = t1,
                    back_end = t2)
        
        bioclimatics_periods[[i]] <- tmp
        
        # save
        writeRaster(bioclimatics_periods[[i]]$back_start,
                    here::here(output_dir,paste(sptogo,Time_periods[i],"BG start bios.tif")), overwrite=TRUE)
        writeRaster(bioclimatics_periods[[i]]$back_end,
                    here::here(output_dir,paste(sptogo,Time_periods[i],"BG end bios.tif")), overwrite=TRUE)
        
    }
}
names(bioclimatics_periods) <- Time_periods

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
        tmp1 <- crop(mask(tmp1,
                          StudyArea[StudyArea$NAME==shift_info$ID[i]]),
                     StudyArea[StudyArea$NAME==shift_info$ID[i]])
        
        tmp2 <- bioclimatics_periods[[shift_info$time_period[i]]]$back_end
        tmp2 <- crop(mask(tmp2,
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
# Reduce collinearity among the predictors

data_clean <- data.table(
    data.frame(
        PresAbs)[,-which(
            names(PresAbs) %in% 
                c("cell","x","y","year","month","tile","pa"))])

# create a SDMtune file
data_test <- SDMtune:::SWD(species = sptogo, 
                           coords = PresAbs[,c('x','y')], 
                           data = data_clean, 
                           pa = PresAbs$pa)
# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(data_test, test = 0.2, only_presence = TRUE)
train_data <- datasets[[1]]
test_data <- datasets[[2]]

# train a model using default parameters
test_m <- SDMtune::train(method = "Maxnet", data = train_data)
# Prepare background locations to test autocorrelation
bg <- SDMtune:::SWD(species = sptogo, 
                    coords = PresAbs[PresAbs$pa==0,c('x','y')], 
                    data = data_clean[PresAbs$pa==0,], 
                    pa = PresAbs$pa[PresAbs$pa==0])
# Remove variables with correlation higher than 0.7 accounting for the AUC,
# in the following example the variable importance is computed as permutation
# importance
vs <- SDMtune::varSel(test_m, 
                      metric = "auc", 
                      bg4cor = bg, 
                      test = test_data, 
                      cor_th = 0.7,
                      permut = 1)
vs <- names(vs@data@data)

# Update PresAbs
PresAbs <- data.table(
    data.frame(
        data.frame(PresAbs)[,which(names(PresAbs) %in% c("cell","x","y","year","month","tile","pa"))], 
        data.frame(data_clean)[,vs]))

# Update bioclimatics_periods BG
for(i in 1:length(bioclimatics_periods)){
    bioclimatics_periods[[i]]$back_start <- bioclimatics_periods[[i]]$back_start[[vs]]
    bioclimatics_periods[[i]]$back_end <- bioclimatics_periods[[i]]$back_end[[vs]]
}

# Update bioclimatics_periods SA
for(i in 1:length(bioclimatics_periods_SA)){
    bioclimatics_periods_SA[[i]]$back_start <- bioclimatics_periods_SA[[i]]$back_start[[vs]]
    bioclimatics_periods_SA[[i]]$back_end <- bioclimatics_periods_SA[[i]]$back_end[[vs]]
}


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
# Fit MaxNet models

PresAbsFull <- data.table(
    data.frame(
        PresAbs)[,-which(
            names(PresAbs) %in% 
                c("cell","x","y","year","month","tile"))])
PresAbsFull <- PresAbsFull[order(PresAbsFull$pa, decreasing = TRUE),]

table(PresAbsFull$pa)

# Optimize model
regMult = seq(0.5, 4, by = 0.5)
classes = c("l","h","lq","lqh","lqhp")

# For testing
# PresAbsFull <- PresAbsFull[c(sample(which(PresAbsFull$pa==1),100), sample(which(PresAbsFull$pa==0),100)),]
# regMult = seq(0.5, 4, by = 0.5)[1:2]
# classes = c("l","h","lq","lqh","lqhp")[1:2]

tic()
MaxNet_output <- trainMaxNet(
    spname = sptogo,
    data=data.frame(PresAbsFull), 
    verbose=FALSE, 
    kfolds = NULL, # k-fold cross validation using k-1 folds for training and 1 fold for testing - equivalent to leaving 20% for testing
    regMult = regMult, 
    classes = classes, 
    out = c("model","tuning"),
    cores = n_jobs,
    output_dir = output_dir,
    check_if_exist = TRUE)
toc()
# 10000 >> 0.48 >> 66 sec >> 40.902 sec elapsed
# 11000 >> 0.94 >> 48.427 sec elapsed
# 12000 >> 0.51 >> 60.439 sec elapsed
# 15000 >> 1.81 >> 91.155 sec elapsed
# 50000 >> 449.4536 sec elapsed == 7.5 min
# 100000 >> 961.5964 sec elapsed == 16 min

# View(MaxNet_output$tuning)

########################
# run cross validation
CV_model <- crossvalSDM(kfolds = 5, 
                        traindata = data.frame(PresAbsFull[,-1]), 
                        pa = as.vector(PresAbsFull$pa), 
                        classes = as.character(MaxNet_output$tuning$classes[1]), 
                        regmult = as.numeric(MaxNet_output$tuning$regMult[1]))

CV_eval <- evalSDM(observation = as.vector(PresAbsFull$pa), 
                   predictions = CV_model)

tmp <- MaxNet_output$tuning[1,]
names(tmp) <- gsub("_Full","",names(tmp))
tmp <- tmp[,names(CV_eval)]

CV_eval <- rbind(data.frame(Mod = "Train", tmp),
                 data.frame(Mod = "CV", CV_eval))

write.csv(CV_eval,
          here::here(output_dir,paste0(sptogo,"_SDM_CV.csv")),
          row.names = FALSE)

########################
# save basic data
basicD <- data.frame(Species = gsub("_"," ",sptogo),
                     N_occ = length(which(PresAbsFull$pa==1)),
                     N_background = length(which(PresAbsFull$pa==0)),
                     N_study_areas = length(StudyArea),
                     N_time_periods = length(Time_periods),
                     Time_periods = paste(Time_periods, collapse = ", "),
                     env_vars = paste(vs, collapse = ", "))

write.csv(basicD,
          here::here(output_dir,paste0(sptogo,"_SDM_info.csv")),
          row.names = FALSE)

########################
# Project models to the Study area and time periods

### project for the BG at each time period

mc.cores <- ifelse(nrow(shift_info) > n_jobs, n_jobs, nrow(shift_info))

mymodels <- lapply(1:nrow(shift_info), function(i){
    xys <- bioclimatics_periods[[i]]$back_start
    m_start <- predictMaxNet(MaxNet_output$model, 
                             bioclimatics_periods[[i]]$back_start, 
                             type = "logistic",clamp = FALSE)
    m_start <- rast(na.omit(data.frame(crds(xys, df=FALSE, na.rm=TRUE),m_start)),type='xyz')
    names(m_start) <- paste(sptogo, 
                            shift_info$ID[i], 
                            strsplit(split = "[-]",as.character(shift_info$time_period[i]))[[1]][1])
    
    xys <- bioclimatics_periods[[i]]$back_end
    m_end <- predictMaxNet(MaxNet_output$model, 
                           bioclimatics_periods[[i]]$back_end, 
                           type = "logistic",clamp = FALSE)
    m_end <- rast(na.omit(data.frame(crds(xys, df=FALSE, na.rm=TRUE),m_end)),type='xyz')
    
    names(m_end) <- paste(sptogo, 
                          shift_info$ID[i], 
                          strsplit(split = "[-]",as.character(shift_info$time_period[i]))[[1]][2])
    
    return(c(m_start,m_end))
})
names(mymodels) <- paste(shift_info$ID, shift_info$time_period)

########## write rasters
plot(mymodels[[1]])

for(i in 1:length(mymodels)) {
    writeRaster(mymodels[[i]], here::here(output_dir,paste0(sptogo," BG ",names(mymodels)[i],".tif")),overwrite=TRUE)
    
    pdf(here::here(output_dir,paste0(sptogo," BG ",names(mymodels)[i],"_SDM.pdf")),
        width = 10, height = 5)
    par(mfrow=c(1,2))
    plot(mymodels[[i]][[1]], main = names(mymodels[[i]][[1]]))
    plot(vect(rnaturalearth::ne_coastline(scale = 110, returnclass = "sp")), add=T)
    plot(mymodels[[i]][[2]], main = names(mymodels[[i]][[2]]))
    plot(vect(rnaturalearth::ne_coastline(scale = 110, returnclass = "sp")), add=T)
    dev.off()
}

### project for each SA and time period

mymodels <- lapply(1:nrow(shift_info), function(i){
    xys <- bioclimatics_periods_SA[[i]]$back_start
    m_start <- predictMaxNet(MaxNet_output$model, 
                             bioclimatics_periods_SA[[i]]$back_start, 
                             type = "logistic",clamp = FALSE)
    m_start <- rast(na.omit(data.frame(crds(xys, df=FALSE, na.rm=TRUE),m_start)),type='xyz')
    names(m_start) <- paste(sptogo, 
                            shift_info$ID[i], 
                            strsplit(split = "[-]",as.character(shift_info$time_period[i]))[[1]][1])
    
    xys <- bioclimatics_periods_SA[[i]]$back_end
    m_end <- predictMaxNet(MaxNet_output$model, 
                           bioclimatics_periods_SA[[i]]$back_end, 
                           type = "logistic",clamp = FALSE)
    m_end <- rast(na.omit(data.frame(crds(xys, df=FALSE, na.rm=TRUE),m_end)),type='xyz')
    names(m_end) <- paste(sptogo, 
                          shift_info$ID[i], 
                          strsplit(split = "[-]",as.character(shift_info$time_period[i]))[[1]][2])
    
    return(c(m_start,m_end))
})
names(mymodels) <- paste(shift_info$ID, shift_info$time_period)

########## write rasters
plot(mymodels[[1]])

for(i in 1:length(mymodels)) {
    writeRaster(mymodels[[i]], here::here(output_dir,paste0(sptogo," SA ",names(mymodels)[i],".tif")),overwrite=TRUE)
    pdf(here::here(output_dir,paste0(sptogo," SA ",names(mymodels)[i],"_SDM.pdf")),
        width = 10, height = 5)
    par(mfrow=c(1,2))
    plot(mymodels[[i]][[1]], main = names(mymodels[[i]][[1]]))
    plot(vect(rnaturalearth::ne_coastline(scale = 110, returnclass = "sp")), add=T)
    plot(mymodels[[i]][[2]], main = names(mymodels[[i]][[2]]))
    plot(vect(rnaturalearth::ne_coastline(scale = 110, returnclass = "sp")), add=T)
    dev.off()
}

