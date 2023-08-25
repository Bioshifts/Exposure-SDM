# --------------------------------------------------------
# title: "Calculate climate change velocity for each study area in the Bioshifts database"
# author: "Brunno F Oliveira & Bioshifts group"
# --------------------------------------------------------

# Gradient-based climate velocities
# 
# Gradient-based climate velocities (gVoCC) are a way to describe the propensity of a species to move in response to climate change. They represent the ratio between the long-term temporal trend in climate conditions by the spatial gradient in climate conditions:
#     
#     gVoCC = long-term trend / spatial gradient
# 
# A high gVoCC indicates that a species in that location would have a high propensity to move in response to climate change. When the gVoCC is high, it indicates that the climate is changing quickly over time relative to the amount of spatial variation in climate. High gVoCCs tend to occur in areas with heterogeneous climates that are experiencing high rates of temporal change, such as deserts. In areas with high gVoCCs, an organism would have to shift quickly and move far in order to reach an analogous climate.
# 
# On the other hand, a low gVoCC indicates that a species in that location would have a high propensity to move in response to climate change. In these areas, the climate is changing slowly relative to the amount of spatial variation in climate conditions. Areas like mountainous regions, where spatial variation in climates is high, tend to have low gVoCCs. In these areas, on organism would be less inclined to move quickly because a similar climate is closer by.

# Setup
rm(list=ls())
gc()

list.of.packages <- c("terra","elevatr","raster","sf","rgdal","Hmisc","dplyr", "data.table","VoCC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)


# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

# source settings
source("R/settings.R")
source("R/bioshiftsFunction.R")

# detect slurm cores
N_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

# Eckert 4 equal-area projection
Eckt <- enmSdmX::getCRS('Eckert 4')

# create dir to store results
if(!dir.exists(velocity_SA_dir)){
    dir.create(velocity_SA_dir)
}

# Load study areas v3
v3_polygons <- list.files(SA_shps_dir,pattern = ".shp",full.names = TRUE)
v3_polygons_names <- gsub(".shp","",list.files(SA_shps_dir,pattern = ".shp"))

# Load Bioshifts v3
Bioshifts_DB <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB))

# Create Polygon IDs
Bioshifts_DB$Polygon <- paste0("A",Bioshifts_DB$ID_v1,"_",Bioshifts_DB$Study_ID_v1)
pos_v2 <- which(is.na(Bioshifts_DB$ID_v1))
Bioshifts_DB$Polygon[pos_v2] <- paste0("B",Bioshifts_DB$ID_v2[pos_v2],"_",Bioshifts_DB$Study_ID_v2[pos_v2])

# Filter Polygons in Study areas v3
Bioshifts_DB <- Bioshifts_DB %>%
    filter(Polygon %in% v3_polygons_names)

keep_poly <- which(v3_polygons_names %in% Bioshifts_DB$Polygon)
v3_polygons_names <- v3_polygons_names[keep_poly]
v3_polygons <- v3_polygons[keep_poly]

dim(Bioshifts_DB)
length(unique(Bioshifts_DB$Polygon))
length(unique(v3_polygons_names))

########################
# Loop through shape files of study areas

i = 1
# parallel::mclapply(1:length(v3_polygons), function(i) { 

SA_i <- terra::vect(v3_polygons[i])

Bioshifts_i <- Bioshifts_DB %>%
    filter(Polygon == v3_polygons_names[i])

# Is it Terrestrial or Marine?
# Check then load temperature data
my_test <- if(any(is.na(SA_i$EleExtentm))){ # it is terrestrial if it has elevation data 
    ECO = "Mar"
} else {
    ECO = "Ter"
}

# shifts start and end
if(grepl("A",v3_polygons_names[i])){
    S_start <- round(unique(Bioshifts_i$START_v1),0)
    S_end <- round(unique(Bioshifts_i$END_v1),0)
    
} else {
    S_start <- Bioshifts_i$START_v2
    S_end <- Bioshifts_i$END_v2
}
S_time <- S_start:S_end

# get layers within time period of shift
# temperature_layers <- list.files(here::here(vars_dir(ECO),paste0("bio_proj_",my_res)))
# temperature_layers_names <- list.files(here::here(vars_dir(ECO),paste0("bio_proj_",my_res)))
temperature_layers <- list.files(here::here(vars_dir(ECO),"bio_proj_5km"),full.names = TRUE)
temperature_layers_names <- list.files(here::here(vars_dir(ECO),"bio_proj_5km"))
temperature_layers_pos <- grepl(paste(S_time,collapse = "|"),temperature_layers_names)
temperature_layers_names <- temperature_layers_names[temperature_layers_pos]
temperature_layers <- temperature_layers[temperature_layers_pos]
temperature_layers <- terra::rast(temperature_layers)
# get MAT only
temperature_layers <- temperature_layers[[which(names(temperature_layers)=="mat")]]
names(temperature_layers) <- S_time
# crop layers to the study area
temperature_layers <- terra::mask(terra::crop(temperature_layers, SA_i), SA_i)

# project to equal-area
temperature_layers <- terra::project(temperature_layers,Eckt)
temperature_layers <- raster::stack(temperature_layers)

# calculate the trend
ttrend = tempTrend(r = temperature_layers,
                   th = 0.25*nlayers(temperature_layers) ## set minimum # obs. to 1/4 time series length
)

# calculate the spatial gradient
spgrad = spatGrad(r = temperature_layers, 
                  projected = TRUE) ## our raster is projected to a coordinate system

## calculate gradient based climate velocity:
gvocc = gVoCC(tempTrend = ttrend, spatGrad = spgrad)

## calculate resulting velocity North
## radians = degree × π/180
gvocc_lat <- gvocc[[1]] * sin(gvocc[[2]]*pi/180) # convert degrees to radians

# change sign if SA in the south hemisphere to reflect a velocity away of the tropics
max_lat <- ext(SA_i)[4]
min_lat <- ext(SA_i)[3]

max_lat <- 40
min_lat <- 0

if(all(sign(c(max_lat,min_lat)) == 1)){
    prop_south <- 0
    prop_north <- 1
} 
if(all(sign(c(max_lat,min_lat)) == -1)){ 
    prop_south <- 1
    prop_north <- 0
}
if(!all(sign(c(max_lat,min_lat)) == 1)){
    total = max_lat - min_lat
    prop_south <- abs(min_lat/total)
    prop_north <- abs(max_lat/total)
} 
prop_north
prop_south
prop_NS <- prop_north-prop_south

if(sign(prop_NS) == -1){ gvocc_lat <- gvocc_lat * -1 }

# summary velocity over study area
sa_vel <- data.frame(t(data.frame(summary(gvocc_lat))))
rownames(sa_vel) <- NULL
sa_vel <- sa_vel[,1:5]
names(sa_vel) <- c("Min","1stQuant","Median","3rdQuant","Max")    
sa_vel$Mean <- as.numeric(terra::global(terra::rast(gvocc_lat),mean,na.rm = TRUE))
sa_vel$ID = v3_polygons_names[i]

write.csv(sa_vel, here::here(velocity_SA_dir, paste0(v3_polygons_names[i],".csv")), row.names = FALSE)

# }, mc.cores = N_cores)


filePath <- system.file("external/", package="climetrics") # path to the dataset folder
pr <- rast(paste0(filePath,'/precip.tif'))
tmean <- rast(paste0(filePath,'/tmean.tif'))
n <- readRDS(paste0(filePath,'/dates.rds')) # corresoinding dates
head(n) # Dates corresponds to the layers in climate variables
####################
# use rts function in the rts package to make a raster time series:
pr.t <- rts(pr,n)
tmean.t <- rts(tmean,n)
###########################
dv <- gVelocity(tmean.t)
plot(dv,main="temp vel package climetrics")
rmote::plot_done()


# calculate the trend
ttrend = tempTrend(r = stack(tmean),
                   th = 0.25*nlayers(temperature_layers) ## set minimum # obs. to 1/4 time series length
)

# calculate the spatial gradient
spgrad = spatGrad(r = stack(tmean)) 

## calculate gradient based climate velocity:
gvocc = gVoCC(tempTrend = ttrend, spatGrad = spgrad)
plot(rast(gvocc[[1]]),main="temp vel package VoCC")
rmote::plot_done()

plot(dv*1000, rast(gvocc[[1]]), xlab = "climetrics", ylab = "VoCC")
rmote::plot_done()


test <- climetrics:::.tempgradTerra(tmean)
plot(rast(ttrend[[1]]), test, xlab = "VoCC trend", ylab = "climetrics trend")
rmote::plot_done()

test2 <- climetrics:::.spatialgradTerra(tmean)
test2_angle <- rast(spgrad[[1]])
test2_angle[test2$icell] <- test2$angle
plot(rast(spgrad[[2]]), test2_angle, xlab = "VoCC angle", ylab = "climetrics angle")
rmote::plot_done()

test2$NS[is.na(test2$NS)] <- 0
test2$WE[is.na(test2$WE)] <- 0
test2$NAsort <- ifelse((abs(test2$NS)+abs(test2$WE)) == 0, NA, 1)
test2$Grad <- test2$NAsort * sqrt((test2$WE^2) + (test2$NS^2))
test2_Grad <- rast(spgrad[[1]])
test2_Grad[test2$icell] <- test2$Grad
plot(rast(spgrad[[1]]), test2_Grad, xlab = "VoCC Grad", ylab = "climetrics Grad")
rmote::plot_done()



