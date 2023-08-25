# devtools::install_github("cbrown5/vocc")
# devtools::install_github("JorGarMol/VoCC", dependencies = TRUE, build_vignettes = FALSE)

list.of.packages <- c("climetrics","vocc","VoCC","terra","ggplot2","GGally","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

# Source velocity functions I adapted from the package climetrics after applying some corrections described below 
source(here::here("R/velocity_functions.R"))

# Use the same set of climate variables
filePath <- system.file("external/", package="climetrics") # path to the dataset folder
tmean <- rast(paste0(filePath,'/tmean.tif'))
n <- readRDS(paste0(filePath,'/dates.rds')) # corresponding dates

# project to equal-area
# Eckert 4 equal-area projection
Eckt <- enmSdmX::getCRS('Eckert 4')
tmean <- terra::project(tmean,Eckt)

# Use rts function in the rts package to make a raster time series:
tmean.t <- rts(tmean,n)

##############
## Climetrics
# I found this package is based on the package vocc (not VoCC!), and that both packages have an issue when calculating the spatial gradient. The issue is that when extracting environmental data from the RasterStack  (in vooc, or SpatRasterTS in climetrics), both package use only the first year to calculate the spatial gradient. Instead, the package VoCC uses the mean of all years (calc(r,mean)), which is the correct thing to do! To prove this issue, I calculate velocity with this package using the both the default (wrong) method in this package, and with the correct procedure (using mean raster).


# trend
trend_climatrics <- climetrics:::.tempgradTerra(tmean)

### Wrong way

start <- Sys.time()

# Velocity
v_climetrics1 <- climetrics:::gVelocity(tmean.t)

climetrics_time <- Sys.time() - start

# spatial gradient (wrong way)
spgrad_out_climetrics1 <- climetrics:::.spatialgradTerra(tmean)

spgrad_climetrics1 <- data.frame(NS = spgrad_out_climetrics1$NS, WE = spgrad_out_climetrics1$WE)
spgrad_climetrics1$NS[is.na(spgrad_climetrics1$NS)] <- 0
spgrad_climetrics1$WE[is.na(spgrad_climetrics1$WE)] <- 0
spgrad_climetrics1$NAsort <- ifelse((abs(spgrad_climetrics1$NS)+abs(spgrad_climetrics1$WE)) == 0, NA, 1)
spgrad_climetrics1 <- spgrad_climetrics1$NAsort * sqrt((spgrad_climetrics1$WE^2) + (spgrad_climetrics1$NS^2))
# Angles
angle_climatrics1 <- spgrad_out_climetrics1$angle


### Right way
# spatial gradient
spgrad_out_climetrics2 <- climetrics:::.spatialgradTerra(mean(tmean))

spgrad_climetrics2 <- data.frame(NS = spgrad_out_climetrics2$NS, WE = spgrad_out_climetrics2$WE)
spgrad_climetrics2$NS[is.na(spgrad_climetrics2$NS)] <- 0
spgrad_climetrics2$WE[is.na(spgrad_climetrics2$WE)] <- 0
spgrad_climetrics2$NAsort <- ifelse((abs(spgrad_climetrics2$NS)+abs(spgrad_climetrics2$WE)) == 0, NA, 1)
spgrad_climetrics2 <- spgrad_climetrics2$NAsort * sqrt((spgrad_climetrics2$WE^2) + (spgrad_climetrics2$NS^2))
# Angles
angle_climatrics2 <- spgrad_out_climetrics2$angle

# Velocity
v_climetrics2 <- climetrics:::.calcvelocity(grad = spgrad_out_climetrics2,
                                            slope = trend_climatrics,
                                            .terra = TRUE)

##############
## vocc
# calculate the trend
# default divisor = 10 gives degrees C per decade 
# Note that the same issue described above happens here. Thus, I calculate velocity with this package using the both the default (wrong) method in this package, and with the correct procedure (using mean raster).

start <- Sys.time()

trend_vocc = vocc::calcslope(rx = stack(tmean), divisor = 1)

### Wrong way
# spatial gradient
spgrad_out_vocc1 = vocc::spatialgrad(rx = stack(tmean), 
                                     y_dist = res(tmean),
                                     y_diff = NA) 
# Angle
angle_vocc1 = spgrad_out_vocc1$angle

# Calculating temperature velocity 
v_vocc1 = vocc::calcvelocity(grad = spgrad_out_vocc1, 
                             slope = trend_vocc)
v_vocc_rast1 <- rast(tmean[[1]])
v_vocc_rast1[] <- v_vocc1$velocity

vocc_time <- Sys.time() - start

### Right way
# spatial gradient
spgrad_out_vocc2 = vocc::spatialgrad(rx = stack(mean(tmean)), 
                                     y_dist = res(tmean),
                                     y_diff = NA) 
# Angle
angle_vocc2 = spgrad_out_vocc2$angle

# Calculating temperature velocity 
v_vocc2 = vocc::calcvelocity(grad = spgrad_out_vocc2, 
                             slope = trend_vocc)
v_vocc_rast2 <- rast(tmean[[2]])
v_vocc_rast2[] <- v_vocc2$velocity


##############
# VoCC

start <- Sys.time()

# calculate the trend
trend_VoCC = VoCC::tempTrend(r = stack(tmean),
                             th = 0.25*nlyr(tmean) ## set minimum # obs. to 1/4 time series length
)

# calculate the spatial gradient
spgrad_VoCC = VoCC::spatGrad(r = stack(tmean), 
                             projected = TRUE) 

# Calculating temperature velocity 
v_VoCC = VoCC::gVoCC(tempTrend = trend_VoCC, 
                     spatGrad = spgrad_VoCC)

# Angle
angle_VoCC = values(spgrad_VoCC[[2]])

VoCC_time <- Sys.time() - start

##############
# Bioshifts

start <- Sys.time()

# trend
trend_bio <- temp_grad(tmean, th = 0.25*nlyr(tmean))

# spatial gradient
spgrad_bio <- spatial_grad(tmean)

# velocity
v_bio <- gVelocity(spgrad_bio,trend_bio)

bio_time <- Sys.time() - start

# # Extra
# # To calculate the velocity N, we have to use the gradient North provided in the output
# # According to the calculations below, we can use the gradient North directly. It is not necessary to velocity rates to the north direction using angles 
# 
# # get data
# v_data <- as.data.frame(spgrad_bio)
# v_data <- na.omit(v_data)
# head(v_data)
# 
# # Convert gradients to the North direction using the angle
# temperature_change_rate <- v_data$Grad  # °C/km
# angle_degrees <- v_data$angle-90        # degrees
# 
# # Convert angle to radians
# angle_radians <- angle_degrees * pi / 180
# 
# # Calculate temperature change in North and West directions
# temperature_change_north <- temperature_change_rate * sin(angle_radians)
# temperature_change_west <- temperature_change_rate * cos(angle_radians)
# 
# plot(v_data$NS,temperature_change_north)
# plot(v_data$WE,temperature_change_west)
# 
# # See! Use NS directly

########################################################
## Compare velocities maps
# climetrics
tmp <- rast(c(climetrics1=v_climetrics1,
              climetrics2=v_climetrics2,
              vocc1=v_vocc_rast1,
              vocc2=v_vocc_rast2,
              VoCC=rast(v_VoCC[[1]]),
              v_bio=v_bio))
# very different values
tmp
plot(tmp)


bioshifts_time <- tictoc::toc()


##############
## Compare velocities correlations
velocities <- data.frame(tmp)
head(velocities)

# climetrics vs VoCC
ggpairs(velocities)

# as we can see, velocity values for climatrics1 and vocc1 are identical (but not perfect because climatrics1 squeezes outliers), but both are different then VoCC1 and VoCC2 because they handle the spatial gradient differently. When we correct the velocities from climatrics and vocc (climatrics2 and vocc2, respectively), the values are identical to the ones in VoCC.

## Compare time duration for calculating velocities
timesFuncs <- data.frame(t(data.frame(climetrics_time[[1]],vocc_time[[1]],VoCC_time[[1]],bio_time[[1]])))
names(timesFuncs) <- "elapsed"
timesFuncs$Package <- gsub("_time..1..","",rownames(timesFuncs))

ggplot(timesFuncs, aes(x = Package, y = elapsed))+
    geom_col() + ylab("Time elapsed (seconds)")

