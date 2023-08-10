# Setup
# remotes::install_github("hughjonesd/ggmagnify")

rm(list=ls())
gc()

library(dplyr)
library(ggplot2)
library(tune)
library(GGally)
library(ggmagnify)
library(lme4)

# computer = "muse"
computer = "personal"

########################
if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    
    work_dir <- getwd()
    
    output_dir <- here::here(work_dir,"Data/SHIFT")
}
if(computer == "personal"){
    
    work_dir <- getwd()
    
    output_dir <- here::here(work_dir,"Data/SHIFT")
}

########################
# load functions
source("R/my_functions.R")
# source settings
source("R/settings.R")

########################
# Load in all shifts SA ens

shifts <- list.files(here::here(output_dir,"Mar"), 
                     pattern = " ens shift SA.csv", 
                     recursive = TRUE, full.names = TRUE)

shifts <- lapply(shifts, function(x) {
    try({read.csv(x)})
})

shifts <- data.table::rbindlist(shifts)
shifts_mar <- data.frame(shifts)
shifts_mar$ECO = "Mar"

shifts <- list.files(here::here(output_dir,"Ter"), 
                     pattern = " ens shift SA.csv", 
                     recursive = TRUE, full.names = TRUE)

shifts <- lapply(shifts, function(x){
    try({read.csv(x)})
})
shifts <- data.table::rbindlist(shifts)
shifts_ter <- data.frame(shifts)
shifts_ter$ECO = "Ter"

shifts <- rbind(shifts_mar,shifts_ter)

# convert from degrees to km
shifts[,grep("shift",names(shifts))] <- shifts[,grep("shift",names(shifts))]*111.1 

# convert from km to km/year
shifts$DUR <- round(shifts$END,1)-round(shifts$START,1)
shifts[,grep("shift",names(shifts))] <- shifts[,grep("shift",names(shifts))]/shifts$DUR 

shifts$new_ID <- paste(shifts$ID,
                       shifts$Species,
                       round(shifts$START,1),
                       round(shifts$END,1),
                       shifts$Type)

dim(shifts)
any(duplicated(shifts$new_ID))

#####################
# get shifts at the range positions

# Predicted 1
# At the LE, use quant 1%
# At the TE, use quant 99%
# At the O, use quant 50%

# Predicted 2
# At the LE, quant 5%
# At the TE, quant 95%
# At the O, quant 50%

# Predicted 3
# At the LE, quant 10%
# At the TE, quant 90%
# At the O, quant 50%

shifts_pred1 <- pbapply::pblapply(
    
    shifts$new_ID, function(x){
        
        # subset shift i
        tmp <- shifts %>% filter(new_ID == x)
        
        # get shifts at range positions
        tmp2 <- data.frame(Param = c("LE","O","TE"),
                           SHIFT_sdm_1 = c(mean(tmp$nsQuantVelocity_quant0p01),
                                           mean(tmp$nsCentroidVelocity),
                                           mean(tmp$nsQuantVelocity_quant0p99)),
                           SHIFT_sdm_2 = c(mean(tmp$nsQuantVelocity_quant0p05),
                                           mean(tmp$nsCentroidVelocity),
                                           mean(tmp$nsQuantVelocity_quant0p95)),
                           SHIFT_sdm_3 = c(mean(tmp$nsQuantVelocity_quant0p1),
                                           mean(tmp$nsCentroidVelocity),
                                           mean(tmp$nsQuantVelocity_quant0p9)))
        
        tmp <- data.frame(tmp2,
                          tmp[,c(1:3,21:23,28)])
        
    })
shifts_pred1 <- data.table::rbindlist(shifts_pred1)
shifts_pred1 <- data.frame(shifts_pred1)

shifts_pred1$new_ID <- paste(shifts_pred1$new_ID, shifts_pred1$Param)
shifts_pred1 <- na.omit(shifts_pred1)
dim(shifts_pred1)
length(unique(shifts_pred1$Species))

# Predicted 2
# At the LE, quant 5%
# At the TE, quant 95%
# At the O, quant 50%
shifts_pred2 <- pbapply::pblapply(
    shifts$new_ID, function(x){
        
        # subset shift i
        tmp <- shifts %>% filter(new_ID == x)
        
        # get shifts at range positions
        LE <- tmp$shift_5.
        O <- tmp$shift_50.
        TE <- tmp$shift_95.
        
        # get latitudinal values at start shift
        Lat_start_LE <- tmp$w_t1_5.
        Lat_start_O <- tmp$w_t1_50.
        Lat_start_TE <- tmp$w_t1_95.
        
        # get latitudinal values at end shift
        Lat_end_LE <- tmp$w_t2_5.
        Lat_end_O <- tmp$w_t2_50.
        Lat_end_TE <- tmp$w_t2_95.
        
        tmp2 <- data.frame(Param = c("LE","O","TE"),
                           SHIFT_sdm_2 = c(LE,O,TE),
                           Lat_sdm_start_2 = c(Lat_start_LE,Lat_start_O,Lat_start_TE),
                           Lat_sdm_end_2 = c(Lat_end_LE,Lat_end_O,Lat_end_TE))
        
        tmp <- data.frame(tmp2,
                          tmp[,25:33])
    })
shifts_pred2 <- data.table::rbindlist(shifts_pred2)
shifts_pred2 <- data.frame(shifts_pred2)

shifts_pred2$new_ID <- paste(shifts_pred2$new_ID, shifts_pred2$Param)
shifts_pred2 <- na.omit(shifts_pred2)
dim(shifts_pred2)
length(unique(shifts_pred2$Species))


# Predicted 2
# At the LE, quant 5%
# At the TE, quant 95%
# At the O, quant 50%
shifts_pred3 <- pbapply::pblapply(
    shifts$new_ID, function(x){
        
        # subset shift i
        tmp <- shifts %>% filter(new_ID == x)
        
        # get shifts at range positions
        LE <- tmp$shift_10.
        O <- tmp$shift_50.
        TE <- tmp$shift_90.
        
        # get latitudinal values at start shift
        Lat_start_LE <- tmp$w_t1_10.
        Lat_start_O <- tmp$w_t1_50.
        Lat_start_TE <- tmp$w_t1_90.
        
        # get latitudinal values at end shift
        Lat_end_LE <- tmp$w_t2_10.
        Lat_end_O <- tmp$w_t2_50.
        Lat_end_TE <- tmp$w_t2_90.
        
        tmp2 <- data.frame(Param = c("LE","O","TE"),
                           SHIFT_sdm_3 = c(LE,O,TE),
                           Lat_sdm_start_3 = c(Lat_start_LE,Lat_start_O,Lat_start_TE),
                           Lat_sdm_end_3 = c(Lat_end_LE,Lat_end_O,Lat_end_TE))
        
        tmp <- data.frame(tmp2,
                          tmp[,25:33])
    })
shifts_pred3 <- data.table::rbindlist(shifts_pred3)
shifts_pred3 <- data.frame(shifts_pred3)

shifts_pred3$new_ID <- paste(shifts_pred3$new_ID, shifts_pred3$Param)
shifts_pred3 <- na.omit(shifts_pred3)
dim(shifts_pred3)
length(unique(shifts_pred3$Species))

any(duplicated(shifts_pred3$new_ID))

any(is.na(shifts_pred3$SHIFT_sdm))


### Group all
all.equal(shifts_pred1$new_ID,shifts_pred2$new_ID,shifts_pred3$new_ID)
shifts_pred <- cbind(shifts_pred1[,1:4],shifts_pred2[,2:4],shifts_pred3[,2:13])
head(shifts_pred)
dim(shifts_pred)

table(shifts_pred$ECO, shifts_pred$Param)

#### Correlation
ggpairs(shifts_pred[-which(shifts_pred$Param=="O"),], columns = c(2,5,8),
        ggplot2::aes(colour=Param))

########################
# Merge with bioshifts
biov1 <- read.csv("Data/Bioshifts/biov1_fixednames.csv", header = T)

biov1 <- biov1 %>% filter(Type == "LAT")
biov1$new_ID <- paste(biov1$ID,
                      biov1$sp_name_std_v1,
                      round(biov1$START,1),
                      round(biov1$END,1),
                      biov1$Type,
                      biov1$Param)

any(duplicated(biov1$new_ID))
# test = which(duplicated(biov1$new_ID))
# biov1$new_ID[test][1]
# biov1[which(biov1$new_ID == biov1$new_ID[test][1]),]

table(biov1$ECO, biov1$Param)

biov1_pred <- merge(biov1, shifts_pred[,c(2:10,19)], 
                    by = "new_ID")
str(biov1_pred)

table(biov1_pred$ECO, biov1_pred$Param)
table(biov1_pred$ECO, biov1_pred$Param, biov1_pred$Class)

########################
# Plot

# Marines


# Optimum
M_O <- ggplot(
    biov1_pred %>%
        filter(ECO == "M" & Param == "O"), 
    aes(x = SHIFT, y = SHIFT_sdm_1, color = Class))+
    ggtitle("Centroid")+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

# predicted 1 = 1% - 99%
M_LE_1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "M" & Param == "LE"), 
    aes(x = SHIFT, y = SHIFT_sdm_1, color = Class))+
    ggtitle("Leading edge")+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

M_TE_1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "M" & Param == "TE"), 
    aes(x = SHIFT, y = SHIFT_sdm_1, color = Class))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()



# predicted 2 = 5% - 95%
M_LE_2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "M" & Param == "LE"), 
    aes(x = SHIFT, y = SHIFT_sdm_2, color = Class))+
    ggtitle("Leading edge")+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

M_TE_2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "M" & Param == "TE"), 
    aes(x = SHIFT, y = SHIFT_sdm_2, color = Class))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()



# predicted 3 = 10% - 90%
M_LE_3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "M" & Param == "LE"), 
    aes(x = SHIFT, y = SHIFT_sdm_3, color = Class))+
    ggtitle("Leading edge")+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

M_TE_3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "M" & Param == "TE"), 
    aes(x = SHIFT, y = SHIFT_sdm_3, color = Class))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

M_O
gridExtra::grid.arrange(M_TE_1,M_LE_1,ncol=2,top="TE = 1%, LE = 99%")
gridExtra::grid.arrange(M_TE_2,M_LE_2,ncol=2,top="TE = 5%, LE = 95%")
gridExtra::grid.arrange(M_TE_3,M_LE_3,ncol=2,top="TE = 10%, LE = 90%")

# just fish
table(biov1_pred$Class,biov1_pred$ECO)

ggplot(
    biov1_pred %>%
        filter(ECO == "M" & !Param == "TE" & Class == "Actinopterygii"), 
    aes(x = SHIFT, y = SHIFT_sdm_2, color = Param, size = DUR))+
    ggtitle("Actinopterygii")+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()+
    facet_wrap(.~Param)

# just Malacostraca
ggplot(
    biov1_pred %>%
        filter(ECO == "M" & !Param == "TE" & Class == "Malacostraca"), 
    aes(x = SHIFT, y = SHIFT_sdm_2, color = Param, size = DUR))+
    ggtitle("Malacostraca")+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()+
    facet_wrap(.~Param)


# Terrestrials

# Optimum
M_O <- ggplot(
    biov1_pred %>%
        filter(ECO == "T" & Param == "O"), 
    aes(x = SHIFT, y = SHIFT_sdm_1, color = Class))+
    ggtitle("Centroid")+
    geom_point()+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

# predicted 1 = 1% - 99%
M_LE_1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "T" & Param == "LE"), 
    aes(x = SHIFT, y = SHIFT_sdm_1, color = Class))+
    ggtitle("Leading edge")+
    geom_point()+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

M_TE_1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "T" & Param == "TE"), 
    aes(x = SHIFT, y = SHIFT_sdm_1, color = Class))+
    ggtitle("Trailing edge")+
    geom_point()+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()



# predicted 2 = 5% - 95%
M_LE_2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "T" & Param == "LE"), 
    aes(x = SHIFT, y = SHIFT_sdm_2, color = Class))+
    ggtitle("Leading edge")+
    geom_point()+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

M_TE_2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "T" & Param == "TE"), 
    aes(x = SHIFT, y = SHIFT_sdm_2, color = Class))+
    ggtitle("Trailing edge")+
    geom_point()+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()



# predicted 3 = 10% - 90%
M_LE_3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "T" & Param == "LE"), 
    aes(x = SHIFT, y = SHIFT_sdm_3, color = Class))+
    ggtitle("Leading edge")+
    geom_point()+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

M_TE_3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "T" & Param == "TE"), 
    aes(x = SHIFT, y = SHIFT_sdm_3, color = Class))+
    ggtitle("Trailing edge")+
    geom_point()+
    geom_point(alpha = .5)+
    ylab("Potential shift (km/yr)")+
    xlab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

M_O
# x0, y0, x1, y1 of target:
from <- c(-20, -20, 20, 20)
# x0, y0, x1, y1 of inset:
to <- c(30, 30, 140,140)
M_O + geom_magnify(from = from, to = to)

gridExtra::grid.arrange(M_TE_1,M_LE_1,ncol=2,top="TE = 1%, LE = 99%")
gridExtra::grid.arrange(M_TE_2,M_LE_2,ncol=2,top="TE = 5%, LE = 95%")
gridExtra::grid.arrange(M_TE_3,M_LE_3,ncol=2,top="TE = 10%, LE = 90%")

################
## Predicted shifts vs Temperature velocities

## Does bioclimatic velocity explains shifts better then temperature velocity?
model_Mar_v <- lmer(SHIFT~v.lat.mean*Param+
                        (1|Article_ID),
                    data = biov1_pred %>%
                        filter(ECO == "M"))
summary(model_Mar_v)
MuMIn::r.squaredGLMM(model_Mar_v)

model_Mar_b <- lmer(SHIFT~SHIFT_sdm_2*Param+
                        (1|Article_ID),
                    data = biov1_pred %>%
                        filter(ECO == "M"))
summary(model_Mar_b)
MuMIn::r.squaredGLMM(model_Mar_b)

model_Ter_v <- lmer(SHIFT~v.lat.mean*Param+
                        (1|Article_ID),
                    data = biov1_pred %>%
                        filter(ECO == "T"))
summary(model_Ter_v)
MuMIn::r.squaredGLMM(model_Ter_v)

model_Ter_b <- lmer(SHIFT~SHIFT_sdm_2*Param+
                        (1|Article_ID),
                    data = biov1_pred %>%
                        filter(ECO == "M"))
summary(model_Ter_b)
MuMIn::r.squaredGLMM(model_Ter_b)


## Marine
p1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "M"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_1, color = Param))+
    ggtitle("TE = 1%, LE = 99%")+
    geom_point()+
    geom_point(alpha = .5)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Potential shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

p2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "M"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_2, color = Param))+
    ggtitle("TE = 5%, LE = 95%")+
    geom_point()+
    geom_point(alpha = .5)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Potential shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

p3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "M"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_3, color = Param))+
    ggtitle("TE = 10%, LE = 90%")+
    geom_point()+
    geom_point(alpha = .5)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Potential shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

gridExtra::grid.arrange(p1,p2,p3,ncol=3,top="Marine")



## Terrestrial
p1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "T"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_1, color = Param))+
    ggtitle("TE = 1%, LE = 99%")+
    geom_point()+
    geom_point(alpha = .5)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Potential shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

p2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "T"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_2, color = Param))+
    ggtitle("TE = 5%, LE = 95%")+
    geom_point()+
    geom_point(alpha = .5)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Potential shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

p3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "T"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_3, color = Param))+
    ggtitle("TE = 10%, LE = 90%")+
    geom_point()+
    geom_point(alpha = .5)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Potential shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

gridExtra::grid.arrange(p1,p2,p3,ncol=3,top="Terrestrials")


################
## Difference >> Predicted shifts vs Observed shift
## Correlate with shift duration



## Difference >> Predicted shifts vs Observed shift
## Correlate with shift start



## Predicted shifts vs Observed shift
## Subset TSS > 0.7


