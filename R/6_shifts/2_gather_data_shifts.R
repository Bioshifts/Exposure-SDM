# Setup
# remotes::install_github("hughjonesd/ggmagnify")
# devtools::install_github(c("yihui/servr", "hafen/rmote")) 

# rmote::Start_rmote() 

rm(list=ls())
gc()

list.of.packages <- c("dplyr","ggplot2","tune","GGally","lme4","terra", "rgdal","maps","readr","tidyr","glmmTMB")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(list.of.packages, require, character.only = TRUE)


###
# set computer
# computer = "matrics"
computer = "personal"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}


source("R/settings.R")
source("R/my_functions.R")
source(here::here("R/fetchGHdata.R"))

work_dir <- getwd()
shift_dir <- here::here(work_dir,"Data/SHIFT")
sdm_dir <- here::here(work_dir,"Data/SDMs")
velocity_SA_dir <- here::here(work_dir,"Data/Velocity_SA")

###
### Load in Bioshifts ----
bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_sdms_selection(bioshifts)
bioshifts$sp_name_std <- gsub(" ","_",bioshifts$sp_name_std)
# N species
length(unique(bioshifts$sp_name_std))


# New ID
bioshifts$new_ID <- paste(bioshifts$ID,
                          bioshifts$sp_name_std,
                          round(bioshifts$Start,0),
                          round(bioshifts$End,0),
                          bioshifts$Type)

# any(duplicated(bioshifts$new_ID))
# test = which(duplicated(bioshifts$new_ID))
# bioshifts$new_ID[test][1]
# bioshifts[which(bioshifts$new_ID == bioshifts$new_ID[test][1]),]

bioshifts$Eco <- ifelse(bioshifts$Eco=="M", "Marine", "Terrestrial")

table(bioshifts$Eco, bioshifts$Param)
dim(bioshifts)


# Plug in new velocity metrics
velocities <- list.files(velocity_SA_dir, full.names = TRUE, pattern = ".csv")
velocities <- lapply(velocities, read.csv)
velocities <- data.table::rbindlist(velocities,fill = TRUE)
velocities <- velocities[,c("ID","vel_mat","v.lat.sd.mat","v.lat.median.sst","v.lat.sd.sst")]
velocities$vel_mat[which(is.na(velocities$vel_mat))] <- velocities$v.lat.median.sst[which(is.na(velocities$vel_mat))]
velocities$vel_mat <- velocities$vel_mat
velocities$vel_mat[which(is.na(velocities$vel_mat))] <- velocities$v.lat.median.sst[which(is.na(velocities$vel_mat))]
velocities$vel_mat_sd <- velocities$v.lat.sd.mat
velocities$vel_mat_sd[which(is.na(velocities$vel_mat_sd))] <- velocities$v.lat.sd.sst[which(is.na(velocities$vel_mat_sd))]

bioshifts <- merge(bioshifts,
                   velocities[,c("ID","vel_mat","vel_mat_sd")],
                   by = "ID",
                   all.x = TRUE)

###
### Load occ data ----
# add group
N_OCC <- read_csv("Data/N_OCC.csv")
N_OCC$scientificName <- gsub(" ","_",N_OCC$scientificName)

###
### Load in SDM CV ----
sdms_CV <- read.csv(here::here("Data","sdms_CV.csv"))
# sdms_CV <- sdms_CV %>% dplyr::filter(metric.eval == "TSS")
# sdms_CV <- sdms_CV %>% dplyr::filter(validation > 0)
sdms_CV$Eco <- sdms_CV$ECO

sdms_CV$Species <- sapply(sdms_CV$full.name, function(x){
    tmp = strsplit(x,"_")[[1]][1]
    gsub("[.]","_",tmp)
})
sdms_CV$Eco <- ifelse(sdms_CV$Eco=="Mar", "Marine", "Terrestrial")

all(sdms_CV$Species %in% N_OCC$scientificName)

# N species
length(unique(sdms_CV$Species))

# N Terrestrial
length(unique(sdms_CV$Species[which(sdms_CV$Eco=="Terrestrial")]))
# N Marine
length(unique(sdms_CV$Species[which(sdms_CV$Eco=="Marine")]))

# add Taxonomic group and N occ
sdms_CV <- merge(sdms_CV, N_OCC[,c("scientificName","class","N_OCC")], 
                 by.x = "Species", 
                 by.y = "scientificName",
                 all.x = TRUE)

par(mfrow=c(2,2))
plot(log(sdms_CV$N_OCC),sdms_CV$validation)
plot(log(sdms_CV$N_OCC),sdms_CV$calibration)
plot(sdms_CV$calibration,sdms_CV$validation)
plot(log(sdms_CV$N_OCC),sdms_CV$sensitivity)
# dev.off()

sdms_CV_plot <- sdms_CV[,c("algo","calibration","validation","class","Eco","metric.eval")] %>%
    tidyr::gather("Type", "value", -c(algo,class,Eco,metric.eval)) 

# Terrestrial
## TSS
p1 <- ggplot(sdms_CV_plot %>%
                 dplyr::filter(Eco=="Terrestrial",
                               metric.eval=="TSS"), 
             aes(x = value, y = class, color = algo))+
    ggtitle("TSS")+
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = c(0,1))+
    facet_grid(Eco~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

## ROC
p2 <- ggplot(sdms_CV_plot %>%
                 dplyr::filter(Eco=="Terrestrial",
                               metric.eval=="ROC"), 
             aes(x = value, y = class, color = algo))+
    ggtitle("ROC")+
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = c(0,1))+
    facet_grid(Eco~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

gridExtra::grid.arrange(p1,p2,ncol=1)

# Marine
## TSS
p1 <- ggplot(sdms_CV_plot %>%
                 dplyr::filter(Eco=="Marine",
                               metric.eval=="TSS"), 
             aes(x = value, y = class, color = algo))+
    ggtitle("TSS")+
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = c(0,1))+
    facet_grid(Eco~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

## ROC
p2 <- ggplot(sdms_CV_plot %>%
                 dplyr::filter(Eco=="Marine",
                               metric.eval=="ROC"), 
             aes(x = value, y = class, color = algo))+
    ggtitle("ROC")+
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = c(0,1))+
    facet_grid(Eco~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

gridExtra::grid.arrange(p1,p2,ncol=1)

# total Terrestrial species
length(unique(sdms_CV$Species[which(sdms_CV$Eco=="Terrestrial")]))
# Proportion Terrestrials species MAXNET TSS > 0.5
sdms_CV %>%
    dplyr::filter(Eco=="Terrestrial",
                  algo=="MAXNET",
                  metric.eval=="TSS") %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(TSS = mean(calibration)) %>%
    dplyr::summarise(length(which(TSS > .5))/length(TSS))

# Proportion Terrestrials species MAXNET TSS > 0.7
sdms_CV %>%
    dplyr::filter(Eco=="Terrestrial",
                  algo=="MAXNET",
                  metric.eval=="TSS") %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(TSS = mean(calibration)) %>%
    dplyr::summarise(length(which(TSS > .7))/length(TSS))

# total Terrestrial species
length(unique(sdms_CV$Species[which(sdms_CV$Eco=="Terrestrial")]))

# Proportion Marine species MAXNET TSS > 0.5
sdms_CV %>%
    dplyr::filter(Eco=="Marine",
                  algo=="MAXNET",
                  metric.eval=="TSS") %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(TSS = mean(calibration)) %>%
    dplyr::summarise(length(which(TSS > .5))/length(TSS))

# Proportion Marine species MAXNET TSS > 0.7
sdms_CV %>%
    dplyr::filter(Eco=="Marine",
                  algo=="MAXNET",
                  metric.eval=="TSS") %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(TSS = mean(calibration)) %>%
    dplyr::summarise(length(which(TSS > .7))/length(TSS))


###
### Load in all shifts SA ens ----

shifts_ens <- list.files(here::here(shift_dir,"Mar"), 
                         pattern = "_SA.csv", 
                         recursive = TRUE, full.names = TRUE)
length(shifts_ens)

shifts_ens <- lapply(shifts_ens, function(x) {
    read.csv(x)
})

shifts_ens <- data.table::rbindlist(shifts_ens, fill = TRUE)

shifts_ens_mar <- data.frame(shifts_ens)
shifts_ens_mar$Eco = "Marine"

length(unique(shifts_ens_mar$Species))

shifts_ens <- list.files(here::here(shift_dir,"Ter"), 
                         pattern = "_SA.csv", 
                         recursive = TRUE, full.names = TRUE)

shifts_ens <- lapply(shifts_ens, function(x) {
    tmp <- read.csv(x)
    tmp$name <- x
    return(tmp)
})


shifts_ens <- data.table::rbindlist(shifts_ens,fill = TRUE)
shifts_ens_ter <- data.frame(shifts_ens)
shifts_ens_ter$Eco = "Terrestrial"

shifts_ens$Start <- shifts_ens$START
shifts_ens$End <- shifts_ens$END

length(unique(shifts_ens_ter$Species))

matching_columns <- match(names(shifts_ens_mar),names(shifts_ens_ter))

shifts_ens <- rbind(shifts_ens_mar,
                    shifts_ens_ter[,names(shifts_ens_mar)])

length(unique(shifts_ens$Species))
length(unique(shifts_ens$Species[which(shifts_ens$Eco=="Terrestrial")]))
length(unique(shifts_ens$Species[which(shifts_ens$Eco=="Marine")]))
# N SAs
length(unique(shifts_ens$ID))

# N species
length(unique(shifts_ens$Species))
# N shifts
length(unique(shifts_ens$ID))
# Eco
table(shifts_ens$Eco)

### Correlation shifts % ----

## Terrestrial
ggpairs(shifts_ens %>%
            filter(Eco == "Terrestrial"), 
        columns = c(23:26, 9, 27:29))

ggpairs(shifts_ens %>%
            filter(Eco == "Marine"), 
        columns = c(23:26, 9, 27:29))


###
### Merge with bioshifts ----
#### plug in shift sdms
shifts_ens$new_ID <- paste(shifts_ens$ID, shifts_ens$Species, round(as.numeric(shifts_ens$START),0), round(as.numeric(shifts_ens$END),0), shifts_ens$Type)

bioshifts_pred <- merge(shifts_ens,
                        data.frame(bioshifts[,-which(names(bioshifts) %in% names(shifts_ens))],
                                   new_ID=bioshifts$new_ID),
                        by = c("new_ID"),
                        all.x = TRUE)

### Plug in TSS ----
avg_TSS <- sdms_CV %>%
    dplyr::filter(metric.eval=="TSS") %>%
    group_by(Species) %>%
    summarise(TSS_vali = median(validation),
              TSS_cali = median(calibration))

all(bioshifts_pred$sp_name_std %in% avg_TSS$Species)

bioshifts_pred <- merge(
    bioshifts_pred, 
    data.frame(avg_TSS[,-which(names(avg_TSS) %in% names(bioshifts_pred))],
               Species=avg_TSS$Species), 
    by.x = "sp_name_std",
    by.y = "Species",
    all.x = TRUE)

###
# Load velocity at edges ----
vel_edge <- read.csv(here::here("Data/vel_edge.csv"))
# N SAs
length(unique(sapply(vel_edge$ID, function(x) strsplit(x, " ")[[1]][1])))
vel_edge$new_ID <- vel_edge$ID
vel_edge <- vel_edge %>% select(-ID)

bioshifts_pred <- merge(bioshifts_pred, 
                        data.frame(vel_edge[,-which(names(vel_edge) %in% names(bioshifts_pred))],
                                   new_ID=vel_edge$new_ID), 
                        by = "new_ID",
                        all.x = TRUE)


### select the observed shift variable  
names(bioshifts_pred)

# N species
length(unique(bioshifts_pred$sp_name_std))

bioshifts_pred %>%
    group_by(Eco) %>%
    summarise(N=length(unique(sp_name_std)))

###
### Compare velocities across edges ----
ggplot(bioshifts_pred %>%
           mutate(TE=bv.lat.1,O=bv.lat.median,LE=bv.lat.9) %>%
           select(TE,O,LE) %>%
           gather(),
       aes(x=value,fill=key))+
    geom_histogram()+
    facet_wrap(.~key,ncol = 2,scales = "free")

ggplot(bioshifts_pred %>%
           mutate(TE=edge1,O=edge5,LE=edge9) %>%
           select(TE,O,LE) %>%
           gather(),
       aes(x=value,fill=key))+
    geom_histogram()+
    facet_wrap(.~key,ncol = 2,scales = "free")

###
### Shift predicted vs observed ----

#### Centroid
ggplot(bioshifts_pred %>%
           filter(Param == "O",
                  # (bv.lat.median > 0 & ShiftRate > 0) | (bv.lat.median < 0 & ShiftRate < 0)
           ), 
       aes(x = bv.lat.median, y = ShiftRate, color = Eco, size = Duration))+
    ggtitle("Centroid shifts")+
    geom_point(alpha = .3)+
    # scale_x_log10()+scale_y_log10()+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    geom_smooth(method = "lm")+
    facet_wrap(.~Eco,ncol = 2,scales = "free")

ggplot(bioshifts_pred %>%
           filter(Param == "O"), 
       aes(x = edge5, y = ShiftRate, color = Eco, size = Duration))+
    ggtitle("Centroid shifts")+
    geom_point(alpha = .3)+
    xlab("Climate velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    geom_smooth(method = "lm")+
    facet_wrap(.~Eco,ncol = 2,scales = "free")

ggplot(bioshifts_pred %>%
           filter(Param == "O"), 
       aes(x = edge5, y = bv.lat.median, color = Eco, size = Duration))+
    ggtitle("Centroid shifts")+
    geom_point(alpha = .3)+
    xlab("Climate velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    geom_smooth(method = "lm")+
    facet_wrap(.~Eco,ncol = 2,scales = "free")


#### Leading edge
ggplot(bioshifts_pred %>%
           filter(Param == "LE"), 
       aes(x = bv.lat.9, y = ShiftRate, color = Eco, size = Duration))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    geom_smooth(method = "lm")+
    facet_wrap(.~Eco,ncol = 2,scales = "free")

ggplot(bioshifts_pred %>%
           filter(Param == "LE"), 
       aes(x = edge9, y = ShiftRate, color = Eco, size = Duration))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Climate velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    geom_smooth(method = "lm")+
    facet_wrap(.~Eco,ncol = 2,scales = "free")

ggplot(bioshifts_pred %>%
           filter(Param == "LE"), 
       aes(x = edge9, y = bv.lat.9, color = Eco, size = Duration))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Climate velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    geom_smooth(method = "lm")+
    facet_wrap(.~Eco, ncol = 2, scales = "free")

#### Trailing edge
ggplot(bioshifts_pred %>%
           filter(Param == "TE"), 
       aes(x = bv.lat.1, y = ShiftRate, color = Eco, size = Duration))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    geom_smooth(method = "lm")+
    facet_wrap(.~Eco,ncol = 2,scales = "free")

ggplot(bioshifts_pred %>%
           filter(Param == "TE"), 
       aes(x = edge1, y = ShiftRate, color = Eco, size = Duration))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Climate velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    geom_smooth(method = "lm")+
    facet_wrap(.~Eco,ncol = 2,scales = "free")

ggplot(bioshifts_pred %>%
           filter(Param == "TE"), 
       aes(x = edge1, y = bv.lat.1, color = Eco, size = Duration))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Climate velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    geom_smooth(method = "lm")+
    facet_wrap(.~Eco,ncol = 2,scales = "free")

###
### Models ----

CE_model <- glmmTMB(ShiftRate ~ bv.lat.median + edge5 + Type + 
                        (1 | Eco), 
                    weights = sensitivity,
                    data = bioshifts_pred %>% filter(Param == "O"))
summary(gau1_allW)

tmp <- bioshifts_pred %>% 
    filter(Param == "O") 


