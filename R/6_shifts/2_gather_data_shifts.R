# Setup
# remotes::install_github("hughjonesd/ggmagnify")
# devtools::install_github(c("yihui/servr", "hafen/rmote")) 

# rmote::start_rmote() 

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
source(here::here("R/fetchGHdata.R"))

work_dir <- getwd()
shift_dir <- here::here(work_dir,"Data/SHIFT")
sdm_dir <- here::here(work_dir,"Data/SDMs")
velocity_SA_dir <- here::here(work_dir,"Data/Velocity_SA")

###
### Load in Bioshifts ----
biov1 <- read.csv(here::here(work_dir,"Data","Bioshifts",Bioshifts_DB_v1), header = T)
# select columns of interest
biov1 <- biov1 %>%
    select(Article_ID, Study_ID, ID, Type, Param, START, END, DUR, ECO, UNIT, 
           SHIFT, Class, group, sp_name_std_v1, v.lat.mean, v.ele.mean)

# Lat only
biov1 <- biov1 %>% 
    mutate(
        Type = case_when(
            Type=="HOR" ~ "LAT",
            TRUE ~ as.character(Type)),
        velocity = ifelse(Type == "LAT", v.lat.mean, v.ele.mean),
        SHIFT = ifelse(UNIT == "m/year", SHIFT/1000, SHIFT) ) %>%
    filter(Type == "LAT") %>%
    select(!v.lat.mean,!v.ele.mean,!UNIT)

# N species
length(unique(biov1$sp_name_std_v1))

# N species
length(unique(biov1$sp_name_std_v1))



# New ID
biov1$new_ID <- paste(biov1$ID,
                      biov1$sp_name_std_v1,
                      round(biov1$START,0),
                      round(biov1$END,0),
                      biov1$Type)

# any(duplicated(biov1$new_ID))
# test = which(duplicated(biov1$new_ID))
# biov1$new_ID[test][1]
# biov1[which(biov1$new_ID == biov1$new_ID[test][1]),]

biov1$ECO <- ifelse(biov1$ECO=="M", "Marine", "Terrestrial")

table(biov1$ECO, biov1$Param)
dim(biov1)


# Plug in new velocity metrics
velocities <- list.files(velocity_SA_dir, full.names = TRUE, pattern = ".csv")
velocities <- lapply(velocities, read.csv)
velocities <- data.table::rbindlist(velocities,fill = TRUE)
velocities <- velocities[,c("ID","v.lat.median.mat","v.lat.sd.mat","v.lat.median.sst","v.lat.sd.sst")]
velocities$v.lat.median.mat[which(is.na(velocities$v.lat.median.mat))] <- velocities$v.lat.median.sst[which(is.na(velocities$v.lat.median.mat))]
velocities$vel_mat <- velocities$v.lat.median.mat
velocities$vel_mat[which(is.na(velocities$vel_mat))] <- velocities$v.lat.median.sst[which(is.na(velocities$vel_mat))]
velocities$vel_mat_sd <- velocities$v.lat.sd.mat
velocities$vel_mat_sd[which(is.na(velocities$vel_mat_sd))] <- velocities$v.lat.sd.sst[which(is.na(velocities$vel_mat_sd))]

biov1 <- merge(biov1,
               velocities[,c("ID","vel_mat","vel_mat_sd")],
               by = "ID",
               all.x = TRUE)

# Update velocity ----
biov1$vel_mat[which(is.na(biov1$vel_mat))] <- biov1$v.lat.mean[which(is.na(biov1$vel_mat))]

###
### Load occ data ----
# add group
n_occ <- read_csv("Data/n_occ2.csv")
n_occ$scientificName <- gsub(" ","_",n_occ$scientificName)

###
### Load in SDM CV ----
sdms_CV <- read.csv(here::here("Data","sdms_CV.csv"))
# sdms_CV <- sdms_CV %>% dplyr::filter(metric.eval == "TSS")
# sdms_CV <- sdms_CV %>% dplyr::filter(validation > 0)

sdms_CV$Species <- sapply(sdms_CV$full.name, function(x){
    tmp = strsplit(x,"_")[[1]][1]
    gsub("[.]","_",tmp)
})
sdms_CV$ECO <- ifelse(sdms_CV$ECO=="Mar", "Marine", "Terrestrial")

all(sdms_CV$Species %in% n_occ$scientificName)

# N species
length(unique(sdms_CV$Species))

# N Terrestrial
length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Terrestrial")]))
# N Marine
length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Marine")]))

# add Taxonomic group and N occ
sdms_CV <- merge(sdms_CV, n_occ[,c("scientificName","Group","n_occ")], 
                 by.x = "Species", 
                 by.y = "scientificName",
                 all.x = TRUE)

par(mfrow=c(2,2))
plot(log(sdms_CV$n_occ),sdms_CV$validation)
plot(log(sdms_CV$n_occ),sdms_CV$calibration)
plot(sdms_CV$calibration,sdms_CV$validation)
plot(log(sdms_CV$n_occ),sdms_CV$sensitivity)
# dev.off()

sdms_CV_plot <- sdms_CV[,c("algo","calibration","validation","Group","ECO","metric.eval")] %>%
    tidyr::gather("Type", "value", -c(algo,Group,ECO,metric.eval)) 

# Terrestrial
## TSS
p1 <- ggplot(sdms_CV_plot %>%
                 dplyr::filter(ECO=="Terrestrial",
                               metric.eval=="TSS"), 
             aes(x = value, y = Group, color = algo))+
    ggtitle("TSS")+
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = c(0,1))+
    facet_grid(ECO~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

## ROC
p2 <- ggplot(sdms_CV_plot %>%
                 dplyr::filter(ECO=="Terrestrial",
                               metric.eval=="ROC"), 
             aes(x = value, y = Group, color = algo))+
    ggtitle("ROC")+
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = c(0,1))+
    facet_grid(ECO~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

gridExtra::grid.arrange(p1,p2,ncol=1)

# Marine
## TSS
p1 <- ggplot(sdms_CV_plot %>%
                 dplyr::filter(ECO=="Marine",
                               metric.eval=="TSS"), 
             aes(x = value, y = Group, color = algo))+
    ggtitle("TSS")+
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = c(0,1))+
    facet_grid(ECO~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

## ROC
p2 <- ggplot(sdms_CV_plot %>%
                 dplyr::filter(ECO=="Marine",
                               metric.eval=="ROC"), 
             aes(x = value, y = Group, color = algo))+
    ggtitle("ROC")+
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = c(0,1))+
    facet_grid(ECO~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

gridExtra::grid.arrange(p1,p2,ncol=1)

# total Terrestrial species
length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Terrestrial")]))
# Proportion Terrestrials species MAXNET TSS > 0.5
sdms_CV %>%
    dplyr::filter(ECO=="Terrestrial",
                  algo=="MAXNET",
                  metric.eval=="TSS") %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(TSS = mean(calibration)) %>%
    dplyr::summarise(length(which(TSS > .5))/length(TSS))

# Proportion Terrestrials species MAXNET TSS > 0.7
sdms_CV %>%
    dplyr::filter(ECO=="Terrestrial",
                  algo=="MAXNET",
                  metric.eval=="TSS") %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(TSS = mean(calibration)) %>%
    dplyr::summarise(length(which(TSS > .7))/length(TSS))

# total Terrestrial species
length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Terrestrial")]))

# Proportion Marine species MAXNET TSS > 0.5
sdms_CV %>%
    dplyr::filter(ECO=="Marine",
                  algo=="MAXNET",
                  metric.eval=="TSS") %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(TSS = mean(calibration)) %>%
    dplyr::summarise(length(which(TSS > .5))/length(TSS))

# Proportion Marine species MAXNET TSS > 0.7
sdms_CV %>%
    dplyr::filter(ECO=="Marine",
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
shifts_ens_mar$ECO = "Marine"

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
shifts_ens_ter$ECO = "Terrestrial"

length(unique(shifts_ens_ter$Species))

matching_columns <- match(names(shifts_ens_mar),names(shifts_ens_ter))

shifts_ens <- rbind(shifts_ens_mar,
                    shifts_ens_ter[,names(shifts_ens_mar)])

length(unique(shifts_ens$Species))
length(unique(shifts_ens$Species[which(shifts_ens$ECO=="Terrestrial")]))
length(unique(shifts_ens$Species[which(shifts_ens$ECO=="Marine")]))
# N SAs
length(unique(shifts_ens$ID))

# N species
length(unique(shifts_ens$Species))
# N shifts
length(unique(shifts_ens$ID))
# ECO
table(shifts_ens$ECO)

### Correlation shifts % ----

## Terrestrial
ggpairs(shifts_ens %>%
            filter(ECO == "Terrestrial"), 
        columns = c(23:26, 9, 27:29))

ggpairs(shifts_ens %>%
            filter(ECO == "Marine"), 
        columns = c(23:26, 9, 27:29))


###
### Merge with bioshifts ----
#### plug in shift sdms
shifts_ens$new_ID <- paste(shifts_ens$ID, shifts_ens$Species, round(as.numeric(shifts_ens$START),0), round(as.numeric(shifts_ens$END),0), shifts_ens$Type)

biov1_pred <- merge(shifts_ens,
                    data.frame(biov1[,-which(names(biov1) %in% names(shifts_ens))],
                               new_ID=biov1$new_ID),
                    by = c("new_ID"),
                    all.x = TRUE)

### Plug in TSS ----
avg_TSS <- sdms_CV %>%
    dplyr::filter(metric.eval=="TSS") %>%
    group_by(Species) %>%
    summarise(TSS_vali = median(validation),
              TSS_cali = median(calibration))

all(biov1_pred$sp_name_std_v1 %in% avg_TSS$Species)

biov1_pred <- merge(
    biov1_pred, 
    data.frame(avg_TSS[,-which(names(avg_TSS) %in% names(biov1_pred))],
               Species=avg_TSS$Species), 
    by.x = "sp_name_std_v1",
    by.y = "Species",
    all.x = TRUE)

###
# Load velocity at edges ----
vel_edge <- read.csv(here::here("Data/vel_edge.csv"))
# N SAs
length(unique(sapply(vel_edge$ID, function(x) strsplit(x, " ")[[1]][1])))
vel_edge$new_ID <- vel_edge$ID
vel_edge <- vel_edge %>% select(-ID)

biov1_pred <- merge(biov1_pred, 
                    data.frame(vel_edge[,-which(names(vel_edge) %in% names(biov1_pred))],
                               new_ID=vel_edge$new_ID), 
                    by = "new_ID",
                    all.x = TRUE)


### select the observed shift variable  
names(biov1_pred)

# N species
length(unique(biov1_pred$sp_name_std_v1))

biov1_pred %>%
    group_by(ECO) %>%
    summarise(N=length(unique(sp_name_std_v1)))

###
### Compare velocities across edges ----
ggplot(biov1_pred %>%
           mutate(TE=bv.lat.1,O=bv.lat.mean,LE=bv.lat.9) %>%
           select(TE,O,LE) %>%
           gather(),
           aes(x=value,fill=key))+
    geom_histogram()+
    facet_wrap(.~key,ncol = 2,scales = "free")

ggplot(biov1_pred %>%
           mutate(TE=edge1,O=edge5,LE=edge9) %>%
           select(TE,O,LE) %>%
           gather(),
       aes(x=value,fill=key))+
    geom_histogram()+
    facet_wrap(.~key,ncol = 2,scales = "free")

###
### Shift predicted vs observed ----

#### Centroid
ggplot(biov1_pred %>%
           filter(Param == "O",
                  # (bv.lat.mean > 0 & SHIFT > 0) | (bv.lat.mean < 0 & SHIFT < 0)
                  ), 
       aes(x = bv.lat.mean, y = SHIFT, color = ECO, size = DUR))+
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
    facet_wrap(.~ECO,ncol = 2,scales = "free")

ggplot(biov1_pred %>%
           filter(Param == "O"), 
       aes(x = v.lat.mean, y = SHIFT, color = ECO, size = DUR))+
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
    facet_wrap(.~ECO,ncol = 2,scales = "free")

ggplot(biov1_pred %>%
           filter(Param == "O"), 
       aes(x = v.lat.mean, y = bv.lat.mean, color = ECO, size = DUR))+
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
    facet_wrap(.~ECO,ncol = 2,scales = "free")


#### Leading edge
ggplot(biov1_pred %>%
           filter(Param == "LE"), 
       aes(x = bv.lat.9, y = SHIFT, color = ECO, size = DUR))+
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
    facet_wrap(.~ECO,ncol = 2,scales = "free")

ggplot(biov1_pred %>%
           filter(Param == "LE"), 
       aes(x = edge9, y = SHIFT, color = ECO, size = DUR))+
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
    facet_wrap(.~ECO,ncol = 2,scales = "free")

ggplot(biov1_pred %>%
           filter(Param == "LE"), 
       aes(x = edge9, y = bv.lat.9, color = ECO, size = DUR))+
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
    facet_wrap(.~ECO, ncol = 2, scales = "free")

#### Trailing edge
ggplot(biov1_pred %>%
           filter(Param == "TE"), 
       aes(x = bv.lat.1, y = SHIFT, color = ECO, size = DUR))+
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
    facet_wrap(.~ECO,ncol = 2,scales = "free")

ggplot(biov1_pred %>%
           filter(Param == "TE"), 
       aes(x = edge1, y = SHIFT, color = ECO, size = DUR))+
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
    facet_wrap(.~ECO,ncol = 2,scales = "free")

ggplot(biov1_pred %>%
           filter(Param == "TE"), 
       aes(x = edge1, y = bv.lat.1, color = ECO, size = DUR))+
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
    facet_wrap(.~ECO,ncol = 2,scales = "free")

###
### Models ----

CE_model <- glmmTMB(SHIFT ~ v.lat.mean + bv.lat.mean + Type + 
                        (1 | Class), 
                     weights = TSS_vali,
                     data = biov1_pred %>% filter(Param == "O"))
summary(gau1_allW)

