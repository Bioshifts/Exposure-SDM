# Setup
# remotes::install_github("hughjonesd/ggmagnify")
# devtools::install_github(c("yihui/servr", "hafen/rmote")) 

# rmote::start_rmote() 

rm(list=ls())
gc()

list.of.packages <- c("dplyr","ggplot2","tune","GGally","ggmagnify","lme4","terra", "rgdal","maps","readr","rmote","tidyr")
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

# Get corrected shifts 
fetchGHdata(gh_account = "Bioshifts", 
            repo = "MethodologicalAdjustment", 
            path = "outputs/biov1_method_corrected_shifts_study_level.csv",
            output = here::here(work_dir,"Data/Bioshifts/biov1_method_corrected_shifts_study_level.csv"))

biov1_corr <- read.csv(
    here::here(work_dir,"Data/Bioshifts/biov1_method_corrected_shifts_study_level.csv"), 
    header = T)

biov1 <- merge(biov1, 
               biov1_corr %>%
                   select(c("Article_ID","Study_ID","Class","Type","SLDiff1")), 
               by = c("Article_ID","Study_ID","Class","Type"),
               all.x = TRUE)

biov1$SHIFT_cor <- abs(biov1$SHIFT) - biov1$SLDiff1


# N species
length(unique(biov1$sp_name_std_v1))



# New ID
biov1$new_ID <- paste(biov1$ID,
                      biov1$sp_name_std_v1,
                      round(biov1$START,0),
                      round(biov1$END,0),
                      biov1$Type,
                      biov1$Param)

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

shifts_ens <- list.files(here::here("Data","SHIFT","Mar"), 
                         pattern = " ens SA.csv", 
                         recursive = TRUE, full.names = TRUE)
length(shifts_ens)

shifts_ens <- lapply(shifts_ens, function(x) {
    read.csv(x)
})

shifts_ens <- data.table::rbindlist(shifts_ens)

shifts_ens_mar <- data.frame(shifts_ens)
shifts_ens_mar$ECO = "Marine"

length(unique(shifts_ens_mar$Species))

shifts_ens <- list.files(here::here("Data","SHIFT","Ter"), 
                         pattern = " ens SA.csv", 
                         recursive = TRUE, full.names = TRUE)

shifts_ens <- lapply(shifts_ens, function(x){
    read.csv(x)
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

# create ID for merge
START <- sapply(shifts_ens$time_period, function(x){
    as.numeric(strsplit(x,"-")[[1]][1])
})
END <- sapply(shifts_ens$time_period, function(x){
    as.numeric(strsplit(x,"-")[[1]][2])
})
shifts_ens$new_ID <- paste(shifts_ens$ID,
                           shifts_ens$Species,
                           START,
                           END,
                           shifts_ens$Type)

# add Taxonomic group
shifts_ens <- merge(shifts_ens, n_occ[,c("scientificName","Group")], 
                    by.x = "Species", 
                    by.y = "scientificName",
                    all.x = TRUE)

dim(shifts_ens)
length(unique(shifts_ens$Species))

# N Species
length(unique(shifts_ens$Species[which(shifts_ens$ECO=="Terrestrial")]))
length(unique(shifts_ens$Species[which(shifts_ens$ECO=="Marine")]))

###
### Get shifts_ens at the range positions ----

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

shifts_ens_ <- list()

all_IDS <- unique(shifts_ens$new_ID)

for(i in 1:length(all_IDS)){ cat(i, "from", length(all_IDS),"\r")
    
    ID_going <- all_IDS[i]
    
    tmp <- try({
        # subset shift i
        tmp <- shifts_ens %>% filter(new_ID == ID_going)
        
        # get avg shifts at range positions
        SHIFT_sdm_1 <- c(mean(tmp$nsQuantVelocity_quant0p01),
                         mean(tmp$nsCentroidVelocity),
                         mean(tmp$nsQuantVelocity_quant0p99))
        
        SHIFT_sdm_22 <- c(mean(tmp$nsQuantVelocity_quant0p05),
                          mean(tmp$nsCentroidVelocity),
                          mean(tmp$nsQuantVelocity_quant0p95))
        
        SHIFT_sdm_3 <- c(mean(tmp$nsQuantVelocity_quant0p1),
                         mean(tmp$nsCentroidVelocity),
                         mean(tmp$nsQuantVelocity_quant0p9))
        
        # get shift as linear coefficient of latitude at edges over time
        SHIFT_sdm_21 <- as.numeric(
            c(lm(tmp$nsQuantLat_quant0p01~tmp$START)[1][[1]][2],
              lm(tmp$nsCentroidLat~tmp$START)[1][[1]][2],
              lm(tmp$nsQuantLat_quant0p99~tmp$START)[1][[1]][2]))
        
        SHIFT_sdm_22 <- as.numeric(
            c(lm(tmp$nsQuantLat_quant0p05~tmp$START)[1][[1]][2],
              lm(tmp$nsCentroidLat~tmp$START)[1][[1]][2],
              lm(tmp$nsQuantLat_quant0p95~tmp$START)[1][[1]][2]))
        
        SHIFT_sdm_23 <- as.numeric(
            c(lm(tmp$nsQuantLat_quant0p1~tmp$START)[1][[1]][2],
              lm(tmp$nsCentroidLat~tmp$START)[1][[1]][2],
              lm(tmp$nsQuantLat_quant0p9~tmp$START)[1][[1]][2]))
        
        data.frame(Param = c("TE","O","LE"),
                   SHIFT_sdm_1 = SHIFT_sdm_1/1000, # from meters to km
                   SHIFT_sdm_22 = SHIFT_sdm_22/1000,
                   SHIFT_sdm_3 = SHIFT_sdm_3/1000,
                   SHIFT_sdm_21 = SHIFT_sdm_21/1000, 
                   SHIFT_sdm_22 = SHIFT_sdm_22/1000,
                   SHIFT_sdm_23 = SHIFT_sdm_23/1000,
                   tmp[1,c("Type", "Species", "time_period","new_ID","ECO","Group")])
    }, silent = TRUE)
    if(class(tmp) == "try-error"){
        shifts_ens_[[i]] <- NULL 
    } else {
        shifts_ens_[[i]] <- tmp 
    }
    
}


shifts_ens <- data.table::rbindlist(shifts_ens_)
shifts_ens <- data.frame(shifts_ens)

shifts_ens$new_ID <- paste(shifts_ens$new_ID, shifts_ens$Param)
# shifts_ens <- na.omit(shifts_ens)
dim(shifts_ens)

# N species
length(unique(shifts_ens$Species))
# N shifts
length(unique(shifts_ens$new_ID))
# ECO
table(shifts_ens$ECO)

table(shifts_ens$ECO, shifts_ens$Param)

### Correlation shifts % ----

## Terrestrial
ggpairs(shifts_ens %>%
            filter(ECO == "Terrestrial") %>%
            filter(Param != "O"), 
        columns = c(2:7),
        ggplot2::aes(colour=Param))

ggpairs(shifts_ens %>%
            filter(ECO == "Marine") %>%
            filter(Param != "O"), 
        columns = c(2:7),
        ggplot2::aes(colour=Param))


###
### Merge with bioshifts ----
#### plug in shift sdms
biov1_pred <- merge(
    biov1, shifts_ens[,c(2:7,11,13)], 
    by = "new_ID")

### Plug in TSS ----
avg_TSS <- sdms_CV %>%
    dplyr::filter(metric.eval=="TSS") %>%
    group_by(Species) %>%
    summarise(TSS_vali = median(validation),
              TSS_cali = median(calibration))

all(biov1_pred$sp_name_std_v1 %in% avg_TSS$Species)

biov1_pred <- merge(
    biov1_pred, avg_TSS, 
    by.x = "sp_name_std_v1",
    by.y = "Species",
    all.x = TRUE)

###
# Load velocity at edges ----
vel_edge <- read.csv(here::here("Data/vel_edge.csv"))
# N SAs
length(unique(sapply(vel_edge$ID, function(x) strsplit(x, " ")[[1]][1])))
vel_edge$new_ID2 <- vel_edge$ID
vel_edge <- vel_edge %>% select(-ID)

biov1_pred$new_ID2 <- sapply(biov1_pred$new_ID, function(x) paste(strsplit(x," ")[[1]][1:5],collapse = " "))

biov1_pred <- merge(biov1_pred, 
                    vel_edge, 
                    by = "new_ID2",
                    all.x = TRUE)
biov1_pred <- biov1_pred %>% select(!new_ID2)

# Fix edge
biov1_pred$clim_vel_edge <- NA 
edges <- sapply(biov1_pred$new_ID, function(x) strsplit(x," ")[[1]][6])
biov1_pred$clim_vel_edge[edges=="LE"]  <- biov1_pred$edge95[edges=="LE"]
biov1_pred$clim_vel_edge[edges=="TE"]  <- biov1_pred$edge05[edges=="TE"]
biov1_pred$clim_vel_edge[edges=="O"]  <- biov1_pred$edge5[edges=="O"]

### select the observed shift variable  
biov1_pred$SHIFT_obs <- biov1_pred$SHIFT
biov1_pred$bioclimatic_vel <- biov1_pred$SHIFT_sdm_22
biov1_pred$clim_vel <- biov1_pred$vel_mat

ggpairs(biov1_pred %>%
            filter(ECO == "Terrestrial"),
        columns = c(41,42,39),
        ggplot2::aes(colour=Param))+
    theme_classic()

ggpairs(biov1_pred %>%
            filter(ECO == "Marine"),
        columns = c(41,42,39),
        ggplot2::aes(colour=Param))+
    theme_classic()


### Stats N species ----
# N species v1
length(unique(biov1$sp_name_std_v1))
# N species sdms
length(unique(shifts_ens$Species))
# N species with SDMs in bioshifts
length(which(unique(shifts_ens$Species) %in% unique(biov1$sp_name_std_v1)))
# N species velocity edge
sp_vel_edge <- unique(sapply(vel_edge$new_ID, function(x) strsplit(x," ")[[1]][2]))
length(unique(sp_vel_edge))
# N species with SDMs in bioshifts
length(which(sp_vel_edge %in% unique(biov1$sp_name_std_v1)))

# shifts in v1
length(unique(biov1$new_ID))
# shifts in SDMs
length(unique(shifts_ens$new_ID))
# N shifts with SDMs in bioshifts
length(which(unique(shifts_ens$new_ID) %in% unique(biov1$new_ID)))
# N species
length(unique(shifts_ens$Species[which(shifts_ens$ECO == "Terrestrial")]))
length(unique(shifts_ens$Species[which(shifts_ens$ECO == "Terrestrial")]))
# N shifts
length(unique(shifts_ens$new_ID[which(shifts_ens$ECO == "Terrestrial")]))
length(unique(shifts_ens$new_ID[which(shifts_ens$ECO == "Terrestrial")]))

spincommon <- intersect(biov1$sp_name_std_v1,shifts_ens$Species)

table(biov1$ECO, biov1$Param)
table(shifts_ens$ECO, shifts_ens$Param)
table(biov1_pred$ECO, biov1_pred$Param)

table(biov1_pred$ECO, biov1_pred$Param, biov1_pred$Class)


### Missing shifts ----
# many cases in which shifts in v1 are not found in sdms

all_shifts_v1 <- unique(biov1$new_ID)
all_shifts_sh <- unique(shifts_ens$new_ID)

missing_shifts <- data.frame(
    new_ID = all_shifts_v1[which(!all_shifts_v1 %in% all_shifts_sh)])
missing_shifts <- merge(missing_shifts,
                        biov1[,c("ECO","new_ID")],
                        by = "new_ID",
                        all.x = TRUE)

dim(missing_shifts)

### Missing shifts sps with sdms ----

all_shifts_v1 <- unique(biov1$new_ID[which(biov1$sp_name_std_v1 %in% spincommon)])
all_shifts_sh <- unique(shifts_ens$new_ID[which(shifts_ens$Species %in% spincommon)])

missing_shifts <- data.frame(
    new_ID = all_shifts_v1[which(!all_shifts_v1 %in% all_shifts_sh)])
missing_shifts <- merge(missing_shifts,
                        biov1[,c("ECO","new_ID")],
                        by = "new_ID",
                        all.x = TRUE)
dim(missing_shifts)



###
### Geographical pattern ----

# mundi <- terra::vect(rnaturalearth::ne_coastline(scale = 10, returnclass = "sp"))
# 
# # get raster bioshifts shp files study areas
# get_raster_bioshifts = "NO"
# if(get_raster_bioshifts=="YES"){
#     my_ext = terra::ext(mundi)
#     my_crs = crs(mundi)
#     
#     # empty raster for study areas
#     rast_biov1 <- terra::rast(my_ext, crs = my_crs, res = 0.5)
#     values(rast_biov1) <- 0
#     
#     # empty raster for N shifts
#     rast_biov1_sh <- terra::rast(my_ext, crs = my_crs, res = 0.5)
#     values(rast_biov1) <- 0
#     
#     fgdb <- "C:/Users/brunn/NextCloud/Bioshifts/Study_Areas.gdb"
#     fc_list <- rgdal::ogrListLayers(fgdb)
#     fc_list <- fc_list[which(fc_list %in% unique(biov1_pred$ID))]
#     
#     for(i in 1:length(fc_list)){ cat("\r",i,"from",length(fc_list))
#         tmp = terra::vect(sf::st_read(fgdb, layer=fc_list[i]))
#         tmp = terra::cells(rast_biov1,tmp)
#         tmp_cell = tmp[,2]
#         tmp_vals = rast_biov1[tmp_cell][,1]
#         rast_biov1[tmp_cell] = tmp_vals+1
#         
#         # get N shifts at the study area
#         N <- nrow(biov1_pred[which(biov1_pred$ID == fc_list[i]),])
#         rast_biov1_sh[tmp_cell] = tmp_vals + N
#     }
#     names(rast_biov1) <- names(rast_biov1_sh) <- "SA"
#     rast_biov1[rast_biov1==0] <- NA
#     rast_biov1_sh[rast_biov1_sh==0] <- NA
#     
#     writeRaster(rast_biov1, "Data/raster_bioshifts_SA.tif", overwrite = TRUE)
#     writeRaster(rast_biov1_sh, "Data/raster_bioshifts_N_SA.tif", overwrite = TRUE)
#     
# } else {
#     rast_biov1 <- terra::rast(here::here("Data/raster_bioshifts_SA.tif"))
#     rast_biov1_sh <- terra::rast(here::here("Data/raster_bioshifts_N_SA.tif"))
# }
# 
# # plot
# mundi_sf <- sf::st_as_sf(mundi)
# 
# test_df <- terra::as.data.frame(rast_biov1, xy = TRUE)
# colnames(test_df) <- c("x", "y", "value")
# 
# p1 = ggplot() +  
#     geom_tile(data=test_df, aes(x=x, y=y, fill=value)) + 
#     ggtitle("N study areas")+
#     geom_sf(data=mundi_sf) +
#     theme_minimal() +
#     scale_fill_viridis_c()+
#     guides(fill = guide_colourbar(title = ""))+
#     xlab("")+ylab("")
# 
# test_df <- terra::as.data.frame(rast_biov1_sh, xy = TRUE)
# colnames(test_df) <- c("x", "y", "value")
# 
# p2 = ggplot() +  
#     geom_tile(data=test_df, aes(x=x, y=y, fill=value)) + 
#     ggtitle("N range shifts")+
#     geom_sf(data=mundi_sf) +
#     theme_minimal() +
#     scale_fill_viridis_c()+
#     guides(fill = guide_colourbar(title = ""))+
#     xlab("")+ylab("")

gridExtra::grid.arrange(p1,p2,ncol=1)

###
### Shift pattern ----

longdata <- biov1_pred[,c("ECO","Param","clim_vel_edge","clim_vel","bioclimatic_vel")] %>%
    tidyr::gather("Type", "SHIFT", -c(Param,ECO))

ggplot(longdata, aes(x=SHIFT, fill = Param))+
    geom_density(alpha=.5)+
    theme_bw()+
    facet_wrap(ECO~Type,scales = "free")

###
### Shift predicted vs observed ----

longdata <- biov1_pred[,c("ECO","Param","clim_vel_edge","clim_vel","bioclimatic_vel","SHIFT_obs","DUR")] %>%
    tidyr::gather("Metric", "value", -c(Param,ECO,SHIFT_obs,DUR))

## all together
ggplot(longdata %>%
           filter(Param == "O"), 
       aes(x = as.numeric(value), y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Centroid shifts")+
    geom_point(alpha = .3)+
    xlab("Velocity metric (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(Metric~ECO,ncol = 2,scales = "free")

ggplot(longdata %>%
           filter(Param == "LE"), 
       aes(x = as.numeric(value), y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Leading edge shifts")+
    geom_point(alpha = .3)+
    xlab("Velocity metric (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(Metric~ECO,ncol = 2,scales = "free")

ggplot(longdata %>%
           filter(Param == "TE"), 
       aes(x = as.numeric(value), y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Traling edge shifts")+
    geom_point(alpha = .3)+
    xlab("Velocity metric (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(Metric~ECO,ncol = 2,scales = "free")


## one by one
ggplot(biov1_pred[,c("ECO","Param","SHIFT_obs","bioclimatic_vel")] %>%
           tidyr::gather("Type", "SHIFT", -c(Param,ECO)), 
       aes(x = Type, y = SHIFT))+
    ylab("Range shift (km/yr)")+
    geom_violin(trim = TRUE, draw_quantiles = c(.05, .5, .95))+
    theme_classic()+
    facet_wrap(.~ECO, scales = "free")

ggplot(biov1_pred[,c("ECO","Param","SHIFT_obs","bioclimatic_vel")] %>%
           tidyr::gather("Type", "SHIFT", -c(Param,ECO)), 
       aes(x = Type, y = SHIFT))+
    ylab("Range shift (km/yr)")+
    geom_violin(trim = TRUE, draw_quantiles = c(.05, .5, .95))+
    theme_classic()+
    facet_wrap(ECO~Param, scales = "free")

ggplot(biov1_pred[,c("ECO","Param","SHIFT_obs","bioclimatic_vel")] %>%
           tidyr::gather("Type", "SHIFT", -c(Param,ECO)), 
       aes(fill = Type, x = SHIFT))+
    xlab("Range shift (km/yr)")+
    geom_density(alpha=.5)+
    theme_classic()+
    facet_wrap(.~ECO, scales = "free")

ggplot(biov1_pred[,c("ECO","Param","SHIFT_obs","bioclimatic_vel")] %>%
           tidyr::gather("Type", "SHIFT", -c(Param,ECO)) , 
       aes(fill = Type, x = SHIFT))+
    xlab("Range shift (km/yr)")+
    geom_density(alpha=.5)+
    theme_classic()+
    facet_wrap(ECO~Param, scales = "free")

ggplot(
    biov1_pred, 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("All shifts")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()+
    facet_wrap(.~ECO)

# Optimum
M_O_bio <- ggplot(
    biov1_pred %>%
        filter(Param == "O"), 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Centroid")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(.~ECO, scales = "free")

M_O_vel <- ggplot(
    biov1_pred %>%
        filter(Param == "O"), 
    aes(x = clim_vel, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Centroid")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(.~ECO, scales = "free")

M_O_vel_ed <- ggplot(
    biov1_pred %>%
        filter(Param == "O"), 
    aes(x = clim_vel_edge, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Centroid")+
    geom_point(alpha = .3)+
    xlab("Climate velocity edge (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(.~ECO, scales = "free")

# LE
M_LE_bio <- ggplot(
    biov1_pred %>%
        filter(Param == "LE" ), 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(.~ECO,scales = "free")

M_LE_vel <- ggplot(
    biov1_pred %>%
        filter(Param == "LE" ), 
    aes(x = clim_vel, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Climate velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(.~ECO,scales = "free")

M_LE_vel_ed <- ggplot(
    biov1_pred %>%
        filter(Param == "LE" ), 
    aes(x = clim_vel_edge, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Climate velocity edge (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(.~ECO,scales = "free")

# TE
M_TE_bio <- ggplot(
    biov1_pred %>%
        filter(Param == "TE"), 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(.~ECO,scales = "free")

M_TE_vel <- ggplot(
    biov1_pred %>%
        filter(Param == "TE" ), 
    aes(x = clim_vel, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Climate velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(.~ECO,scales = "free")

M_TE_vel_ed <- ggplot(
    biov1_pred %>%
        filter(Param == "TE" ), 
    aes(x = clim_vel_edge, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Climate velocity edge (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    # tune::coord_obs_pred()+
    facet_wrap(.~ECO,scales = "free")

M_O_bio
M_O_vel
M_O_vel_ed
M_LE_bio
M_LE_vel
M_LE_vel_ed
M_TE_bio
M_TE_vel
M_TE_vel_ed

## N species shifting same direction of bioclimatic velocity

table(sign(biov1_pred$SHIFT_obs) == sign(biov1_pred$bioclimatic_vel), biov1_pred$ECO)
table(sign(biov1_pred$SHIFT_obs) == sign(biov1_pred$clim_vel), biov1_pred$ECO)
table(sign(biov1_pred$SHIFT_obs) == sign(biov1_pred$clim_vel_edge), biov1_pred$ECO)



# Optimum
M_O <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial" & Param == "O"), 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = Group, size = DUR))+
    ggtitle("Centroid")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic()+
    # tune::coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

# LE
M_LE <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial" & Param == "LE"), 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = Group, size = DUR))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic()+
    # tune::coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

# TE
M_TE <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial" & Param == "TE"), 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = Group, size = DUR))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic()+
    # tune::coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")


M_O
M_LE
M_LE + lims(x=c(-5,2.5),y=c(-10,25))
M_TE

# gridExtra::grid.arrange(grobs = list(M_O,
#                                      M_LE + lims(x=c(-2.5,2.5),y=c(-1,22)),
#                                      M_TE),ncol = 3)

# just class with > 20 shifts
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>20 & tmp$Var2=="Terrestrial",]
tmp <- tmp[order(tmp$Var1),]

plots <- list()
for(i in 1:length(tmp$Var1)){
    
    plots[[i]] <- ggplot(
        biov1_pred %>%
            filter(ECO == "Terrestrial" & Param == tmp$Var3[i] & Group == tmp$Var1[i]), 
        aes(x = bioclimatic_vel, y = SHIFT_obs, color = Class))+
        ggtitle(paste(tmp$Var1[i],tmp$Var3[i]))+
        geom_point(alpha = .5)+
        xlab("Bioclimatic velocity (km/yr)")+
        ylab("Observed shift (km/yr)")+
        theme_classic() +
        # tune::coord_obs_pred()+
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
        geom_hline(yintercept = 0, linetype = "dashed")+
        geom_vline(xintercept = 0, linetype = "dashed")
}
names(plots) <- paste(tmp$Var1,tmp$Var3)

gridExtra::grid.arrange(grobs = plots,ncol = 2)



# Terrestrials

tmp <- biov1_pred %>%
    filter(ECO == "Terrestrial")

table(sign(tmp$SHIFT_obs) == sign(tmp$bioclimatic_vel), tmp$Param)

# Optimum
T_O <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial" & Param == "O"), 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = Group))+
    ggtitle("Centroid")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    # tune::coord_obs_pred()+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

# LE
T_LE <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial" & Param == "LE"), 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = Group))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    # tune::coord_obs_pred()+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

# TE
T_TE <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial" & Param == "TE"), 
    aes(x = bioclimatic_vel, y = SHIFT_obs, color = Group))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    # tune::coord_obs_pred()+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")


T_O
T_LE
T_TE

# gridExtra::grid.arrange(grobs = list(T_O,T_LE,T_TE),ncol = 3)

# just class with > 10 shifts
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10 & tmp$Var2=="Terrestrial",]
tmp <- tmp[order(tmp$Var1),]


plots <- list()
for(i in 1:length(tmp$Var1)){
    
    plots[[i]] <- ggplot(
        biov1_pred %>%
            filter(ECO == "Terrestrial" & Param == tmp$Var3[i] & Group == tmp$Var1[i]), 
        aes(x = bioclimatic_vel, y = SHIFT_obs, color = Class))+
        ggtitle(paste(tmp$Var1[i],tmp$Var3[i]))+
        geom_point(alpha = .8)+
        xlab("Bioclimatic velocity (km/yr)")+
        ylab("Observed shift (km/yr)")+
        theme_classic() +
        # tune::coord_obs_pred()+
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
        geom_hline(yintercept = 0, linetype = "dashed")+
        geom_vline(xintercept = 0, linetype = "dashed")
}
names(plots) <- paste(tmp$Var1,tmp$Var3)

gridExtra::grid.arrange(grobs = plots,ncol = 2)


###
### calculate lags ----
biov1_pred$bioclimatic_vel_lag <- biov1_pred$bioclimatic_vel - biov1_pred$SHIFT_obs
biov1_pred$clim_vel_lag <- biov1_pred$clim_vel - biov1_pred$SHIFT_obs
biov1_pred$clim_vel_edge_lag <- biov1_pred$clim_vel_edge - biov1_pred$SHIFT_obs


# fix lag
# same direction: change lag sign when vel is negative
pos <- which(biov1_pred$bioclimatic_vel<0)
biov1_pred$bioclimatic_vel_lag[pos] <- biov1_pred$bioclimatic_vel_lag[pos] * -1

pos <- which(biov1_pred$clim_vel<0)
biov1_pred$clim_vel_lag[pos] <- biov1_pred$clim_vel_lag[pos] * -1

pos <- which(biov1_pred$clim_vel_edge<0)
biov1_pred$clim_vel_edge_lag[pos] <- biov1_pred$clim_vel_edge_lag[pos] * -1

# opposite direction: lag = vel
biov1_pred$bioclimatic_vel_lag2 <- biov1_pred$bioclimatic_vel_lag
pos <- which((biov1_pred$bioclimatic_vel>0 & biov1_pred$SHIFT_obs<0)|
                 biov1_pred$bioclimatic_vel<0 & biov1_pred$SHIFT_obs>0)
biov1_pred$bioclimatic_vel_lag2[pos] <- biov1_pred$bioclimatic_vel[pos]
biov1_pred$clim_vel_lag2[pos] <- biov1_pred$bioclimatic_vel[pos]
biov1_pred$clim_vel_edge_lag2[pos] <- biov1_pred$bioclimatic_vel[pos]


ggpairs(biov1_pred %>%
            filter(ECO == "Terrestrial"),
        columns = c("bioclimatic_vel_lag","clim_vel_lag","clim_vel_edge_lag"),
        ggplot2::aes(colour=Param))+
    theme_classic()

ggpairs(biov1_pred %>%
            filter(ECO == "Marine"),
        columns = c("bioclimatic_vel_lag","clim_vel_lag","clim_vel_edge_lag"),
        ggplot2::aes(colour=Param))+
    theme_classic()


## % species shifting faster vs slower than expected
table(biov1_pred$bioclimatic_vel_lag > 0, biov1_pred$ECO)
table(biov1_pred$clim_vel_lag > 0, biov1_pred$ECO)
table(biov1_pred$clim_vel_edge_lag > 0, biov1_pred$ECO)


## Terrestrials
# just class with > 10 shifts
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10 & tmp$Var2=="Terrestrial",]
plot_data <- biov1_pred %>%
    dplyr::filter(ECO == "Terrestrial" & Group %in% tmp$Var1) %>%
    dplyr::select(Group,Param,bioclimatic_vel_lag,clim_vel_lag,clim_vel_edge_lag) %>%
    gather("LagMetric", "Value", -c(Group,Param))

# limits <- biov1_pred %>%
#     dplyr::filter(ECO == "Terrestrial" & (Group %in% tmp$Var1)) %>%
#     dplyr::reframe(lag=quantile(bioclimatic_vel_lag,c(0.01, 0.99))) 

ggplot(plot_data, 
       aes(x = Value, y = Group, color = LagMetric))+
    # geom_violin(trim = TRUE, draw_quantiles = c(0.05,0.5,0.95)) +
    ggtitle("Terrestrial")+
    geom_point(alpha=.1)+
    geom_boxplot(outlier.alpha = 1, alpha=.5)+
    # scale_x_continuous(limits = limits[,1])+
    xlab("Lag = Velocity metric - Observed range shift")+
    ylab("")+
    geom_vline(xintercept = 0, linetype = "longdash")+
    facet_wrap(Param~., scales = "free")+
    theme_classic()

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & Group %in% tmp$Var1), 
       aes(x = bioclimatic_vel_lag, fill = Group, color = Group))+
    # geom_violin(trim = TRUE, draw_quantiles = c(0.05,0.5,0.95)) +
    geom_density(alpha = 0.3)+
    scale_x_continuous(limits = limits[,1])+
    xlab("Lag = Bioclimatic velocity - Range shift velocity")+
    ylab("Density")+
    geom_vline(xintercept = 0, linetype = "longdash")+
    facet_wrap(Param~., scales = "free")+
    theme_classic()

## Terrestrials
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10 & tmp$Var2=="Terrestrial",]
plot_data <- biov1_pred %>%
    dplyr::filter(ECO == "Terrestrial" & Group %in% tmp$Var1) %>%
    dplyr::select(Group,Param,bioclimatic_vel_lag,clim_vel_lag,clim_vel_edge_lag) %>%
    gather("LagMetric", "Value", -c(Group,Param))

# limits <- biov1_pred %>%
#     dplyr::filter(ECO == "Terrestrial" & (Group %in% tmp$Var1)) %>%
#     dplyr::reframe(lag=quantile(bioclimatic_vel_lag,c(0.01, 0.99))) 

ggplot(plot_data, 
       aes(x = Value, y = Group, color = LagMetric))+
    # geom_violin(trim = TRUE, draw_quantiles = c(0.05,0.5,0.95)) +
    ggtitle("Terrestrial")+
    geom_point(alpha=.1)+
    geom_boxplot(outlier.alpha = 1, alpha=.5)+
    # scale_x_continuous(limits = limits[,1])+
    xlab("Lag = Velocity metric - Observed range shift")+
    ylab("")+
    geom_vline(xintercept = 0, linetype = "longdash")+
    facet_wrap(Param~., scales = "free")+
    theme_classic()

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & Group %in% tmp$Var1), 
       aes(x = bioclimatic_vel_lag, fill = Group, color = Group))+
    # geom_violin(trim = TRUE, draw_quantiles = c(0.05,0.5,0.95)) +
    geom_density(alpha = 0.3)+
    scale_x_continuous(limits = limits[,1])+
    xlab("Lag = Bioclimatic velocity - Range shift velocity")+
    ylab("")+
    geom_vline(xintercept = 0, linetype = "longdash")+
    facet_wrap(~Param, scales = "free")+
    theme_classic()

###
# Summarise velocities in study areas ----

plot_data <- biov1_pred %>%
    group_by(ID,ECO) %>%
    dplyr::summarise(bioclimatic_vel = mean(bioclimatic_vel,na.rm=TRUE),
                     clim_vel = mean(clim_vel,na.rm=TRUE),
                     clim_vel_edge = mean(clim_vel_edge,na.rm=TRUE)) %>%
    gather("Study", "Value", -c(ID,ECO))

ggplot(plot_data, 
       aes(x = Value, y = Study, color = LagMetric))+
    # geom_violin(trim = TRUE, draw_quantiles = c(0.05,0.5,0.95)) +
    ggtitle("Terrestrial")+
    geom_point(alpha=.1)+
    geom_boxplot(outlier.alpha = 1, alpha=.5)+
    # scale_x_continuous(limits = limits[,1])+
    xlab("Lag = Velocity metric - Observed range shift")+
    ylab("")+
    geom_vline(xintercept = 0, linetype = "longdash")+
    facet_wrap(Param~., scales = "free")+
    theme_classic()




###
## Model predict range shift
all_models <- expand.grid(ECOs = c("Terrestrial", "Marine"),
                          Param = c("TE","LE","O"))

models <- lapply(1:nrow(all_models), function(x){
    tmp = lmer(SHIFT_obs~
                   bioclimatic_vel + clim_vel_edge + clim_vel + (1|Article_ID),
               data = biov1_pred %>%
                   filter(ECO == all_models$ECOs[x],
                          Param == all_models$Param[x]))
    tmp_conf <- confint(tmp)
    tmp <- summary(tmp)
    tmp <- data.frame(tmp$coefficients)
    tmp <- cbind(tmp, tmp_conf[rownames(tmp),])
    data.frame(tmp, 
               vars = rownames(tmp), 
               ECO = all_models$ECOs[x], Param = all_models$Param[x])
})
models <- data.table::rbindlist(models)
models$sig <- 0.5
models$sig[which((models$X2.5.. > 0 & models$X97.5.. > 0) | (models$X2.5.. < 0 & models$X97.5.. < 0))] <- 1

ggplot(models %>%
           filter(!vars == "(Intercept)"),
       aes(x = Estimate, y = vars, color = Param, alpha = sig))+
    geom_point()+
    geom_linerange(aes(xmin = X2.5.., xmax = X97.5..))+
    theme_bw()+
    geom_vline(xintercept = 0)+
    facet_wrap(ECO~Param,scales = "free_x")


models <- lapply(1:nrow(all_models), function(x){
    tmp = lmer(SHIFT_obs~
                   bioclimatic_vel + clim_vel_edge + clim_vel + (1|Article_ID),
               data = biov1_pred %>%
                   filter(ECO == all_models$ECOs[x],
                          Param == all_models$Param[x]))
    tmp_conf <- confint(tmp)
    tmp <- summary(tmp)
    tmp <- data.frame(tmp$coefficients)
    tmp <- cbind(tmp, tmp_conf[rownames(tmp),])
    data.frame(tmp, 
               vars = rownames(tmp), 
               ECO = all_models$ECOs[x], Param = all_models$Param[x])
})
models <- data.table::rbindlist(models)
models$sig <- 0.5
models$sig[which((models$X2.5.. > 0 & models$X97.5.. > 0) | (models$X2.5.. < 0 & models$X97.5.. < 0))] <- 1

ggplot(models %>%
           filter(!vars == "(Intercept)"),
       aes(x = Estimate, y = vars, color = Param, alpha = sig))+
    geom_point()+
    geom_linerange(aes(xmin = X2.5.., xmax = X97.5..))+
    theme_bw()+
    geom_vline(xintercept = 0)+
    facet_wrap(ECO~Param,scales = "free_x")




###
### Lag vs duration ----

### Terrestrial
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10 & tmp$Var2=="Terrestrial",]
plot_data <- biov1_pred %>%
    dplyr::filter(ECO == "Terrestrial" & Group %in% tmp$Var1) %>%
    dplyr::select(Group,Param,bioclimatic_vel_lag,clim_vel_lag,clim_vel_edge_lag,DUR) %>%
    gather("LagMetric", "Value", -c(Group,Param,DUR))

ggplot(plot_data, aes(x=DUR,y=Value))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("Duration (years)")+
    theme_bw()+
    facet_wrap(LagMetric~Param)


ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial"),
       aes(x = DUR, y = bioclimatic_vel_lag))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("Duration (years)")+
    theme_bw()+
    facet_wrap(.~Param)

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial"),
       aes(x = START, y = bioclimatic_vel_lag))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("Start year")+
    theme_bw()+
    facet_wrap(.~Param)

### Terrestrial

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & bioclimatic_vel_lag>-100),
       aes(x = DUR, y = bioclimatic_vel_lag, size=abs(SHIFT_obs)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("Duration (years)")+
    theme_bw()+
    facet_wrap(.~Param)

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & bioclimatic_vel_lag>-100),
       aes(x = START, y = bioclimatic_vel_lag, size=abs(SHIFT_obs)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("Start year")+
    theme_bw()+
    facet_wrap(.~Param)



### Lag vs TSS ----

### Terrestrial

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial"),
       aes(x = TSS_vali, y = abs(bioclimatic_vel_lag),size=abs(bioclimatic_vel)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("TSS validation")+
    theme_bw()+
    facet_wrap(.~Param)

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial"),
       aes(x = TSS_cali, y = abs(bioclimatic_vel_lag),size=abs(bioclimatic_vel)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("TSS calibration")+
    theme_bw()+
    facet_wrap(.~Param)

### Terrestrial

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial"),
       aes(x = TSS_vali, y = abs(bioclimatic_vel_lag),size=abs(bioclimatic_vel)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("TSS validation")+
    theme_bw()+
    facet_wrap(.~Param, scales = "free")

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial"),
       aes(x = TSS_cali, y = abs(bioclimatic_vel_lag),size=abs(bioclimatic_vel)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("TSS calibration")+
    theme_bw()+
    facet_wrap(.~Param, scales = "free")


###
### Predicted shifts vs Temperature velocities ----

## Does bioclimatic velocity explains shifts better then temperature velocity?

### Terrestrials
# all groups with > > 10 shifts
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10,]

mods <- data.frame()

for(i in 1:nrow(tmp)){
    
    mod_vel_i <- lm(
        SHIFT_obs~velocity,
        data = biov1_pred %>%
            filter(ECO == tmp$Var2[i] & 
                       Param == tmp$Var3[i] & 
                       Group == tmp$Var1[i]))
    
    # mod_vel_i <- lm(
    #     SHIFT_obs~vel_var,
    #     data = biov1_pred %>%
    #         filter(ECO == tmp$Var2[i] & 
    #                    Param == tmp$Var3[i] & 
    #                    Group == tmp$Var1[i]) %>%
    #         mutate(vel_var = ifelse(tmp$Var2[i]=="Terrestrial", 
    #                                 v.lat.median.sst,
    #                                 v.lat.median.mat)))
    
    mod_vel_i_sum <- summary(mod_vel_i)
    
    mod_bio_i <- lm(
        SHIFT_obs~bioclimatic_vel,
        data = biov1_pred %>%
            filter(ECO == tmp$Var2[i] & 
                       Param == tmp$Var3[i] & 
                       Group == tmp$Var1[i]))
    
    mod_bio_i_sum <- summary(mod_bio_i)
    
    R2_vel <- mod_vel_i_sum$r.squared
    R2_bio <- mod_bio_i_sum$r.squared
    coeff_vel <- try(mod_vel_i_sum$coefficients[2,1])
    coeff_bio <- mod_bio_i_sum$coefficients[2,1]
    
    mods <- rbind(mods,
                  data.frame(R2_vel= R2_vel,
                             R2_bio = R2_bio,
                             coeff_vel = coeff_vel,
                             coeff_bio = coeff_bio))
}
mods$coeff_vel <- as.numeric(mods$coeff_vel)
mods <- cbind(mods,tmp[,1:3])

ggplot(mods, aes(x = R2_vel, y = R2_bio))+
    geom_point()+
    facet_wrap(.~Var2, scales = "free")+
    theme_classic()

ggplot(mods, aes(x = coeff_vel, y = coeff_bio))+
    geom_point()+
    facet_wrap(.~Var2, scales = "free")+
    theme_classic()

## Terrestrial
p1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_21, color = Param))+
    ggtitle("TE = 1%, LE = 99%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic()+
    # tune::coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

p2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial"), 
    aes(x = v.lat.mean, y = bioclimatic_vel, color = Param))+
    ggtitle("TE = 5%, LE = 95%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic()+
    # tune::coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

p3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_23, color = Param))+
    ggtitle("TE = 10%, LE = 90%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic()+
    # tune::coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

gridExtra::grid.arrange(p1,p2,p3,ncol=3,top="Terrestrial")



## Terrestrial
p1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_21, color = Param))+
    ggtitle("TE = 1%, LE = 99%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic()+
    # tune::coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

p2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial"), 
    aes(x = v.lat.mean, y = bioclimatic_vel, color = Param))+
    ggtitle("TE = 5%, LE = 95%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic()+
    # tune::coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

p3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_23, color = Param))+
    ggtitle("TE = 10%, LE = 90%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic()+
    # tune::coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

gridExtra::grid.arrange(p1,p2,p3,ncol=3,top="Terrestrials")


ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial", Param == "O"), 
    aes(x = velocity, y = bioclimatic_vel))+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

#### compare bioclimatic vel with temp vel over the entire range

# expectations
# climate velocity have been used. how it is calculated. can be applied by across species but does not account for species specific sensitivities >>> SDMs

# 1) better predictions for Terrestrial than terrestrials because landscape is more connected - dispersal constrains in land. Terrestrial more in eq with the environment
# 2) but we would expect lags associated with dispersal and persistance (demography, life history)

# sdms methods

# bioshifts

# plot all param >>> only LE, TE O

# expectations by group >> fishes should track better than crustacea
