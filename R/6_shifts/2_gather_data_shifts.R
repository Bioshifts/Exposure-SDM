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
library(terra)
library(rgdal)
library(maps)
library(readr)

###
# computer = "muse"
computer = "personal"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}

work_dir <- getwd()
shift_dir <- here::here(work_dir,"Data/SHIFT")
sdm_dir <- here::here(work_dir,"Data/SDMs")

###
### Load in Bioshifts ----
biov1 <- read.csv("Data/Bioshifts/biov1_fixednames.csv", header = T)

# N species
length(unique(biov1$sp_name_std_v1))

biov1 <- biov1 %>% 
    mutate(
        Type = case_when(
            Type=="HOR" ~ "LAT",
            TRUE ~ as.character(Type)),
        velocity = ifelse(Type == "LAT", v.lat.mean, v.ele.mean)) %>%
    filter(Type == "LAT")

# N species
length(unique(biov1$sp_name_std_v1))

biov1$new_ID <- paste(biov1$ID,
                      biov1$sp_name_std_v1,
                      round(biov1$START,0),
                      round(biov1$END,0),
                      biov1$Type,
                      biov1$Param)

any(duplicated(biov1$new_ID))
# test = which(duplicated(biov1$new_ID))
# biov1$new_ID[test][1]
# biov1[which(biov1$new_ID == biov1$new_ID[test][1]),]

biov1$ECO <- ifelse(biov1$ECO=="M", "Marine", "Terrestrial")

table(biov1$ECO, biov1$Param)
dim(biov1)

###
### Load occ data ----
# add group
n_occ <- read_csv("Data/n_occ2.csv")
n_occ$scientificName <- gsub(" ","_",n_occ$scientificName)

###
### Load in SDM CV ----
sdms_CV <- read.csv(here::here(sdm_dir,"sdms_CV.csv"))
sdms_CV <- sdms_CV %>% dplyr::filter(metric.eval == "TSS")
# sdms_CV <- sdms_CV %>% dplyr::filter(validation > 0)

sdms_CV$Species <- sapply(sdms_CV$full.name, function(x){
    tmp = strsplit(x,"_")[[1]][1]
    gsub("[.]","_",tmp)
})
sdms_CV$ECO <- ifelse(sdms_CV$ECO=="Mar", "Marine", "Terrestrial")

all(sdms_CV$Species %in% n_occ$scientificName)

# N species
length(unique(sdms_CV$Species))

# add Taxonomic group and N occ
sdms_CV <- merge(sdms_CV, n_occ[,c("scientificName","Group","n_occ")], 
                 by.x = "Species", 
                 by.y = "scientificName",
                 all.x = TRUE)

plot(log(sdms_CV$n_occ),sdms_CV$validation)
plot(log(sdms_CV$n_occ),sdms_CV$calibration)
plot(sdms_CV$calibration,sdms_CV$validation)
plot(log(sdms_CV$n_occ),sdms_CV$sensitivity)

sdms_CV_plot <- sdms_CV[,c("algo","calibration","validation","Group","ECO")] %>%
    tidyr::gather("Type", "TSS", -c(algo,Group,ECO)) 
    
# Marine
ggplot(sdms_CV_plot %>%
           dplyr::filter(ECO=="Marine"), 
       aes(x = TSS, y = Group, color = algo))+
    geom_boxplot(outlier.alpha = 0)+
    facet_grid(ECO~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

ggplot(sdms_CV_plot %>%
           dplyr::filter(ECO=="Marine"), aes(x = TSS, y = algo, color = algo))+
    geom_violin(draw_quantiles = c(0.05,0.5,0.95)) +
    facet_grid(~Type,scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

# Terrestrial
ggplot(sdms_CV_plot %>%
           dplyr::filter(ECO=="Terrestrial"), 
       aes(x = TSS, y = Group, color = algo))+
    geom_boxplot(outlier.alpha = 0)+
    facet_grid(ECO~Type, scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

ggplot(sdms_CV_plot %>%
           dplyr::filter(ECO=="Terrestrial"), aes(x = TSS, y = algo, color = algo))+
    geom_violin(draw_quantiles = c(0.05,0.5,0.95)) +
    facet_grid(~Type,scales = "free")+
    geom_vline(xintercept = .5,color='red',linetype = "dashed")+
    geom_vline(xintercept = .7,color='darkgreen',linetype = "dashed")+
    theme_classic()

# total Marine species
length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Marine")]))
# total Marines species MAXNET TSS > 0.5
length(unique(sdms_CV$Species[which(sdms_CV$validation > 0.5 & sdms_CV$ECO=="Marine")]))/length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Marine")]))

# total Marines species MAXNET TSS > 0.7
length(unique(sdms_CV$Species[which(sdms_CV$validation > 0.7 & sdms_CV$ECO=="Marine")]))/length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Marine")]))

# total Terrestrial species
length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Terrestrial")]))
# total Terrestrial species MAXNET TSS > 0.5
length(unique(sdms_CV$Species[which(sdms_CV$validation > 0.5 & sdms_CV$ECO=="Terrestrial")]))/length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Terrestrial")]))

# total Terrestrial species MAXNET TSS > 0.7
length(unique(sdms_CV$Species[which(sdms_CV$validation > 0.7 & sdms_CV$ECO=="Terrestrial")]))/length(unique(sdms_CV$Species[which(sdms_CV$ECO=="Terrestrial")]))


###
### Load in all shifts SA ens ----

shifts_ens <- list.files(here::here(shift_dir,"Mar"), 
                         pattern = " ens shift SA.csv", 
                         recursive = TRUE, full.names = TRUE)

shifts_ens <- lapply(shifts_ens, function(x) {
    read.csv(x)
})

shifts_ens <- data.table::rbindlist(shifts_ens)
shifts_ens_mar <- data.frame(shifts_ens)
shifts_ens_mar$ECO = "Marine"

shifts_ens <- list.files(here::here(shift_dir,"Ter"), 
                         pattern = " ens shift SA.csv", 
                         recursive = TRUE, full.names = TRUE)

shifts_ens <- lapply(shifts_ens, function(x){
    read.csv(x)
})

shifts_ens <- data.table::rbindlist(shifts_ens)
shifts_ens_ter <- data.frame(shifts_ens)
shifts_ens_ter$ECO = "Terrestrial"

shifts_ens <- rbind(shifts_ens_mar,shifts_ens_ter)

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

shifts_ens <- pbapply::pblapply(unique(shifts_ens$new_ID), function(x){
    
    # subset shift i
    tmp <- shifts_ens %>% filter(new_ID == x)
    
    # get shifts at range positions
    SHIFT_sdm_1 <- c(mean(tmp$nsQuantVelocity_quant0p01),
                     mean(tmp$nsCentroidVelocity),
                     mean(tmp$nsQuantVelocity_quant0p99))
    
    SHIFT_sdm_2 <- c(mean(tmp$nsQuantVelocity_quant0p05),
                     mean(tmp$nsCentroidVelocity),
                     mean(tmp$nsQuantVelocity_quant0p95))
    
    SHIFT_sdm_3 <- c(mean(tmp$nsQuantVelocity_quant0p1),
                     mean(tmp$nsCentroidVelocity),
                     mean(tmp$nsQuantVelocity_quant0p9))
    
    tmp <- data.frame(Param = c("TE","O","LE"),
                      SHIFT_sdm_1 = SHIFT_sdm_1/1000,
                      SHIFT_sdm_2 = SHIFT_sdm_2/1000,
                      SHIFT_sdm_3 = SHIFT_sdm_3/1000,
                      tmp[1,c("Type", "Species", "time_period","new_ID","ECO","Group")])
    return(tmp)
})

shifts_ens <- data.table::rbindlist(shifts_ens)
shifts_ens <- data.frame(shifts_ens)

shifts_ens$new_ID <- paste(shifts_ens$new_ID, shifts_ens$Param)
shifts_ens <- na.omit(shifts_ens)
dim(shifts_ens)

# N species
length(unique(shifts_ens$Species))
# N shifts
length(unique(shifts_ens$new_ID))

table(shifts_ens$ECO, shifts_ens$Param)

### Correlation shifts % ----

## Marine
ggpairs(shifts_ens %>%
            filter(ECO == "Marine") %>%
            filter(Param != "O"), 
        columns = c(2:4),
        ggplot2::aes(colour=Param))

ggpairs(shifts_ens %>%
            filter(ECO == "Terrestrial") %>%
            filter(Param != "O"), 
        columns = c(2:4),
        ggplot2::aes(colour=Param))

###
### Merge with bioshifts ----
biov1_pred <- merge(
    biov1, shifts_ens[,c(2:4,8,10)], 
    by = "new_ID")

## Add TSS
avg_TSS <- sdms_CV %>%
    group_by(Species) %>%
    summarise(TSS_vali = mean(validation),
              TSS_cali = mean(calibration))

all(biov1_pred$sp_name_std_v1 %in% avg_TSS$Species)

biov1_pred <- merge(
    biov1_pred, avg_TSS, 
    by.x = "sp_name_std_v1",
    by.y = "Species",
    all.x = TRUE)

biov1_pred$SHIFT_obs <- biov1_pred$SHIFT
biov1_pred$SHIFT_pred <- biov1_pred$SHIFT_sdm_2

### Stats N species ----
# N species v1
length(unique(biov1$sp_name_std_v1))
# N species sdms
length(unique(shifts_ens$Species))
# N species with SDMs in bioshifts
length(which(unique(shifts_ens$Species) %in% unique(biov1$sp_name_std_v1)))

# shifts in v1
length(unique(biov1$new_ID))
# shifts in SDMs
length(unique(shifts_ens$new_ID))
# N shifts with SDMs in bioshifts
length(which(unique(shifts_ens$new_ID) %in% unique(biov1$new_ID)))

length(unique(shifts_ens$Species[which(shifts_ens$ECO == "Marine")]))
length(unique(shifts_ens$Species[which(shifts_ens$ECO == "Terrestrial")]))

spincommon <- intersect(biov1$sp_name_std_v1,shifts_ens$Species)

table(biov1$ECO, biov1$Param)
table(shifts_ens$ECO, shifts_ens$Param)
table(biov1_pred$ECO, biov1_pred$Param)

table(biov1_pred$ECO, biov1_pred$Param, biov1_pred$Group)


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

mundi <- terra::vect(rnaturalearth::ne_coastline(scale = 10, returnclass = "sp"))

# get raster bioshifts shp files study areas
get_raster_bioshifts = "NO"
if(get_raster_bioshifts=="YES"){
    my_ext = terra::ext(mundi)
    my_crs = crs(mundi)
    
    # empty raster for study areas
    rast_biov1 <- terra::rast(my_ext, crs = my_crs, res = 0.5)
    values(rast_biov1) <- 0
    
    # empty raster for N shifts
    rast_biov1_sh <- terra::rast(my_ext, crs = my_crs, res = 0.5)
    values(rast_biov1) <- 0
    
    fgdb <- "C:/Users/brunn/NextCloud/Bioshifts/Study_Areas.gdb"
    fc_list <- rgdal::ogrListLayers(fgdb)
    fc_list <- fc_list[which(fc_list %in% unique(biov1_pred$ID))]
    
    for(i in 1:length(fc_list)){ cat("\r",i,"from",length(fc_list))
        tmp = terra::vect(sf::st_read(fgdb, layer=fc_list[i]))
        tmp = terra::cells(rast_biov1,tmp)
        tmp_cell = tmp[,2]
        tmp_vals = rast_biov1[tmp_cell][,1]
        rast_biov1[tmp_cell] = tmp_vals+1
        
        # get N shifts at the study area
        N <- nrow(biov1_pred[which(biov1_pred$ID == fc_list[i]),])
        rast_biov1_sh[tmp_cell] = tmp_vals + N
    }
    names(rast_biov1) <- names(rast_biov1_sh) <- "SA"
    rast_biov1[rast_biov1==0] <- NA
    rast_biov1_sh[rast_biov1_sh==0] <- NA
    
    writeRaster(rast_biov1, "Data/raster_bioshifts_SA.tif", overwrite = TRUE)
    writeRaster(rast_biov1_sh, "Data/raster_bioshifts_N_SA.tif", overwrite = TRUE)
    
} else {
    rast_biov1 <- terra::rast("Data/raster_bioshifts_SA.tif")
    rast_biov1_sh <- terra::rast("Data/raster_bioshifts_N_SA.tif")
}

# plot
mundi_sf <- sf::st_as_sf(mundi)

test_df <- terra::as.data.frame(rast_biov1, xy = TRUE)
colnames(test_df) <- c("x", "y", "value")

p1 = ggplot() +  
    geom_tile(data=test_df, aes(x=x, y=y, fill=value)) + 
    ggtitle("N study areas")+
    geom_sf(data=mundi_sf) +
    theme_minimal() +
    scale_fill_viridis_c()+
    guides(fill = guide_colourbar(title = ""))+
    xlab("")+ylab("")

test_df <- terra::as.data.frame(rast_biov1_sh, xy = TRUE)
colnames(test_df) <- c("x", "y", "value")

p2 = ggplot() +  
    geom_tile(data=test_df, aes(x=x, y=y, fill=value)) + 
    ggtitle("N range shifts")+
    geom_sf(data=mundi_sf) +
    theme_minimal() +
    scale_fill_viridis_c()+
    guides(fill = guide_colourbar(title = ""))+
    xlab("")+ylab("")

gridExtra::grid.arrange(p1,p2,ncol=1)

###
### Shift predicted vc observed ----

ggplot(biov1_pred[,c("ECO","Param","SHIFT_obs","SHIFT_pred")] %>%
           tidyr::gather("Type", "SHIFT", -c(Param,ECO)), 
       aes(x = Type, y = SHIFT))+
    ylab("Range shift (km/yr)")+
    geom_violin(trim = TRUE, draw_quantiles = c(.05, .5, .95))+
    theme_classic()+
    facet_wrap(.~ECO, scales = "free")

ggplot(biov1_pred[,c("ECO","Param","SHIFT_obs","SHIFT_pred")] %>%
           tidyr::gather("Type", "SHIFT", -c(Param,ECO)), 
       aes(x = Type, y = SHIFT))+
    ylab("Range shift (km/yr)")+
    geom_violin(trim = TRUE, draw_quantiles = c(.05, .5, .95))+
    theme_classic()+
    facet_wrap(ECO~Param, scales = "free")

ggplot(biov1_pred[,c("ECO","Param","SHIFT_obs","SHIFT_pred")] %>%
           tidyr::gather("Type", "SHIFT", -c(Param,ECO))%>% 
           dplyr::filter(SHIFT < 100 & SHIFT > -30), 
       aes(fill = Type, x = SHIFT))+
    xlab("Range shift (km/yr)")+
    geom_density(alpha=.5)+
    theme_classic()+
    facet_wrap(.~ECO, scales = "free")

ggplot(biov1_pred[,c("ECO","Param","SHIFT_obs","SHIFT_pred")] %>%
           tidyr::gather("Type", "SHIFT", -c(Param,ECO)) %>% 
           dplyr::filter(SHIFT < 50 & SHIFT > -30), 
       aes(fill = Type, x = SHIFT))+
    xlab("Range shift (km/yr)")+
    geom_density(alpha=.5)+
    theme_classic()+
    facet_wrap(ECO~Param, scales = "free")

ggplot(
    biov1_pred %>%
        filter(SHIFT_obs < 50), 
    aes(x = SHIFT_pred, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("All shifts")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

# Optimum
M_O <- ggplot(
    biov1_pred %>%
        filter(Param == "O" & SHIFT_obs < 50), 
    aes(x = SHIFT_pred, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Centroid")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

# LE
M_LE <- ggplot(
    biov1_pred %>%
        filter(Param == "LE" & SHIFT_obs < 50), 
    aes(x = SHIFT_pred, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

# TE
M_TE <- ggplot(
    biov1_pred %>%
        filter(Param == "TE" & SHIFT_pred < 50 & SHIFT_pred > -20), 
    aes(x = SHIFT_pred, y = SHIFT_obs, color = ECO, size = DUR))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()


M_O
M_LE
M_TE

# Marines ----


## % species shifting same direction of bioclimatic velocity
tmp <- biov1_pred %>%
    filter(ECO == "Marine")
           
table(sign(tmp$SHIFT_obs) == sign(tmp$SHIFT_pred), tmp$Param)

# Optimum
M_O <- ggplot(
    biov1_pred %>%
        filter(ECO == "Marine" & Param == "O"), 
    aes(x = SHIFT_pred, y = SHIFT_obs, color = Group, size = DUR))+
    ggtitle("Centroid")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

# LE
M_LE <- ggplot(
    biov1_pred %>%
        filter(ECO == "Marine" & Param == "LE" & SHIFT_obs < 25), 
    aes(x = SHIFT_pred, y = SHIFT_obs, color = Group, size = DUR))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

# TE
M_TE <- ggplot(
    biov1_pred %>%
        filter(ECO == "Marine" & Param == "TE" & SHIFT_pred>-20), 
    aes(x = SHIFT_pred, y = SHIFT_obs, color = Group, size = DUR))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()


M_O
M_LE
M_TE

gridExtra::grid.arrange(grobs = list(M_O,
                                     M_LE,
                                     M_TE),ncol = 3)

# just class with > 10 shifts
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10 & tmp$Var2=="Marine",]
tmp <- tmp[order(tmp$Var1),]

plots <- list()
for(i in 1:length(tmp$Var1)){
    
    plots[[i]] <- ggplot(
        biov1_pred %>%
            filter(ECO == "Marine" & Param == tmp$Var3[i] & Group == tmp$Var1[i]), 
        aes(x = SHIFT_pred, y = SHIFT_obs, color = Class))+
        ggtitle(paste(tmp$Var1[i],tmp$Var3[i]))+
        geom_point(alpha = .5)+
        xlab("Bioclimatic velocity (km/yr)")+
        ylab("Observed shift (km/yr)")+
        theme_classic() +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
        geom_hline(yintercept = 0, linetype = "dashed")+
        geom_vline(xintercept = 0, linetype = "dashed")+
        tune::coord_obs_pred()
}
names(plots) <- paste(tmp$Var1,tmp$Var3)

gridExtra::grid.arrange(grobs = plots,ncol = 2)



# Terrestrials

tmp <- biov1_pred %>%
    filter(ECO == "Terrestrial")

table(sign(tmp$SHIFT_obs) == sign(tmp$SHIFT_pred), tmp$Param)

# Optimum
T_O <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial" & Param == "O" & SHIFT_obs < 100 & SHIFT_obs > -20), 
    aes(x = SHIFT_pred, y = SHIFT_obs, color = Group))+
    ggtitle("Centroid")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

# LE
T_LE <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial" & Param == "LE"), 
    aes(x = SHIFT_sdm_3, y = SHIFT_obs, color = Group))+
    ggtitle("Leading edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

# TE
T_TE <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial" & Param == "TE"), 
    aes(x = SHIFT_pred, y = SHIFT_obs, color = Group))+
    ggtitle("Trailing edge")+
    geom_point(alpha = .3)+
    xlab("Bioclimatic velocity (km/yr)")+
    ylab("Observed shift (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()


T_O
T_LE
T_TE

gridExtra::grid.arrange(grobs = list(T_O,T_LE,T_TE),ncol = 3)

# just class with > 10 shifts
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10 & tmp$Var2=="Terrestrial",]
tmp <- tmp[order(tmp$Var1),]


plots <- list()
for(i in 1:length(tmp$Var1)){
    
    plots[[i]] <- ggplot(
        biov1_pred %>%
            filter(ECO == "Terrestrial" & Param == tmp$Var3[i] & Group == tmp$Var1[i]), 
        aes(x = SHIFT_pred, y = SHIFT_obs, color = Class))+
        ggtitle(paste(tmp$Var1[i],tmp$Var3[i]))+
        geom_point(alpha = .8)+
        xlab("Bioclimatic velocity (km/yr)")+
        ylab("Observed shift (km/yr)")+
        theme_classic() +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
        geom_hline(yintercept = 0, linetype = "dashed")+
        geom_vline(xintercept = 0, linetype = "dashed")+
        tune::coord_obs_pred()
}
names(plots) <- paste(tmp$Var1,tmp$Var3)

gridExtra::grid.arrange(grobs = plots,ncol = 4)



###
### calculate lags ----
biov1_pred$biolag1 <- biov1_pred$SHIFT_pred - biov1_pred$SHIFT_obs

# fix lag
# same direction: change lag sign when vel is negative
pos <- which(biov1_pred$SHIFT_pred<0 & biov1_pred$SHIFT_obs<0)
biov1_pred$biolag1[pos] <- biov1_pred$biolag1[pos] * sign(biov1_pred$SHIFT_pred[pos])

# opposite direction: lag = vel
biov1_pred$biolag2 <- biov1_pred$biolag1
pos <- which((biov1_pred$SHIFT_pred>0 & biov1_pred$SHIFT_obs<0)|
                 biov1_pred$SHIFT_pred<0 & biov1_pred$SHIFT_obs>0)
biov1_pred$biolag2[pos] <- biov1_pred$SHIFT_pred[pos]


## % species shifting faster vs slower than expected
table(biov1_pred$biolag1 > 0, biov1_pred$ECO)


## Marines
# just class with > 10 shifts
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10 & tmp$Var2=="Marine",]

limits <- biov1_pred %>%
    dplyr::filter(ECO == "Marine" & Group %in% tmp$Var1) %>%
    summarise(quantile(biolag1,c(0.05, 0.95))) 

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Marine" & Group %in% tmp$Var1), 
       aes(x = biolag1, y = Group, color = Param))+
    # geom_violin(trim = TRUE, draw_quantiles = c(0.05,0.5,0.95)) +
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = limits[,1])+
    xlab("Lag = Bioclimatic velocity - Range shift velocity")+
    ylab("")+
    geom_vline(xintercept = 0, linetype = "longdash")+
    facet_grid(Param~., scales = "free", space = "free")+
    theme_classic()

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Marine" & Group %in% tmp$Var1), 
       aes(x = biolag1, fill = Group, color = Group))+
    # geom_violin(trim = TRUE, draw_quantiles = c(0.05,0.5,0.95)) +
    geom_density(alpha = 0.5)+
    scale_x_continuous(limits = limits[,1])+
    xlab("Lag = Bioclimatic velocity - Range shift velocity")+
    ylab("Density")+
    geom_vline(xintercept = 0, linetype = "longdash")+
    facet_wrap(~Param, scales = "free")+
    theme_classic()

## Terrestrials
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10 & tmp$Var2=="Terrestrial",]

limits <- biov1_pred %>%
    dplyr::filter(ECO == "Marine" & Group %in% tmp$Var1) %>%
    summarise(quantile(biolag1,c(0.05, 0.95))) 

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & Group %in% tmp$Var1), 
       aes(x = biolag1, y = Group, color = Param))+
    # geom_violin(trim = TRUE, draw_quantiles = c(0.05,0.5,0.95)) +
    geom_boxplot(outlier.alpha = 0)+
    scale_x_continuous(limits = limits[,1])+
    xlab("Lag = Bioclimatic velocity - Range shift velocity")+
    ylab("")+
    geom_vline(xintercept = 0, linetype = "longdash")+
    facet_grid(Param~., scales = "free", space = "free")+
    theme_classic()

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & Group %in% tmp$Var1), 
       aes(x = biolag1, fill = Group, color = Group))+
    # geom_violin(trim = TRUE, draw_quantiles = c(0.05,0.5,0.95)) +
    geom_density(alpha = 0.3)+
    scale_x_continuous(limits = limits[,1])+
    xlab("Lag = Bioclimatic velocity - Range shift velocity")+
    ylab("")+
    geom_vline(xintercept = 0, linetype = "longdash")+
    facet_wrap(~Param, scales = "free")+
    theme_classic()

###
### Lag vs duration ----

### Marine

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Marine"),
       aes(x = DUR, y = biolag1, size=abs(SHIFT_obs)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("Duration (years)")+
    theme_bw()+
    facet_wrap(.~Param)

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Marine"),
       aes(x = START, y = biolag1, size=abs(SHIFT_obs)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("Start year")+
    theme_bw()+
    facet_wrap(.~Param)

### Terrestrial

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & biolag1>-100),
       aes(x = DUR, y = biolag1, size=abs(SHIFT_obs)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("Duration (years)")+
    theme_bw()+
    facet_wrap(.~Param)

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & biolag1>-100),
       aes(x = START, y = biolag1, size=abs(SHIFT_obs)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("Start year")+
    theme_bw()+
    facet_wrap(.~Param)



### Lag vs TSS ----

### Marine

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Marine"),
       aes(x = TSS_vali, y = abs(biolag1),size=abs(SHIFT_pred)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("TSS validation")+
    theme_bw()+
    facet_wrap(.~Param)

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Marine"),
       aes(x = TSS_cali, y = abs(biolag1),size=abs(SHIFT_pred)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("TSS calibration")+
    theme_bw()+
    facet_wrap(.~Param)

### Terrestrial

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & biolag1>-100),
       aes(x = TSS_vali, y = abs(biolag1),size=abs(SHIFT_pred)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("TSS validation")+
    theme_bw()+
    facet_wrap(.~Param, scales = "free")

ggplot(biov1_pred %>%
           dplyr::filter(ECO == "Terrestrial" & biolag1>-100),
       aes(x = TSS_cali, y = abs(biolag1),size=abs(SHIFT_pred)))+
    geom_point(alpha=.3)+
    geom_smooth(method = 'lm')+
    ylab("Lag = Bioclimatic velocity - Range shift velocity")+
    xlab("TSS calibration")+
    theme_bw()+
    facet_wrap(.~Param, scales = "free")


###
### Predicted shifts vs Temperature velocities ----

## Does bioclimatic velocity explains shifts better then temperature velocity?

### Marines
# all groups with > > 10 shifts
tmp <- data.frame(table(biov1_pred$Group,biov1_pred$ECO, biov1_pred$Param))
tmp <- tmp[tmp$Freq>10,]

mods <- data.frame()

for(i in 1:nrow(tmp)){
    
    mod_vel_i <- lm(
        SHIFT_obs~v.lat.mean,
        data = biov1_pred %>%
            filter(ECO == tmp$Var2[i] & 
                       Param == tmp$Var3[i] & 
                       Group == tmp$Var1[i]))
    
    mod_vel_i_sum <- summary(mod_vel_i)
    
    mod_bio_i <- lm(
        SHIFT_obs~SHIFT_pred,
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

## Marine
p1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Marine"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_1, color = Param))+
    ggtitle("TE = 1%, LE = 99%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

p2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Marine"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_2, color = Param))+
    ggtitle("TE = 5%, LE = 95%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

p3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Marine"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_3, color = Param))+
    ggtitle("TE = 10%, LE = 90%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

gridExtra::grid.arrange(p1,p2,p3,ncol=3,top="Marine")



## Terrestrial
p1 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_1, color = Param))+
    ggtitle("TE = 1%, LE = 99%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

p2 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_2, color = Param))+
    ggtitle("TE = 5%, LE = 95%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

p3 <- ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial"), 
    aes(x = v.lat.mean, y = SHIFT_sdm_3, color = Param))+
    ggtitle("TE = 10%, LE = 90%")+
    geom_point(alpha = .3)+
    xlab("Temperature velocity (km/yr)")+
    ylab("Bioclimatic velocity (km/yr)")+
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    tune::coord_obs_pred()

gridExtra::grid.arrange(p1,p2,p3,ncol=3,top="Terrestrials")


ggplot(
    biov1_pred %>%
        filter(ECO == "Terrestrial", Param == "O"), 
    aes(x = velocity, y = SHIFT_pred))+
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

# 1) better predictions for marine than terrestrials because landscape is more connected - dispersal constrains in land. marine more in eq with the environment
# 2) but we would expect lags associated with dispersal and persistance (demography, life history)

# sdms methods

# bioshifts

# plot all param >>> only LE, TE O

# expectations by group >> fishes should track better than crustacea

## Plot lag by group

## next steps
# account for methodological factors
# connectivity
# microclimatic variability
# traits

# funny pic
# announcement about interaction with the group Wednesday after to movie
