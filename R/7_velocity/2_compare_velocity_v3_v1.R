# Compare temperature velocity values in Bioshifts v1 with the new velocities calculated for Bioshifts v3

### P.S.: If running this on the cluster:
### Use library rmote for visualization of graphs
### ssh -L 4321:localhost:4321 -L 8100:localhost:8100 oliveirab@muse-login02.meso.umontpellier.fr 
### On terminal, open new screen, load R
### load package rmote and run rmote::start_rmote() !!!

library(data.table)
library(dplyr)

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    
    work_dir <- getwd()
}

# source settings
source("R/settings.R")

# Load in v3 velocities
SAs <- list.files(velocity_SA_dir)
SAs <- gsub(".csv","",SAs)
SA_got <- list.files(velocity_SA_dir, full.names = TRUE)
SA_got <- lapply(1:length(SA_got), function(i) {
    tmp <- read.csv(SA_got[i])
    tmp$ID = SAs[i]
    tmp
})
SA_got <- rbindlist(SA_got, fill = TRUE)
names(SA_got)
SA_got <- SA_got %>% dplyr::select(c(mat_gVel_Median,mat_gVelLat_Median,Type_grad,ID))
SA_got$vel_v3 <- SA_got$mat_gVelLat_Median
SA_got$vel_v3[which(SA_got$Type_grad == "ELE")] <- SA_got$mat_gVel_Median[which(SA_got$Type_grad == "ELE")]
SA_got <- SA_got %>% dplyr::select(c(vel_v3,Type_grad,ID))

# load in Bioshifts v1
Bioshifts_DB <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v1))

# Filter Polygons in Study areas v3
Bioshifts_DB <- Bioshifts_DB[Bioshifts_DB$ID %in% unique(SA_got$ID),]

# Latitude shifts
Bioshifts_DB_LAT <-merge(Bioshifts_DB %>% filter(Type == "LAT"),
                         SA_got %>% filter(Type_grad == "LAT"),
                         by = "ID")
dim(Bioshifts_DB_LAT)
names(Bioshifts_DB_LAT)
unique(Bioshifts_DB_LAT$UNIT)

outliers <- summary(na.omit(Bioshifts_DB_LAT$vel_v3))[c(2,5)]

plot(vel_v3~v.lat.mean,
     Bioshifts_DB_LAT %>% 
         filter(vel_v3 > outliers[1] & vel_v3 < 4000 & !is.infinite(vel_v3)),
     xlab = "Velocity v1",
     ylab = "Velocity v3")
rmote::plot_done()

par(mfrow=c(1,2))
plot(SHIFT~v.lat.mean,
     Bioshifts_DB_LAT %>% 
         filter(vel_v3 > outliers[1] & vel_v3 < 4000 & !is.infinite(vel_v3)),
     xlab = "Velocity v1",
     ylab = "SHIFT v1")
# rmote::plot_done()

plot(SHIFT~vel_v3,
     Bioshifts_DB_LAT %>% 
         filter(vel_v3 > outliers[1] & vel_v3 < 4000 & !is.infinite(vel_v3)),
     xlab = "Velocity v3",
     ylab = "SHIFT v1")
rmote::plot_done()


par(mfrow=c(1,2))
tmp <- Bioshifts_DB_LAT %>% 
    filter(vel_v3 > outliers[1] & vel_v3 < 4000 & !is.infinite(vel_v3))
hist(tmp$vel_v3)
hist(tmp$v.lat.mean)
dev.off()


# Elevation shifts
Bioshifts_DB_ELE <- merge(Bioshifts_DB %>% filter(Type == "ELE"),
                          SA_got %>% filter(Type_grad == "ELE"),
                          by = "ID")
names(Bioshifts_DB_ELE)
unique(Bioshifts_DB_ELE$UNIT)

plot(mat_gVel_Mean~v.ele.mean,
     Bioshifts_DB_ELE,
     xlab = "Velocity v1",
     ylab = "Velocity v3")
rmote::plot_done()

plot(SHIFT~v.ele.mean,
     Bioshifts_DB_ELE,
     xlab = "Velocity v1",
     ylab = "SHIFT v1")
rmote::plot_done()

plot(SHIFT~mat_gVel_Mean,
     Bioshifts_DB_ELE,
     xlab = "Velocity v3",
     ylab = "SHIFT v1")
rmote::plot_done()
