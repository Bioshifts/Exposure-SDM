library(raster)
library(terra)
library(enmSdmX)


# source funcion
source("R/bioshiftsFunction.R")

# simulate a shift in species distribution:

rbase <- raster(extent(c(0,1,0,1)), nrow = 50, ncol = 50)
rx1 <- rx2 <- rbase
rx1[] <- xFromCell(rbase, 1:ncell(rbase))
rx2[] <- yFromCell(rbase, 1:ncell(rbase))

## t1:
rt1 <- 2*rx1  - 2*rx1^2- 3*rx2
rt1 <- range01raster(rast(rt1))

## t2:
rt2 <- 2*rx1  - 2*rx1^3 - rx2/2
rt2 <- range01raster(rast(rt2))

## view
par(mfrow=c(1,2))

plot(rt1, col = the_colors)
plot(rt2, col = the_colors)

####
## 1) simulate a shift Northwards

r1 = rt1
r2 = rt2
shift_plot(r1, r2)

## 2) simulate a shift Southwards

r1 = rt2
r2 = rt1
shift_plot(r1, r2)

## 3) simulate a shift Northwards in the South Hemisphere

r1 = rt2
r2 = rt1

r1 <- as.data.frame(r1, xy = TRUE)
r1$y <- r1$y * -1
r1 <- rast(r1, type = c("xyz"))

r2 <- as.data.frame(r2, xy = TRUE)
r2$y <- r2$y * -1
r2 <- rast(r2, type = c("xyz"))

shift_plot(r1, r2, times = c(1990,1999))

## 4) simulate a shift Southwards in the South Hemisphere

r1 = rt1
r2 = rt2

r1 <- as.data.frame(r1, xy = TRUE)
r1$y <- r1$y * -1
r1 <- rast(r1, type = c("xyz"))

r2 <- as.data.frame(r2, xy = TRUE)
r2$y <- r2$y * -1
r2 <- rast(r2, type = c("xyz"))

shift_plot(r1, r2)
