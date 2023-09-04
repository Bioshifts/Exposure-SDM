library(raster)
library(rgdal)
library(maptools)
library(fields)
library(RColorBrewer)


#==============================
# Spatial gradient
# WORLDCLIM - 10 min. 1950-2000
#==============================
rast=raster("bio1_5m.bil")  # Annual Mean Temperature
image(rast, col = rev(brewer.pal(10, "RdBu")))

resol=0.08333333
e1=extent(-180-resol,180+resol,-60-resol,90+resol)      
e=extent(-180-2*resol,180+2*resol,-60-2*resol,90+2*resol)

r_south <- extend(shift(rast, x=0, y=-resol),e1)
r_north <- extend(shift(rast, x=0, y=+resol),e1)
r_NS=r_south-r_north

delta_ag=extend(shift(r_NS,x=resol,y=0),e)
delta_bh=extend(r_NS,e)
delta_ci=extend(shift(r_NS,x=-resol,y=0),e)

moy_lat=(delta_ag+2*delta_bh+delta_ci)/8     # 10°C/pixel
image(moy_lat, col = rev(brewer.pal(10, "RdBu")))
writeRaster(moy_lat,"grad_latitude",overwrite=TRUE)

liste=list.files("C:/Users/Lise/AppData/Local/Temp/R_raster_Lise/")
file.remove(paste("C:/Users/Lise/AppData/Local/Temp/R_raster_Lise/", liste, sep="/"))

#===================
# Climate Velocities
#===================
lapse.rate=raster("grad_latitude.grd")      # 10°C/pixel
lapse.rate.km = lapse.rate / (10*9.2766)    # °C/km
temporal.trend=raster("temporal_trend.grd")	# °C/year
temporal.trend = resample(temporal.trend, lapse.rate)

velo = temporal.trend / lapse.rate.km       # km/year

# Correcting velocities in North Hemisphere
e2=extent(-180.5, 180.5, 0, 90.58333)
velo.nord = crop(velo, e2)
velo.nord = -1*velo.nord
e3=extent(-180.5, 180.5, -60.58333, 0)
velo.sud = crop(velo, e3)
vel = merge(velo.nord, velo.sud)

writeRaster(vel,"climate_velocity_new",overwrite=TRUE)	

#vel.fig = vel
#vel.fig[vel.fig>10] = 10
#vel.fig[vel.fig< -10] = -10
#image(vel.fig, col = rev(brewer.pal(10, "RdBu")), zlim=c(-10,10))
#my_breaks=c(-10,-7,-4,-1,0,0.01,1,4,7,10)
#my_col= rev(brewer.pal(10, "RdBu"))
#plot(vel.fig, breaks = my_breaks, col = my_col)

#image(vel.fig, col = my_col, zlim=c(-10,10))

#=============================
# Extract data from shapefiles
#=============================
setwd("C:/Users/Lise/Desktop/SOTM 2016/ANALYSES/SHP_ALL")
list_shp <- list.files(path=".", pattern="*.shp", full.names = TRUE, recursive = TRUE, include.dirs = FALSE)
shp_objects <- lapply(list_shp, function(x) {readOGR(dsn=x, layer=ogrListLayers(x))})

setwd("C:/Users/Lise/Desktop/SOTM 2016/ANALYSES/CLIMATE VELOCITY") 
vel = raster("climate_velocity_new.grd")
temporal.trend = raster("temporal_trend.grd") 
vel.shape = data.frame(matrix(nrow=length(shp_objects), ncol=4))
colnames(vel.shape) = c("Article_Study_ID", "vel.lat.median", "vel.ele.median", "trend")

for(i in 1:length(shp_objects)){
shape = shp_objects[[i]]
v.lat = crop(vel, shape)
trend = crop(temporal.trend, shape)
vel.shape[i,1] = strsplit(strsplit(list_shp[i], split = "[.]")[[1]][2],"/")[[1]][2]
vel.shape[i,2] = median(v.lat, na.rm=T)
vel.shape[i,4] = median(trend, na.rm=T)
vel.shape[i,3] = vel.shape[i,4]/(6.5/1000)
print(i)}

write.table(vel.shape, "velocities.shapes.txt", sep="\t")
