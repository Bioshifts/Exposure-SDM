plot_SA_location <- function(sptogo, 
                             BA_shp,
                             realm, 
                             sp_occ,
                             StudyArea,
                             maps.dir = NULL,
                             show.plot = TRUE){
    
    maps.dir <- here::here(maps.dir)
    
    if(!dir.exists(maps.dir)){
        dir.create(maps.dir,recursive = T)
    }
    
    sp_data <- sp_occ[which(sp_occ$pa==1),]
    sp_occ <- terra::vect(sp_data[,c('x','y')], geom = c('x','y'))
    crs(sp_occ) <- crs(BA_shp)
    
    map_globe <- vect(rnaturalearth::ne_coastline(scale = 110, returnclass = "sf"))
    
    SA_square <- lapply(1:length(StudyArea), function(x) {
        tmp = as.polygons(terra::ext(StudyArea[x])+5)
        crs(tmp) <- crs(BA_shp)
        return(tmp)
    })
    SA_square <- vect(SA_square)
    
    if(!show.plot){
        pdf(here::here(maps.dir,paste0(paste(sptogo,realm, sep = "_"),".pdf")),
            width = 10, height = 10)
    }
    
    par(mfrow = c(2,1))
    plot(map_globe, main = gsub("_"," ",sptogo), lwd = .5, col = "darkgray", type = "none", axes = FALSE)
    plot(BA_shp, alpha = .3, add=T, legend=FALSE, col = "green", border = NULL)
    plot(SA_square, border = NULL, alpha = .5, add=T, type = "none")
    plot(sp_occ, col = "red", cex = .3, border = NULL, alpha = .3, add=T, type = "none")
    plot(StudyArea, col = "blue", border = NULL, alpha = .3, add=T, type = "none")
    
    plot(SA_square, main = "Study area", border = NULL, alpha = .5, type = "none", box = FALSE)
    plot(map_globe, lwd = .5, col = "darkgray", alpha = .5, add=T, type = "none")
    plot(BA_shp, alpha = .3, add=T, legend=FALSE, col = "green", border = NULL)
    plot(StudyArea, col = "blue", border = NULL, alpha = .3, add=T, type = "none")
    plot(sp_occ, col = "red", border = NULL, cex = 1, alpha = .3, add=T, type = "none")
    
    if(!show.plot){
        dev.off()
    }
}