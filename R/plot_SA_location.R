plot_SA_location <- function(sptogo, 
                             BA,
                             realm, 
                             PresAbs,
                             StudyArea,
                             maps.dir = "/media/seagate/boliveira/SDMs/Maps_SA",
                             show.plot = TRUE){
    
    maps.dir <- here::here(maps.dir)
    
    if(!dir.exists(maps.dir)){
        dir.create(maps.dir,recursive = T)
    }
    
    sp_data <- PresAbs[which(PresAbs$pa==1),]
    xyk <- terra::vect(sp_data[,c('x','y')], geom = c('x','y'))
    crs(xyk) <- crs(BA)
    
    map_back <- vect(rnaturalearth::ne_coastline(scale = 110, returnclass = "sf"))
    
    back_back <- lapply(1:length(StudyArea), function(x) {
        tmp = as.polygons(terra::ext(StudyArea[x])+5)
        crs(tmp) <- crs(BA)
        return(tmp)
    })
    back_back <- vect(back_back)
    
    if(!show.plot){
        pdf(here::here(maps.dir,paste0(paste(sptogo,realm, sep = "_"),".pdf")),
            width = 10, height = 10)
    }
    
    par(mfrow = c(2,1))
    plot(map_back, main = gsub("_"," ",sptogo), lwd = .5, col = "darkgray", type = "none")
    plot(BA, alpha = .3, add=T, legend=FALSE)
    plot(back_back, border = NULL, alpha = .5, add=T, type = "none")
    plot(xyk, col = "red", cex = .3, border = NULL, alpha = .3, add=T, type = "none")
    plot(StudyArea, col = "blue", border = NULL, alpha = .3, add=T, type = "none")
    
    plot(back_back, main = "Study area", border = NULL, alpha = .5, type = "none")
    plot(map_back, lwd = .5, col = "darkgray", alpha = .5, add=T, type = "none")
    plot(BA, alpha = .3, add=T, legend=FALSE)
    plot(xyk, col = "red", border = NULL, cex = 1, alpha = .3, add=T, type = "none")
    plot(StudyArea, col = "blue", border = NULL, alpha = .3, add=T, type = "none")
    
    if(!show.plot){
        dev.off()
    }
}