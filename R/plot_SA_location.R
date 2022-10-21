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
    
    pdf(here::here(maps.dir,paste0(paste(sptogo,realm, sep = "_"),".pdf")),
        width = 10, height = 10)
    
    par(mfrow = c(2,1))
    
    sp_data <- PresAbs[which(PresAbs$pa==1),]
    xyk <- terra::vect(SpatialPoints(data.frame(sp_data[,c('x','y')])))
    crs(xyk) <- crs(BA)
    
    map_back <- vect(rnaturalearth::ne_coastline(scale = 110, returnclass = "sp"))
    
    back_back <- lapply(1:length(StudyArea), function(x) {
        tmp = as.polygons(terra::ext(StudyArea[x])+5)
        crs(tmp) <- crs(BA)
        return(tmp)
    })
    back_back <- vect(back_back)
    
    plot(map_back, main = gsub("_"," ",sptogo), lwd = .5, col = "darkgray")
    plot(BA, alpha = .3, add=T, legend=FALSE)
    plot(back_back, border = NULL, alpha = .5, add=T)
    plot(xyk, col = "red", cex = .3, border = NULL, alpha = .3, add=T)
    plot(StudyArea, col = "blue", border = NULL, alpha = .3, add=T)
    
    plot(back_back, main = "Study area", border = NULL, alpha = .5)
    plot(map_back, lwd = .5, col = "darkgray", alpha = .5, add=T)
    plot(BA, alpha = .3, add=T, legend=FALSE)
    plot(xyk, col = "red", border = NULL, cex = 1, alpha = .3, add=T)
    plot(StudyArea, col = "blue", border = NULL, alpha = .3, add=T)
    
    dev.off()
    
    if(show.plot){
        par(mfrow = c(2,1))
        
        sp_data <- PresAbs[which(PresAbs$pa==1),]
        xyk <- terra::vect(SpatialPoints(data.frame(sp_data[,c('x','y')])))
        crs(xyk) <- crs(BA)
        
        map_back <- vect(rnaturalearth::ne_coastline(scale = 110, returnclass = "sp"))
        
        back_back <- lapply(1:length(StudyArea), function(x) {
            tmp = as.polygons(terra::ext(StudyArea[x])+5)
            crs(tmp) <- crs(BA)
            return(tmp)
        })
        back_back <- vect(back_back)
        
        plot(map_back, main = gsub("_"," ",sptogo), lwd = .5, col = "darkgray")
        plot(BA, alpha = .3, add=T, legend=FALSE)
        plot(back_back, border = NULL, alpha = .5, add=T)
        plot(xyk, col = "red", cex = .3, border = NULL, alpha = .3, add=T)
        plot(StudyArea, col = "blue", border = NULL, alpha = .3, add=T)
        
        plot(back_back, main = "Study area", border = NULL, alpha = .5)
        plot(map_back, lwd = .5, col = "darkgray", alpha = .5, add=T)
        plot(BA, alpha = .3, add=T, legend=FALSE)
        plot(xyk, col = "red", border = NULL, cex = 1, alpha = .3, add=T)
        plot(StudyArea, col = "blue", border = NULL, alpha = .3, add=T)
    }
}
