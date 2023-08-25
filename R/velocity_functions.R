# Authors: Brunno Oliveira
# Last update :  Oct 2023
#--------
# This is based on the functions from the climetrics package (https://github.com/shirintaheri/climetrics)

# temporal gradients using time series:
# threshold of minimun N obs are considered to get slope 
temp_gradFun <- function(x, th) {
    x <- x[!is.na(x)]
    if (length(x) > th) {
        s <- lm(x~c(1:length(x)))
        s$coefficients[2]
    } else NA
}
#--------------
temp_grad <- function(x, th) {
    tmp <- app(x,temp_gradFun,th=th)
    names(tmp) <- "Trend"
    return(tmp)
}

#--------------

spatial_grad <- function(rx, y_diff = 1) {
    
    if(nlyr(rx) > 1){ rx <- mean(rx,na.rm = TRUE) }
    
    if (.getProj(rx) == 'longlat') y_dist <- res(rx) * c(111.325, 111.325)
    
    if (!.is_package_installed("dplyr") || !.is_package_installed('tidyr')) stop('The packages dplyr and tidyr are needed for this metric; Please make sure they are installed!')
    
    else {
        y_dist <- res(rx)
        y_diff <- NA
    }
    
    nlats <- nrow(rx)
    nlons <- ncol(rx)
    y <- data.frame(adjacent(rx, cells=1:ncell(rx), directions=8,pairs=TRUE))
    y <- y[order(y$from, y$to),]
    y <- na.omit(y)
    y$sst <- rx[y$to][,1]
    y$sy <- rowFromCell(rx, y$from)-rowFromCell(rx, y$to)
    y$sx <- colFromCell(rx, y$to)-colFromCell(rx, y$from)
    y$sx[y$sx > 1] <- -1
    y$sx[y$sx < -1] <- 1
    y$code <- paste(y$sx, y$sy)
    
    y$code1 <- eval(parse(text='dplyr::recode(y$code,
                           `1 0` = "sstE",
                           `-1 0` = "sstW",
                           `-1 1` = "sstNW",
                           `-1 -1` = "sstSW",
                           `1 1` = "sstNE",
                           `1 -1` = "sstSE",
                           `0 1` = "sstN",
                           `0 -1` = "sstS")'),envir =environment())
    # y$code1 <- dplyr::recode(y$code,
    #                          `1 0` = "sstE",
    #                          `-1 0` = "sstW",
    #                          `-1 1` = "sstNW",
    #                          `-1 -1` = "sstSW",
    #                          `1 1` = "sstNE",
    #                          `1 -1` = "sstSE",
    #                          `0 1` = "sstN",
    #                          `0 -1` = "sstS")
    
    y3b <- eval(parse(text="dplyr::select(y,from, code1, sst)"),envir =environment())
    y3b <- eval(parse(text="tidyr::spread(y3b,code1, sst)"),envir =environment())
    y3b$sstFocal <- rx[y3b$from][,1]
    y3b$LAT <- yFromCell(rx, y3b$from)
    
    if(!is.na(y_diff)) {
        y3b <- eval(parse(text="dplyr::mutate(y3b,
                         latpos = cos(.rad(LAT + y_diff)),
                         latneg = cos(.rad(LAT - y_diff)),
                         latfocal = cos(.rad(LAT)))"),envir =environment())
    } else {
        
        y3b <- eval(parse(text="dplyr::mutate(y3b,
                         latpos = 1,
                         latneg = 1,
                         latfocal = 1)"),envir =environment())
    }
    
    y3c <- "dplyr::mutate(y3b,
                       gradWE1 = (sstN-sstNW)/
                         (latpos *  y_dist[1]),
                       gradWE2 = (sstFocal - sstW)/(latfocal * y_dist[1]),
                       gradWE3 = (sstS-sstSW)/(latneg * y_dist[1]),
                       gradWE4 = (sstNE-sstN)/(latpos * y_dist[1]),
                       gradWE5 = (sstE-sstFocal)/(latfocal * y_dist[1]),
                       gradWE6 = (sstSE-sstS)/(latneg*y_dist[1]),
                       gradNS1 = (sstNW-sstW)/y_dist[2],
                       gradNS2 = (sstN-sstFocal)/y_dist[2],
                       gradNS3 = (sstNE-sstE)/y_dist[2],
                       gradNS4 = (sstW-sstSW)/y_dist[2],
                       gradNS5 = (sstFocal-sstS)/y_dist[2],
                       gradNS6 = (sstE-sstSE)/y_dist[2])" 
    
    y3c <- eval(parse(text=y3c),envir=environment())
    
    y3c <- eval(parse(text="dplyr::rowwise(y3c)"),envir=environment())
    y3c <- eval(parse(text="dplyr::mutate(y3c,
      WEgrad = .mnwm(gradWE1, gradWE2, gradWE3, gradWE4, gradWE5, gradWE6),
      NSgrad = .mnwm(gradNS1, gradNS2, gradNS3, gradNS4, gradNS5, gradNS6),
      angle = .ang(WEgrad, NSgrad))"),envir=environment())
    
    y3c <- eval(parse(text="dplyr::select(y3c,icell = from, WE = WEgrad, NS = NSgrad, angle = angle)"),envir=environment())
    
    NS <- y3c$NS
    WE <- y3c$WE
    NS[is.na(NS)] <- 0
    WE[is.na(WE)] <- 0
    NAsort <- ifelse((abs(NS)+abs(WE)) == 0, NA, 1)
    y3c$Grad <- NAsort * sqrt((WE^2) + (NS^2))
    
    return(y3c)
}

#----------
# squeeze is for bouding max and min values to upper (95%) and lower (5%) quantiles, respectively
gVelocity <- function(grad, slope, squeeze=FALSE) {
    
    v <- rast(slope)
    
    v[grad$icell] <- slope[grad$icell] / grad$Grad
    
    grad$slope <- NA
    grad$slope[grad$icell] <- as.numeric(slope[grad$icell])
    
    if(squeeze){
        .o <- as.matrix(global(v,fun=quantile,probs=c(0.05,0.95),na.rm=TRUE))[1,]
        v[v < .o[1]] <- .o[1]
        v[v > .o[2]] <- .o[2] 
    }
    
    names(v) <- "GradVel"
    return(list(GradVel = v,
                grad))
}


###########################
# Utils
#----
.is_package_installed <- function(n) {
    names(n) <- n
    sapply(n, function(x) length(unlist(lapply(.libPaths(), function(lib) find.package(x, lib, quiet=TRUE, verbose=FALSE)))) > 0)
}
#----
.getProj <- function(x) {
    if (inherits(x,'Raster')) {
        if (!is.na(projection(x))) strsplit(strsplit(projection(x),'\\+proj=')[[1]][2],' ')[[1]][1]
        else {
            if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
            else 'projected'
        }
    } else {
        if (!is.na(crs(x))) strsplit(strsplit(crs(x,proj=TRUE),'\\+proj=')[[1]][2],' ')[[1]][1]
        else {
            if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
            else 'projected'
        }
    }
    
}
#-----
.mnwm <- function(d1, d2, d3, d4, d5, d6){
    X <- sum(c(d1, d2*2, d3, d4, d5*2, d6), na.rm = T)
    w <- sum(c(1,2,1,1,2,1) * is.finite(c(d1, d2, d3, d4, d5, d6)))
    return(X/w)
}
#-----
.ang <- function(dx, dy){
    ifelse(dy < 0, 180 + .deg(atan(dx/dy)),
           ifelse(dx < 0, 360 + .deg(atan(dx /dy )), .deg(atan(dx/dy))))
}
#---
.is.projected <- function(x) {
    if (inherits(x,'SpatRaster')) {
        e <- as.vector(terra::ext(x))
    } else e <- as.vector(extent(x))
    
    !all(e >= -180 & e <= 180)
}
#---
.rad <- function (degree) {
    (degree * pi) / 180
}
#---
.deg <-  function (radian) {
    (radian * 180) / pi
}
