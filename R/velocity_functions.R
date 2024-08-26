# Authors: Brunno Oliveira
# Last update :  Oct 2023
#--------
# This is based on the functions from the climetrics package (https://github.com/shirintaheri/climetrics)


#--------------
# temporal gradients using time series:
# threshold of minimun N obs are considered to get slope 
temp_gradFun <- function(x, th) {
  x <- x[!is.na(x)]
  if (length(x) > th) {
    s <- lm(x~c(1:length(x)))
    s$coefficients[2]
  } else NA
}
temp_grad <- function(x, th, ncores=NULL, tempfile="",overwrite=TRUE) {
  if(is.null(ncores)){
    tmp <- terra::app(x,temp_gradFun,th=th,filename=tempfile,overwrite=overwrite)
  } else {
    tmp <- terra::app(x,temp_gradFun,th=th,cores=ncores,filename=tempfile,overwrite=overwrite)
  }
  names(tmp) <- "Trend"
  return(tmp)
}

#--------------
# spatial gradient
spatial_grad <- function(rx, y_diff = 1, unit_out = "km") {
  
  if(!(unit_out == "km" | unit_out == "m")){
    stop("unit_out should be 'km' or 'm'")
  }
  
  if(nlyr(rx) > 1){ rx <- mean(rx,na.rm = TRUE) }
  
  if (.getProj(rx) == 'longlat') {
    if(unit_out=="km"){
      y_dist <- d2km(res(rx)) # from degrees to km
    }
    if(unit_out=="m"){
      y_dist <- d2km(res(rx))*1000 # from degrees to m
    }
    
  } else {
    if(unit_out=="km"){
      y_dist <- res(rx) / 1000 # from meters to km
    }
    if(unit_out=="m"){
      y_dist <- res(rx)
    }
    y_diff <- NA
  }
  
  if (!.is_package_installed("dplyr") || !.is_package_installed('tidyr')) stop('The packages dplyr and tidyr are needed for this metric; Please make sure they are installed!')
  
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
  
  y3b <- eval(parse(text="dplyr::select(y,from, code1, sst)"),envir =environment())
  y3b <- eval(parse(text="tidyr::spread(y3b,code1, sst)"),envir =environment())
  y3b$sstFocal <- rx[y3b$from][,1]
  y3b$LAT <- yFromCell(rx, y3b$from)
  
  if(!is.na(y_diff)) {
    y3b <- eval(parse(text="dplyr::mutate(y3b,
                         latpos = cos(deg_to_rad(LAT + y_diff)),
                         latneg = cos(deg_to_rad(LAT - y_diff)),
                         latfocal = cos(deg_to_rad(LAT)))"),envir =environment())
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
  
  rx <- c(rx, rx, rx, rx)
  names(rx) <- names(y3c)[-1]
  rx[[1]][y3c$icell] <- y3c$WE
  rx[[2]][y3c$icell] <- y3c$NS
  rx[[3]][y3c$icell] <- y3c$angle
  rx[[4]][y3c$icell] <- y3c$Grad
  
  return(rx)
}

#----------
# truncate is for bounding max and min values to upper (95%) and lower (5%) quantiles, respectively
gVelocity <- function(grad, slope, grad_col = "Grad", truncate=FALSE) {
  
  v <- slope
  g <- grad[grad_col]
  v_ang <- grad["angle"]
  
  # velocity angles have opposite direction to the spatial climatic gradient if warming and same direction (cold to warm) if cooling
  ind <- cells(v > 0)
  v_ang[ind] <- (v_ang[ind] + 180) %% 360
  
  # calculate velocity
  v <- v / g
  
  if(truncate){
    .o <- as.matrix(global(v,fun=quantile,probs=c(0.05,0.95),na.rm=TRUE))[1,]
    v[v < .o[1]] <- .o[1]
    v[v > .o[2]] <- .o[2] 
  }
  
  output <- c(v,v_ang)
  names(output) <- c("Vel", "Ang")
  return(output)
}

#----
# Utils
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
  ifelse(dy < 0, 180 + rad_to_deg(atan(dx/dy)),
         ifelse(dx < 0, 360 + rad_to_deg(atan(dx /dy )), rad_to_deg(atan(dx/dy))))
}
#---
deg_to_rad <- function (degree) {
  (degree * pi) / 180
}
#---
rad_to_deg <-  function (radian) {
  (radian * 180) / pi
}
#---
d2km <- function (d, base.latitude = 1) 
{
  if (!requireNamespace("fields")) 
    stop("Required fields package is missing.")
  onerad_to_degree.dist <- fields::rdist.earth(matrix(c(0, base.latitude), ncol = 2), 
                                               matrix(c(1, base.latitude), ncol = 2), 
                                               miles = FALSE)[,1]
  out <- d * onerad_to_degree.dist
  return(out)
}
#---
angle_map <- function(x, main = ""){
  
  myramp1 <- colorRampPalette(colors = c("red","green","blue","yellow","red"))(360)
  
  the_palette_fc <- leaflet::colorNumeric(
    palette = myramp1, 
    domain = c(0,360),
    reverse = FALSE)
  
  myramp <- the_palette_fc(seq(0, 360, length.out = 50))
  
  layout(matrix(c(1,1,1,2), nrow = 1, ncol = 4, byrow=T))
  plot(x, col = myramp, main = main)
  plotrix::polar.plot(
    start = 90,
    lengths = c(rnorm(360,mean = 1, sd = 0.001)),
    polar.pos = seq(0,360,by=1),
    clockwise = TRUE,
    cex=.1,
    show.grid.labels=0,
    line.col=myramp1)
}
#---
velocity_map <- function(x, main = "", ..., ggplot=FALSE){
  
  x_range <- terra::minmax(x)
  the_palette_fc <- leaflet::colorNumeric(
    palette = "RdBu", 
    domain = c(-max(abs(x_range)),max(abs(x_range))),
    reverse = TRUE)
  
  the_colors <- the_palette_fc(seq(min(x_range), max(x_range), length.out = 50))
  
  breaks <- c(seq(x_range[1],0,length.out=4),seq(0,x_range[2],length.out=3))
  breaks <- breaks[-which(duplicated(breaks))]
  
  if(ggplot){
    ggplot() +
      geom_spatraster(data = x) + 
      ggtitle(label = main) +
      scale_fill_gradient2(name = "",
                           low = the_colors[1],
                           mid = "white",
                           high = the_colors[length(the_colors)], 
                           breaks = round(breaks,0),
                           limits = round(c(x_range[1], x_range[2]), 0),
                           na.value = "transparent") +
      theme_light() +
      theme(
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.box.background = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))
  } else {
    plot(x, main = main, col = the_colors, ...)
  }
}
#---
circ_mean <- function(x) { 
  x=na.omit(x)
  x = x * pi/180 
  x = atan2(mean(sin(x)), mean(cos(x))) * 180/pi
  return(ifelse(x < 0, x+360,x))
}