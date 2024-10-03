
range_shift <- function(x,
                        periods = NULL,
                        SA = NULL,
                        proj_equal = TRUE,
                        probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                        raster_size_tolerance=NULL){
    
    if(is.null(periods)){
        stop("Periods associated with time series must be provided!")
    }
    
    crs_x <- crs(x)
    
    if(class(SA)=="SpatVector"){
        random_points <- terra::spatSample(SA, raster_size_tolerance)
    } 
    
    if(is.null(raster_size_tolerance)){
        x <- terra::as.data.frame(x, na.rm = TRUE, xy = TRUE)
    } else {
        if(terra::ncell(x) > raster_size_tolerance){
            if(class(SA)=="SpatVector"){
                x <- terra::extract(x, random_points, xy = TRUE, ID = FALSE)
            } else {
                x <- terra::spatSample(x, raster_size_tolerance, xy = TRUE, na.rm = TRUE)
            }
        } else {
            x <- terra::as.data.frame(x, na.rm = TRUE, xy = TRUE)
        }
    }
    x <- terra::vect(x, geom=c("x","y"), crs = crs_x)
    
    if(proj_equal){
        ## Project to equal area for more accurate statistics
        x <- terra::project(x, Eckt)
    }
    x <- terra::as.data.frame(x, geom = "XY")
    
    duration <- as.numeric(periods)
    duration <- max(duration) - min(duration)
    ####################################
    # calculate centroid shift
    # get centroid coordinates for each year
    centroids <- lapply(periods, function(i){
        sui_xy <- na.omit(data.frame(x=x$x,y=x$y,z=x[,i]))
        data.frame(year=i, 
                   centroid_x = weighted.mean(sui_xy[,1],sui_xy[,3]), 
                   centroid_y = weighted.mean(sui_xy[,2],sui_xy[,3]))
    })
    centroids <- do.call(rbind,centroids)
    centroids$year = as.numeric(centroids$year)
    centroids <- na.omit(centroids)
    
    centroids_xy <- terra::vect(centroids, geom=c("centroid_x","centroid_y"), crs = crs_x)
    # plot(x[[1]]);plot(centroids_xy,add=TRUE);dev.off()
    
    # centroid shift - the magnitude of the resulting vector from lat long shifts
    centroid_shift_lat_lm_m <- lm(centroid_y~as.numeric(year), centroids)
    centroid_shift_lat_lm_m <- summary(centroid_shift_lat_lm_m)
    centroid_shift_lat_lm <- centroid_shift_lat_lm_m$coef[2,1]
    
    centroid_shift_long_lm_m <- lm(centroid_x~as.numeric(year), centroids)
    centroid_shift_long_lm_m <- summary(centroid_shift_long_lm_m)
    centroid_shift_long_lm <- centroid_shift_long_lm_m$coef[2,1]
    
    centroid_shift_xy_lm <- sqrt(centroid_shift_lat_lm^2 + centroid_shift_long_lm^2)
    centroid_shift_xy_angle <- .ang(centroid_shift_long_lm, centroid_shift_lat_lm)
    
    # centroid shift - The total geographical distance between centroids in t1 and tn (tn is the last year)
    centroid_shift_xy_total <- as.numeric(terra::distance(centroids_xy[c(1,length(centroids_xy))])) / duration
    # centroid latitudinal shift  - The total latitudinal distance between centroids in t1 and tn (tn is the last year)
    centroid_shift_lat_total <- abs(centroids$y[nrow(centroids_xy)] - centroids$y[1]) / duration
    
    
    # additional parameters of the centroid latitudinal shift
    centroid_shift_lat_lm_error <- centroid_shift_lat_lm_m$coef[2,2]
    centroid_shift_lat_lm_r2 <- centroid_shift_lat_lm_m$r.squared
    centroid_shift_lat_lm_pvalue <- centroid_shift_lat_lm_m$coef[2,4]
    
    centroid_shift <- data.frame(centroid_shift_xy_lm,
                                 centroid_shift_xy_angle,
                                 centroid_shift_lat_lm,
                                 centroid_shift_lat_lm_error,
                                 centroid_shift_lat_lm_r2,
                                 centroid_shift_lat_lm_pvalue,
                                 centroid_shift_xy_total,
                                 centroid_shift_lat_total)
    
    periods <- as.character(centroids$year)
    ####################################
    # calculate range edge shift
    # get edge coordinates for each year
    edges <- lapply(periods, function(i){
        tmp_y <- wtd.quantile(x$y, x[,i], probs = probs, na.rm = TRUE) %>% t %>% data.frame()
        tmp_x <- wtd.quantile(x$x, x[,i], probs = probs, na.rm = TRUE) %>% t %>% data.frame()
        
        names(tmp_x) <- paste0("x_",as.character(probs))
        names(tmp_y) <- paste0("y_",as.character(probs))
        
        data.frame(tmp_x, tmp_y)
        
    })
    edges <- do.call(rbind,edges)
    
    # divide edges by prob
    edges <- lapply(probs, function(i){
        tmp <- edges[,grep(i,names(edges))]
        names(tmp) <- c("edge_x","edge_y")
        data.frame(year = as.numeric(periods), tmp)
    })
    names(edges) <- as.character(probs)
    
    
    # calculate range edge shift for each prob
    edge_shift <- lapply(probs, function(i){
        
        i <- as.character(i)
        edges_i <- edges[[i]]
        edges_xy <- terra::vect(edges_i, geom=c("edge_x","edge_y"), crs = crs_x)
        
        # edge shift - the magnitude of the resulting vector from lat long shifts
        edge_shift_lat_lm_m <- lm(edge_y~as.numeric(year), edges_i)
        edge_shift_lat_lm_m <- summary(edge_shift_lat_lm_m)
        edge_shift_lat_lm <- edge_shift_lat_lm_m$coef[2,1]
        
        edge_shift_long_lm_m <- lm(edge_x~as.numeric(year), edges_i)
        edge_shift_long_lm_m <- summary(edge_shift_long_lm_m)
        edge_shift_long_lm <- edge_shift_long_lm_m$coef[2,1]
        
        edge_shift_xy_lm <- sqrt(edge_shift_lat_lm^2 + edge_shift_long_lm^2)
        edge_shift_xy_angle <- .ang(edge_shift_long_lm, edge_shift_lat_lm)
        
        # edge shift - The total geographical distance between edges in t1 and tn (tn is the last year)
        edge_shift_xy_total <- as.numeric(terra::distance(edges_xy[c(1,length(edges_xy))])) / duration
        # edge latitudinal shift  - The total latitudinal distance between edges in t1 and tn (tn is the last year)
        edge_shift_lat_total <- abs(edges_i$y[nrow(edges_xy)] - edges_i$y[1]) / duration
        
        
        # additional parameters of the edge latitudinal shift
        edge_shift_lat_lm_error <- edge_shift_lat_lm_m$coef[2,2]
        edge_shift_lat_lm_r2 <- edge_shift_lat_lm_m$r.squared
        edge_shift_lat_lm_pvalue <- edge_shift_lat_lm_m$coef[2,4]
        
        tmp <- data.frame(edge_shift_xy_lm,
                          edge_shift_xy_angle,
                          edge_shift_lat_lm,
                          edge_shift_lat_lm_error,
                          edge_shift_lat_lm_r2,
                          edge_shift_lat_lm_pvalue,
                          edge_shift_xy_total,
                          edge_shift_lat_total)
        names(tmp) <- paste(names(tmp),i,sep = "_")
        return(tmp)
    })
    
    # merge all results
    edge_shift <- do.call(cbind,edge_shift)
    
    all_shifts <- cbind(centroid_shift,edge_shift)
    
    all_edges <- lapply(1:length(edges), function(i){
        tmp <- edges[[i]][,2:3]
        names(tmp) <- paste(names(tmp),probs[i],sep = "_")
        return(tmp)
    })
    all_edges <- do.call(cbind,all_edges)
    all_edges <- data.frame(centroids,all_edges)
    
    return(list(all_shifts,
                all_edges))
}
