
# Get bioclimatic variables from species occurrences and backgrounds records

# sp_name = species name. Used for saving output with the species name
# occ = species occurrence records

# sp_name=sptogo
# occu=PresAbs
# layers=tiles
# min_year_layers = min_year_layers[[realm]]
# realm = realm
# dir_env_vars = varsdir[[realm]]
# myvars = myvars[[realm]]
# n_yr_bioclimatic = n_yr_bioclimatic
# limit_occ_recs = NULL
# output_dir = output_dir
# map_sps_occ = "/media/seagate/boliveira/SDMs/Maps_occs"
# check_if_file_exists = FALSE
# return_data = FALSE

Get_Bioclimatics <- function(sp_name,
                             occu,
                             layers,
                             min_year_layers,
                             realm,
                             dir_env_vars, 
                             myvars,
                             n_yr_bioclimatic = 2,
                             limit_occ_recs = NULL,
                             output_dir = "/media/seagate/boliveira/SDMs/Env_data", 
                             map_sps_occ = NULL,
                             check_if_file_exists = F,
                             save_data = TRUE,
                             return_data = TRUE){
    
    if(!save_data & !return_data){
        stop("You must choose either saving or returning the data or both.")
    }
    
    tic.clearlog()
    tic()
    
    require(sp)
    require(raster)
    require(terra)
    require(snow)
    require(parallel)
    require(doParallel)
    require(pbapply)
    require(data.table)
    
    if(is.null(output_dir)){
        stop("Output directory (output_dir) should be provided")
    }
    if(!dir.exists(output_dir)){
        dir.create(output_dir,recursive = T)
    }
    if(check_if_file_exists){
        file.test <- here::here(output_dir,paste(gsub(" ","_",sp_name),realm,"bio.RDS",sep="_"))
        if(!file.exists(file.test)){
            WTD <- "run"
        } else {
            WTD <- "load"
        }
    } else {
        WTD <- "run"
    } 
    if(WTD == "load"){
        if(return_data){
            cat('load')
            bios <- readRDS(file.test)
            return(bios)
        } else {
            cat('Data for this species has been collected. Data is not returned because return_data = FALSE.')
        }
    }
    if(WTD == "run"){
        cat('\nCollecting data for', sp_name)
        n_tiles <- length(layers)
        # keep only occurrences with date >= min_year_layers + n_yr_bioclimatic
        min_date <- min_year_layers + n_yr_bioclimatic
        if(any(which(occu$year < min_date))){
            occu <- occu[-which(occu$year < min_date),]
        }
        # Fix dates
        dtf <- unlist(mclapply(occu$month, function(x) sprintf("%02d", x), mc.cores = n_tiles))
        occu$month <-  dtf # fix month
        dtf <- unlist(mclapply(1:nrow(occu), function(i) {paste(occu$year[i], occu$month[i], sep = "_")}, mc.cores = n_tiles))
        occu$id <- dtf
        # group occurrences by tiles
        cat("\nGroup occurrences by tile")
        xys <- vect(occu, geom = c("x","y"))
        tile_ext <- lapply(layers, function(x) terra::ext(rast(x[1]))) # get extent of each tile
        xys <- lapply(tile_ext, function(x) crop(xys, x)) # crop occurrences by tile
        # Remove tiles without occurrences
        cat("\nRemove tiles without occurrences")
        rem <- which(sapply(xys, function(x) length(geom(x))==0))
        if(any(rem)){
            xys <- xys[-rem]
            layers <- layers[-rem]
            n_tiles <- length(layers)
        }
        # remove cells with NAs in env data
        cat("\nRemove cells with NAs in env data")
        xys <- lapply(1:n_tiles, function(i){
            xy2go <- xys[[i]]
            layers_i <- layers[[i]]
            layers_i_tmp <- sapply(myvars, function(x) grep(x,layers_i)[1])
            layers_i_tmp <- layers_i[layers_i_tmp]
            layers_i_tmp <- terra::rast(layers_i_tmp)
            xyk <- terra::geom(xy2go)[,c("x","y")]
            if(class(xyk)[1]=="numeric"){
                xyk <- data.frame(x = xyk[1], y = xyk[2])
            } else {
                xyk <- data.frame(xyk)
            }
            tmpcells <- terra::cellFromXY(layers_i_tmp[[1]],xy=xyk)
            if(length(tmpcells)>0){
                env_tile_tmp <- layers_i_tmp[tmpcells]
                rem <- which(apply(env_tile_tmp, 1, function(x) any(is.na(x)))) 
                if(any(rem)){
                    xy2go <- xy2go[-rem]
                }
            }
            return(xy2go)
        })
        # reduce number of occurrences 
        if(!is.null(limit_occ_recs)){
            if(nrow(occu) > limit_occ_recs){
                cat("Reduce the number of occurrences \n")
                occu <- lapply(xys, terra::as.data.frame, geom = "XY")
                occu <- rbindlist(occu)
                set.seed(666)
                occu <- occu[sample(1:nrow(occu), limit_occ_recs), ]
                # re-group occurrences by tiles
                cat("\nGroup occurrences by tile")
                xys <- vect(occu, geom = c("x","y"))
                tile_ext <- lapply(layers, function(x) terra::ext(rast(x[1]))) # get extent of each tile
                xys <- lapply(tile_ext, function(x) crop(xys, x)) # crop occurrences by tile
            }
        }
        # Save map occurrences sp i
        if(!is.na(map_sps_occ)){
            if(!dir.exists(map_sps_occ)){
                dir.create(map_sps_occ,recursive = T)
            }
            map_back <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sp")
            occs_tmp <- rbindlist(lapply(xys, as.data.frame, geom = "XY"))
            occs_tmp <- vect(occs_tmp, geom = c("x","y"))
            pdf(here::here(map_sps_occ,paste0(gsub(" ","_",sp_name),"_",realm,".pdf")),width = 20, height = 10)
            plot(map_back, main = sp_name)
            plot(subset(occs_tmp, occs_tmp$pa==1),add = TRUE,cex=.5, col = "red")
            plot(subset(occs_tmp, occs_tmp$pa==0),add = TRUE,cex=.5, col = "gray")
            dev.off()
        }
        # Select layers within the time-frame of interest at each tile
        cat("\nSelect layers within the time-frame of interest")
        layers <- lapply(1:n_tiles, function(i){
            # get all possible dates
            dates <- as.data.frame(xys[[i]])
            if(any(duplicated(dates$id))){
                dates2go <- dates[-which(duplicated(dates$id)),]
            } else {
                dates2go <- dates
            }
            # get temporal range for each possible date
            timerange <- GetTimeRange(dates2go[, c("year","month")], n_yr_bioclimatic)
            timeframe <- unique(unlist(timerange))
            timeframe <- paste(timeframe,collapse = "|")
            # Select layers within the time-frame of interest
            layers_i <- layers[[i]]
            pos <- grep(timeframe, names(layers_i))
            return(layers_i[pos])
        })
        # get env data from occurrences & calculate bioclimatics
        cat("\nCalculate bioclimatics")
        bios <- mclapply(1:n_tiles, function(i){
            xy2go <- xys[[i]]
            layers_i <- layers[[i]]
            layers_i_names <- names(layers_i)
            # get cells
            layers_i_tmp <- layers_i[1]
            layers_i_tmp <- terra::rast(layers_i_tmp)
            xyk <- terra::geom(xy2go)[,c("x","y")]
            if(class(xyk)[1]=="numeric"){
                xyk <- data.frame(x = xyk[1], y = xyk[2])
            } else {
                xyk <- data.frame(xyk)
            }
            tmpcells <- terra::cellFromXY(layers_i_tmp[[1]],xy=xyk)
            # extract env data
            layers_i <- terra::rast(layers_i)
            names(layers_i) <- layers_i_names
            env_tile <- layers_i[tmpcells]
            # Get possible dates
            possibledates <- data.frame(xy2go)
            if(any(duplicated(possibledates$id))){
                possibledates <- possibledates[-which(duplicated(possibledates$id)),]
            }
            # calculate bioclimatics for each possible date
            biosclim <- lapply(1:nrow(possibledates), function(k) {
                possibledates_range_cell <- unlist(GetTimeRange(possibledates[k,], n_yr_bioclimatic))
                cells2get <- which(xy2go$id %in% possibledates[k,'id'])
                lyr2get <- sapply(names(env_tile), function(x) paste(strsplit(x,"_")[[1]][c(2,3)],collapse = "_"))
                lyr2get <- which(lyr2get %in% possibledates_range_cell)
                envtimecells <- env_tile[cells2get,lyr2get]
                xyk <- terra::geom(xy2go)[cells2get,c("x","y")]
                if(class(xyk)[1]=="numeric"){
                    xyk <- data.frame(x = xyk[1], y = xyk[2])
                } else {
                    xyk <- data.frame(xyk)
                }
                # calculate bioclimatics
                if(realm == "Mar"){
                    mybios <- bioclimatics_ocean(envtimecells, 
                                                 myvars = myvars, 
                                                 n_yr_bioclimatic = n_yr_bioclimatic)
                } else {
                    mybios <- bioclimatics_land(envtimecells, 
                                                n_yr_bioclimatic = n_yr_bioclimatic)
                }
                data.frame(cell = xy2go$cell[cells2get], 
                           xyk, 
                           year = xy2go$year[cells2get], 
                           month = xy2go$month[cells2get], 
                           tile = i, 
                           pa = xy2go$pa[cells2get], 
                           mybios)
            })
            biosclim <- rbindlist(biosclim)
            return(biosclim)
            
        }, mc.cores = n_tiles)
        bios <- rbindlist(bios)
        
        if(save_data){
            saveRDS(bios, here::here(output_dir,paste(gsub(" ","_",sp_name),realm,"bio.RDS",sep="_")))
        }
        
        # time log
        tmp <- toc(quiet = TRUE)
        tmp <- tmp$toc-tmp$tic
        if(tmp<=60){
            cat("\n", tmp, 'seconds elapsed')
        } else {
            cat("\n", tmp/60, 'minutes elapsed')
        }
        
        if(return_data){
            return(bios)
        }
        
    } 
}


