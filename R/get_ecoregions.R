# load ecoregions
get_ecoregions <- function(realm, PresAbs, mask.ras, return.shp = FALSE, varsdir){
    
    if(realm=="Ter"){
        
        # check if shape was downloaded and if not download it
        if(dir.exists(here::here(varsdir,"WWF_ecoregions"))){
            BA_shp <- terra::vect(here::here(varsdir,"WWF_ecoregions/official/wwf_ecorealms.shp"))
        } else {
            BA_shp <- speciesgeocodeR::WWFload(here::here(varsdir))
            BA_shp <- terra::vect(BA_shp)
            # Create the EcoRealm layer
            # but first fix an issue >> there is conflict with the code NA for Neartic. Fix to avoid confusion with real.NA.
            tmp <- as.character(BA_shp$REALM)
            tmp[which(BA_shp$REALM=="NA")] <- "NEA" 
            BA_shp$EcoRealm <- paste(BA_shp$ECO_ID, tmp)
            writeVector(BA_shp, here::here(varsdir,"WWF_ecoregions/official/wwf_ecorealms.shp"), overwrite = TRUE)
            
            # plot(BA_shp, "EcoRealm", col=rainbow(12))
            # dev.off()
            
        }
        # check if raster file exists, if not create and save
        if(file.exists(here::here(varsdir,"WWF_ecoregions/official/EcoRealm.tif"))){
            ECO <- terra::rast(here::here(varsdir,"WWF_ecoregions/official/EcoRealm.tif"))
        } else {
            ECO <- terra::rasterize(BA_shp, mask.ras, field = "EcoRealm")
            # remove Antarctica from mask, which is then be used to remove from BA 
            ext.tmp <- ext(mask.ras)
            ext.tmp[4] <- -55
            rem <- cells(mask.ras,ext.tmp)
            mask.ras[rem] <- 0
            # This removes cells without environmental data & Antarctica
            ECO[mask.ras==0] <- NA
            # save
            writeRaster(ECO, here::here(varsdir,"WWF_ecoregions/official/EcoRealm.tif"), overwrite = TRUE)
        }
    }
    
    if(realm=="Mar"){
        
        # create dir to save ecoregion 
        if(!dir.exists(here::here(varsdir,"MEOW"))){
            dir.create(here::here(varsdir,"MEOW"),recursive = TRUE)
        }
        
        # check if shape was downloaded and if not download it
        if(file.exists(here::here(varsdir,"MEOW","MEOW.zip"))){
            BA_shp <- terra::vect(here::here(varsdir,"MEOW", "DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp"))
        } else {
            download.file(url = "https://datadownload-production.s3.amazonaws.com/WCMC036_MEOW_PPOW_2007_2012_v1.zip",
                          destfile = here::here(varsdir,"MEOW","MEOW.zip"))
            ## Unzip
            decompress_file(directory = here::here(varsdir,"MEOW"), 
                            file = here::here(varsdir,"MEOW","MEOW.zip"))
            BA_shp <- terra::vect(here::here(varsdir,"MEOW", "DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp"))
        }
        # check if raster file exists, if not create and save
        if(file.exists(here::here(varsdir,"MEOW/MEOW.tif"))){
            ECO <- terra::rast(here::here(varsdir,"MEOW/MEOW.tif"))
        } else {
            ECO <- terra::rasterize(BA_shp, mask.ras, field = "PROVINC")
            # Removes cells without environmental data
            ECO[mask.ras==0] <- NA
            # save
            writeRaster(ECO, here::here(varsdir,"MEOW/MEOW.tif"), overwrite = TRUE)
        }
    }
    
    # Get a pres/abs BA
    tmpcells <- terra::cellFromXY(ECO, xy=data.frame(PresAbs[,c('x','y')]))
    # which ecoregions/realms the species occur?
    if(realm=="Ter"){
        # raster
        ecoreg <- ECO[tmpcells]
        ecoreg <- as.character(ecoreg$EcoRealm)
        ecoreg <- na.omit(unique(ecoreg))
        # Crop rasters to the ecoregions/realms the species occur
        ECO <- ECO %in% ecoreg
        
        # shapes
        if(return.shp){
            ecoreg_shp <- terra::as.polygons(ECO) %>%
                filter(EcoRealm == 1)
        }
    }
    if(realm=="Mar"){
        # raster
        ecoreg <- ECO[tmpcells]
        ecoreg <- as.character(ecoreg$PROVINC)
        ecoreg <- na.omit(unique(ecoreg))
        # Crop rasters to the ecoregions/realms the species occur
        ECO <- ECO %in% ecoreg
        
        # shapes
        if(return.shp){
            ecoreg_shp <- terra::as.polygons(ECO) %>%
                filter(PROVINC == 1)
        }
    }
    # from TRUE/FALSE to 1/0 
    ECO <- ifel(ECO,1,0)
    names(ECO) <- "BA"
    
    if(return.shp){
        return(list(ECO, ecoreg_shp))
    }
    else{
        return(ECO)
    }
}

