# load ecoregions
get_ecoregions <- function(realm, PresAbs){
    if(realm=="Ter"){
        
        varsdir = "/media/seagate/boliveira/Land"
        
        ter.ras <- terra::rast(here::here(varsdir,"model_raster_ter.tif"))
        # check if shape was downloaded and if not download it
        if(dir.exists(here::here("Data/WWF_ecoregions"))){
            BA_shp <- terra::vect(here::here("Data/WWF_ecoregions/official/wwf_terr_ecos.shp"))
        } else {
            BA_shp <- WWFload(here::here("Data"))
        }
        # check if raster file exists, if not create and save
        if(file.exists(here::here("Data/WWF_ecoregions/WWF_ecoregions.tif"))){
            BA <- terra::rast(here::here("Data/WWF_ecoregions/WWF_ecoregions.tif"))
        } else {
            BA <- terra::rasterize(BA_shp, ter.ras, field = "ECO_NUM")
            writeRaster(BA, here::here("Data/WWF_ecoregions/WWF_ecoregions.tif"), overwrite = TRUE)
        }
    }
    
    if(realm=="Mar"){
        
        varsdir = "/media/seagate/boliveira/Marine"
        
        mar.ras <- terra::rast(here::here(varsdir,"model_raster_mar.tif"))
        # check if shape was downloaded and if not download it
        if(file.exists(here::here("Data/MEOW","MEOW.zip"))){
            BA_shp <- terra::vect(here::here("Data/MEOW", "DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp"))
        } else {
            download.file(url = "https://datadownload-production.s3.amazonaws.com/WCMC036_MEOW_PPOW_2007_2012_v1.zip",
                          destfile = here::here("Data/MEOW","MEOW.zip"))
            ## Unzip
            decompress_file(directory = "Data/MEOW", 
                            file = here::here("Data/MEOW","MEOW.zip"))
            BA_shp <- terra::vect(here::here("Data/MEOW", "DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp"))
        }
        # check if raster file exists, if not create and save
        if(file.exists(here::here("Data/MEOW/MEOW.tif"))){
            BA <- terra::rast(here::here("Data/MEOW/MEOW.tif"))
        } else {
            BA <- terra::rasterize(BA_shp, mar.ras, field = "PROVINC")
            writeRaster(BA, here::here("Data/MEOW/MEOW.tif"), overwrite = TRUE)
        }
    }
    
    tmpcells <- terra::cellFromXY(BA,xy=data.frame(PresAbs[,c('x','y')]))
    ecoreg <- BA[tmpcells]
    if(realm=="Ter"){ecoreg <- ecoreg$ECO_ID}
    if(realm=="Mar"){ecoreg <- ecoreg$PROVINC}
    ecoreg <- na.omit(unique(ecoreg))
    if(realm=="Ter"){BA <- BA %in% ecoreg}
    if(realm=="Mar"){BA <- BA %in% as.character(ecoreg)}
    BA <- ifel(BA,1,0)
    names(BA) <- "BA"
    
    return(BA)    
    
}
