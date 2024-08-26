# transfer files from meso to matrics

# folders to transfer
meso <- c("/lustre/oliveirab/SDMs/Mar",
          "/lustre/oliveirab/Data/Marine/oras/bio_proj",
          "/lustre/oliveirab/Data/Marine/oras/bio_proj_SA",
          "/lustre/oliveirab/Data/Marine/oras/model_raster_mar.tif",
          "/lustre/oliveirab/Data/Land/cruts/bio_proj_1km",
          "/lustre/oliveirab/Data/Land/cruts/model_raster_ter_1km.tif",
          "/lustre/oliveirab/SDMs/Ter")

matrics <- c("/scratch/boliveira/SDMs",
             "/scratch/boliveira/Data/Marine/oras",
             "/scratch/boliveira/Data/Marine/oras",
             "/scratch/boliveira/Data/Marine/oras",
             "/scratch/boliveira/Data/Land/cruts",
             "/scratch/boliveira/Data/Land/cruts",
             "/scratch/boliveira/SDMs")

downloads <- data.frame(meso,matrics)

for(i in 1:nrow(downloads)){
    system(paste(paste0("scp -rp oliveirab@muse-login02.meso.umontpellier.fr:", downloads$meso[i]), downloads$matrics[i]))
}

