# Script: Script to do the decimation process ALS (decrease of the density of points to make
# it comparable between dates).
# 
# Created by Gregoire Vincent and Claudia Huertas



# Libraries
library(raster)
library(lidR)


### Get DTM
# We start from a MNT of 50 cm.
MNT_50cm=raster("X:/lidar/ALS/Paracou/2019/MNT_0_5m/MNT50cm_Par2019.tif")
# Resample to 1 m
MNT_1m=aggregate(MNT_50cm, fac=2)
MNT_1m=crop(MNT_1m, c(284800, 287700,581200, 584300))
# MNT_50cm=crop(MNT_50cm, c(284700, 287800,581200, 584300))
# Specify the working directory
setwd("x:/lidar/ALS/Paracou/2019/Las/")

# Read the files LIDAR formatted in .las
file_list=dir(, pattern="las", full.names = F)

# Save the MNT model
writeRaster(MNT_1m, "d:/temp/ALSPara2019/CHM1m", format="GTiff", overwrite=T)

# Loop that gets the las file and returns it to 10 points per square meter.
csm_list=list()
for (f in 1:length(file_list))
{
  las = readLAS(file_list[f], select = "xyz", filter = "-keep_first -drop_abs_scan_angle_above 20")
  las<-lasclip(las, extent(MNT_1m))
  HN=lidR::lasnormalize(las, MNT_1m)
  HN@data[Z<0]$Z <-0
  HN@data<-HN@data[Z<60] #get rid of high points
  thinned10ppm = lasfilterdecimate(HN, random(10))
  csm_list[[f]]<-grid_canopy(thinned10ppm, res = 1, algorithm = p2r())
}

csm_full <- do.call(merge, csm_list)

# Output raster
writeRaster(csm_full, "d:/temp/ALSPara2019/CHM1m_10ppm_fo", format="GTiff", overwrite=T)
