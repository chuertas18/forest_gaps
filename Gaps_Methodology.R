rm(list = ls())

# Abstract: Different methodologies used for the detection of mortality from gap 
# dynamics analysis. (Benoist, 2009; Brokaw, 1985; Leitold et al., 2018; Senécal et al., 2018) 
# and Hunter et al. (2015)
# 
# Description: Brokaw (1985) Aperture extending at all levels up to an average height of 2 
# meters above the ground.
# Benoist (2009) Developed an algorithm to detect and delimit gaps. 
# "The threshold of 5m² minimum makes it possible to exclude low points that are 
# too isolated and that correspond to old gaps that are in the process of aggradation. 
# On the other hand, the 10m height criterion restores their shape well, and it is 
# reasonable to accept a high point in the middle of a hole." 
# Senécal et al. (2018). Height reduction areas (HRA) groups of pixels in which a 
# reduction of at minimum 1 m in height and a minimum size of 5 m2 were observed (six years). 
# Within the HRA, new gap pixels were detected as CHM year1 pixels in the HRA of height 3 m or 
# less, and for year0 heights above 3 m, while canopy height erosion pixels were detected as 
# CHM year1 pixels of height above 3 m after the height loss.
# Leitold et al., (2018), threshold >= 4 m² and height losses >= 3 m.
#
# Created by Claudia Huertas 

################################################################################################
################################################################################################
### FUNCTIONS
################################################################################################

# Variables
# chm0 - Canopy model year 1
# chm1 - Canopy model year 2
# outdir - Directory where the data will be stored
# outname - Output file name



######### new_gaps_brokaw ####################
#### Calculate the gaps, following Brokaw's (1982) methodology
# Brokaw (1982) -2 meters above the ground.
# Brokaw, N.V.L. (1985), Gap-Phase Regeneration in a Tropical Forest. Ecology, 
# 66: 682-687. https://doi.org/10.2307/1940529

new_gaps_brokaw<-function(chm0,chm1,outdir,outname){
  require(GISTools)
  require(rgdal)
  require(rgeos)
  chm1=raster(chm1)
  chm1_g=chm1
  chm1_g[chm1_g<=2]<-1
  chm1_g[chm1_g>2]<-0
  chm0=raster(chm0)
  chm0_g=chm0
  chm0_g[chm0_g<=2]<-0
  chm0_g[chm0_g>2]<-1
  delta<-chm1-chm0
  gaps_s<-chm1_g*chm0_g
  gaps_v<-delta*gaps_s
  #plot(gaps_s)
  # gaps_pol <- rasterToPolygons(gaps_s, fun=function(x){x==1},dissolve=TRUE) # Poligonizar y disolver
  # gaps_pol <- disaggregate(gaps_pol) # Disaggregate polygons
  # #length(gaps_pol)
  writeRaster(gaps_s, filename=paste0(outdir,"/",outname,"_surface"), format="GTiff", overwrite=TRUE)
  # writeRaster(gaps_v, filename=paste0(outdir,"/",outname,"_volumetric"), format="GTiff", overwrite=TRUE)
  # writeOGR(obj=gaps_pol, dsn=outdir, layer=outname, driver="ESRI Shapefile", overwrite=TRUE) # this is in geographical projection
  # 
}

## new_gaps_senecal = chm1 <=3 m chm0 >3 m
## -	New gap pixels were detected as CHM year1 pixels in the HRA of height 3 m or less, 
# and for year0 heights above 3 m
# Senécal, J.-F., Doyon, F., Messier, C., 2018. Tree Death Not Resulting in Gap 
# Creation: An Investigation of Canopy Dynamics of Northern Temperate Deciduous 
# Forests. Remote Sens. 10, 121. https://doi.org/10.3390/rs10010121

new_gaps_senecal<-function(chm0,chm1,outdir,outname){
  require(GISTools)
  require(rgdal)
  require(rgeos)
  chm1=raster(chm1)
  chm1_g=chm1
  chm1_g[chm1_g<=3]<-1
  chm1_g[chm1_g>3]<-0
  chm0=raster(chm0)
  chm0_g=chm0
  chm0_g[chm0_g<=3]<-0
  chm0_g[chm0_g>3]<-1
  gaps<-chm1_g*chm0_g
  #plot(gaps)
  delta<-chm1-chm0
  gaps_s<-chm1_g*chm0_g
  gaps_v<-delta*gaps_s
  #plot(gaps_s)
  # gaps_pol <- rasterToPolygons(gaps_s, fun=function(x){x==1},dissolve=TRUE) # Poligonizar y disolver
  # gaps_pol <- disaggregate(gaps_pol) # Disaggregate polygons
  # #length(gaps_pol)
  writeRaster(gaps_s, filename=paste0(outdir,"/",outname,"_surface"), format="GTiff", overwrite=TRUE)
  # writeRaster(gaps_v, filename=paste0(outdir,"/",outname,"_volumetric"), format="GTiff", overwrite=TRUE)
  # writeOGR(obj=gaps_pol, dsn=outdir, layer=outname, driver="ESRI Shapefile", overwrite=TRUE) # this is in geographical projection
}

##### Vincent Benoist - 5m² , 10m height criterion
## new_gaps_benoist =   delta_area >=  5m² , chm1 <=10 m
# Benoist, V., 2009. Étude de la mortalité des arbres en forêt tropicale humide par laser aéroporté.
new_gaps_benoist<-function(chm0,chm1,outdir,outname){
  require(GISTools)
  require(rgdal)
  require(rgeos)
  chm1=raster(chm1)
  chm1_g=chm1
  chm1_g[chm1_g<=10]<-1 # new_gaps_benoist =   chm1 <=10 m
  chm1_g[chm1_g>10]<-0
  chm0=raster(chm0)
  chm0_g=chm0
  r=chm1-chm0 # delta
  r[r>=0]<-0
  #plot(r)
  
  # extend r with a number of rows and culomns (at each side)
  # to isolate clumps adjacents to plot axes 
  r2<-extend(r, c(1,1))
  plot(r2)
  rc <- clump(r2, directions = 8) 
  #plot(rc)
  
  # get frequency table    
  f<-freq(rc)
  # save frequency table as data frame
  f<-as.data.frame(f)
  
  clumps_pixels=5 # Vincent Benoist - 5m² 
  
  # which rows of the data.frame are only represented by clumps under 9pixels?
  str(which(f$count <= clumps_pixels))
  # which values do these correspond to?
  str(f$value[which(f$count <= clumps_pixels)])
  # put these into a vector of clump ID's to be removed
  excludeID <- f$value[which(f$count <= clumps_pixels)]
  
  # make a new raster to be sieved
  formaskSieve <- rc
  # assign NA to all clumps whose IDs are found in excludeID
  formaskSieve[rc %in% excludeID] <- NA
  
  #plot(formaskSieve)
  formaskSieve[is.na(formaskSieve)] <- 0
  formaskSieve[formaskSieve>0]<-1
  
  raster_ip=r*formaskSieve
  raster_ip[is.na(raster_ip)] <- 0
  #plot(raster_ip)
  # return(raster_ip)
  
  # gaps_pol <- rasterToPolygons(formaskSieve, fun=function(x){x==1},dissolve=TRUE) # Poligonizar y disolver
  # gaps_pol <- disaggregate(gaps_pol) # Disaggregate polygons
  # 
  writeRaster(formaskSieve, filename=paste0(outdir,"/",outname,"_surface"), format="GTiff", overwrite=TRUE)
  # writeRaster(raster_ip, filename=paste0(outdir,"/",outname,"_volumetric"), format="GTiff", overwrite=TRUE)
  # writeOGR(obj=gaps_pol, dsn=outdir, layer=outname, driver="ESRI Shapefile", overwrite=TRUE) # this is in geographical projection
  # 
}


#### hra_senecal = delta_area >=  5 m² and  delta_h >= 1 m (height losses between chm)
# Senécal, J.-F., Doyon, F., Messier, C., 2018. Tree Death Not Resulting in Gap 
# Creation: An Investigation of Canopy Dynamics of Northern Temperate Deciduous 
# Forests. Remote Sens. 10, 121. https://doi.org/10.3390/rs10010121
hra_senecal<-function(chm0,chm1,outdir,outname){
  require(GISTools)
  require(rgdal)
  require(rgeos)
  chm1=raster(chm1)
  chm1_g=chm1
  chm0=raster(chm0)
  chm0_g=chm0
  r=chm1-chm0 # delta
  r[r>=-1]<-0 # delta_h >= 1 m (height losses between chm)
  #plot(r)
  
  # extend r with a number of rows and culomns (at each side)
  # to isolate clumps adjacents to plot axes 
  r2<-extend(r, c(1,1))
  #plot(r2)
  rc <- clump(r2, directions = 8) 
  #plot(rc)
  
  # get frequency table    
  f<-freq(rc)
  # save frequency table as data frame
  f<-as.data.frame(f)
  
  clumps_pixels=5 # Vincent Benoist - 5m² 
  
  # which rows of the data.frame are only represented by clumps under 9pixels?
  str(which(f$count <= clumps_pixels))
  # which values do these correspond to?
  str(f$value[which(f$count <= clumps_pixels)])
  # put these into a vector of clump ID's to be removed
  excludeID <- f$value[which(f$count <= clumps_pixels)]
  
  # make a new raster to be sieved
  formaskSieve <- rc
  # assign NA to all clumps whose IDs are found in excludeID
  formaskSieve[rc %in% excludeID] <- NA
  
  #plot(formaskSieve)
  formaskSieve[is.na(formaskSieve)] <- 0
  formaskSieve[formaskSieve>0]<-1
  
  raster_ip=r*formaskSieve
  raster_ip[is.na(raster_ip)] <- 0
  #plot(raster_ip)
  # return(raster_ip)
  
  # gaps_pol <- rasterToPolygons(formaskSieve, fun=function(x){x==1},dissolve=TRUE) # Poligonizar y disolver
  # gaps_pol <- disaggregate(gaps_pol) # Disaggregate polygons
  # 
  writeRaster(formaskSieve, filename=paste0(outdir,"/",outname,"_surface"), format="GTiff", overwrite=TRUE)
  # writeRaster(raster_ip, filename=paste0(outdir,"/",outname,"_volumetric"), format="GTiff", overwrite=TRUE)
  # writeOGR(obj=gaps_pol, dsn=outdir, layer=outname, driver="ESRI Shapefile", overwrite=TRUE) # this is in geographical projection
  # 
}

#### gaps_leitold 
# gaps_leitold = delta_area >=  4 m² delta_h >= 3 m (height losses between chm) . 
# Note for the definition of the different types of clearings 
# Leitold uses volume and makes calculation from vector
# Leitold, V., Morton, D.C., Longo, M., dos-Santos, M.N., Keller, M., Scaranello, 
# M., 2018. El Niño drought increased canopy turnover in Amazon forests. New Phytol. 
# 219, 959-971. https://doi.org/10.1111/nph.15110
gaps_leitold <-function(chm0,chm1,outdir,outname){
  require(GISTools)
  require(rgdal)
  require(rgeos)
  chm1=raster(chm1)
  chm1_g=chm1
  chm0=raster(chm0)
  chm0_g=chm0
  r=chm1-chm0 # delta
  r[r>=-3]<-0 # delta_h >= 3 m (height losses between chm)
  #plot(r)
  
  # extend r with a number of rows and culomns (at each side)
  # to isolate clumps adjacents to plot axes 
  r2<-extend(r, c(1,1))
  #plot(r2)
  rc <- clump(r2, directions = 8) 
  #plot(rc)
  
  # get frequency table    
  f<-freq(rc)
  # save frequency table as data frame
  f<-as.data.frame(f)
  
  clumps_pixels=4 # 4m² 
  
  # which rows of the data.frame are only represented by clumps under 9pixels?
  str(which(f$count <= clumps_pixels))
  # which values do these correspond to?
  str(f$value[which(f$count <= clumps_pixels)])
  # put these into a vector of clump ID's to be removed
  excludeID <- f$value[which(f$count <= clumps_pixels)]
  
  # make a new raster to be sieved
  formaskSieve <- rc
  # assign NA to all clumps whose IDs are found in excludeID
  formaskSieve[rc %in% excludeID] <- NA
  
  #plot(formaskSieve)
  formaskSieve[is.na(formaskSieve)] <- 0
  formaskSieve[formaskSieve>0]<-1
  
  raster_ip=r*formaskSieve
  raster_ip[is.na(raster_ip)] <- 0
  #plot(raster_ip)
  # return(raster_ip)
  
  # gaps_pol <- rasterToPolygons(formaskSieve, fun=function(x){x==1},dissolve=TRUE) # Poligonizar y disolver
  # gaps_pol <- disaggregate(gaps_pol) # Disaggregate polygons
  # 
  writeRaster(formaskSieve, filename=paste0(outdir,"/",outname,"_surface"), format="GTiff", overwrite=TRUE)
  # writeRaster(raster_ip, filename=paste0(outdir,"/",outname,"_volumetric"), format="GTiff", overwrite=TRUE)
  # writeOGR(obj=gaps_pol, dsn=outdir, layer=outname, driver="ESRI Shapefile", overwrite=TRUE) # this is in geographical projection
  
}

############################################################################################################################

# The function defined in Hunter's article is used.
# Hunter, M.O., Keller, M., Morton, D., Cook, B., Lefsky, M., Ducey, M., Saleska, 
# S., de Oliveira, R.C., Schietti, J., 2015. Structural Dynamics of Tropical Moist 
# Forest Gaps. PLOS ONE 10, e0132144. https://doi.org/10.1371/journal.pone.0132144


gaps_hunter<-function(chm,threshold_hunter,size_hunter,outdir,outname){
  require(GISTools)
  require(rgdal)
  require(rgeos)
  require(ForestGapR)
  chm=raster(chm)
  raster_hunter<-getForestGaps(chm_layer=chm, threshold=threshold_hunter, size=size_hunter)
  writeRaster(raster_hunter, filename=paste0(outdir,"/",outname), format="GTiff", overwrite=TRUE)
  formaskSieve<-raster_hunter
  formaskSieve[is.na(formaskSieve)] <- 0
  formaskSieve[formaskSieve>0]<-1
  gaps_pol <- rasterToPolygons(formaskSieve, fun=function(x){x==1},dissolve=TRUE) # Poligonizar y disolver
  gaps_pol <- disaggregate(gaps_pol)
  writeOGR(obj=gaps_pol, dsn=outdir, layer=outname, driver="ESRI Shapefile", overwrite=TRUE) # this is in geographical projection
  
}



clump_pixel<-function (chm0,chm1,dif_elev,clumps_pixels){
  library(raster)
  # chm0 =paste0(dir_chm,"/",nom_chm,"_FO_PAR_",year0,"_Plot_",i,".tif")
  chm0 =raster(chm0)
  #plot(chm0)
  range(chm0[],na.rm = T)
  
  # chm1 =paste0(dir_chm,"/",nom_chm,"_FO_PAR_",year1,"_Plot_",i,".tif")
  chm1 =raster(chm1)
  range(chm1[],na.rm = T)
  #plot(chm1)
  
  #reproducible example
  delta=chm1-chm0
  #plot(delta)
  r=delta
  r[r>=(dif_elev*-1)]<-NA
  #plot(r)
  
  # extend r with a number of rows and culomns (at each side)
  # to isolate clumps adjacents to plot axes 
  r2<-extend(r, c(1,1))
  #plot(r2)
  rc <- clump(r2, directions = 8) 
  #plot(rc)
  
  # get frequency table    
  f<-freq(rc)
  # save frequency table as data frame
  f<-as.data.frame(f)
  
  # which rows of the data.frame are only represented by clumps under 9pixels?
  str(which(f$count <= clumps_pixels))
  # which values do these correspond to?
  str(f$value[which(f$count <= clumps_pixels)])
  # put these into a vector of clump ID's to be removed
  excludeID <- f$value[which(f$count <= clumps_pixels)]
  
  # make a new raster to be sieved
  formaskSieve <- rc
  # assign NA to all clumps whose IDs are found in excludeID
  formaskSieve[rc %in% excludeID] <- NA
  
  #plot(formaskSieve)
  formaskSieve[is.na(formaskSieve)] <- 0
  formaskSieve[formaskSieve>0]<-1
  
  raster_ip=delta*formaskSieve
  raster_ip[is.na(raster_ip)] <- 0
  #plot(raster_ip)
  return(raster_ip)
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################
mnc0_max="Y:/users/ClaudiaHuertas/Mortality/Data/raster/CHM/decimation/CHM2015_10ppm_fo.tif"
mnc1_max="Y:/users/ClaudiaHuertas/Mortality/Data/raster/CHM/decimation/CHM2019_10ppm_fo.tif"
year0=2015
year1=2019
dec=TRUE

new_gaps_brokaw(mnc0_max,mnc1_max,"Y:/users/ClaudiaHuertas/Mortality/Data/GAPS/dec",paste0("HmaxCHM_FO_dec",dec,"_PAR_",year0,"_",year1,"_Brokaw"))
new_gaps_senecal(mnc0_max,mnc1_max,"Y:/users/ClaudiaHuertas/Mortality/Data/GAPS/dec",paste0("HmaxCHM_FO_dec",dec,"_PAR_",year0,"_",year1,"_Senecal"))
new_gaps_benoist(mnc0_max,mnc1_max,"Y:/users/ClaudiaHuertas/Mortality/Data/GAPS/dec",paste0("HmaxCHM_FO_dec",dec,"_PAR_",year0,"_",year1,"_Benoist"))
hra_senecal(mnc0_max,mnc1_max,"Y:/users/ClaudiaHuertas/Mortality/Data/GAPS/dec",paste0("HmaxCHM_FO_dec",dec,"_PAR_",year0,"_",year1,"_hra_Senecal"))
gaps_leitold(mnc0_max,mnc1_max,"Y:/users/ClaudiaHuertas/Mortality/Data/GAPS/dec",paste0("HmaxCHM_FO_dec",dec,"_PAR_",year0,"_",year1,"_Leitold"))


mnc0_max="Y:/users/ClaudiaHuertas/Mortality/Data/raster/CHM/decimation/CHM2009_10ppm_avsepoct.tif"
mnc1_max="Y:/users/ClaudiaHuertas/Mortality/Data/raster/CHM/decimation/CHM2015_10ppm_fo.tif"
year0=2009
year1=2015
dec=TRUE
