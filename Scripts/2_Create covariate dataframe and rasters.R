##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
## 1: Create raster covariates for JAGS multispp occ model ##
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++##

## Billy Geary
## March 2021
## Upper Warren Mammal Project

library(sf)
library(fasterize)
library(raster)
library(dplyr)
#### Step 1: Read in Trap Points ####
trap.points = read_sf("~/Dropbox/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Eradicat/Camera_Locations.shp")
trap.points = trap.points %>% filter(Treatment == "Ground") 

# Get dates of inital trapping and clean up covariate table
ct.table.site = read.csv("~/Dropbox/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Eradicat data/camera.operational.csv")
ct.cam.list = read.csv("~/Dropbox/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Eradicat data/Site.List.csv")
trap.points = left_join(trap.points, ct.cam.list, by=c("LocationNa" = "LocationName"))
trap.points = left_join(trap.points, ct.table.site, by=c("Study.AreaName" = "Site"))

# Constrain data to just Upper Warren
site.subset = c("Balban","Boyicup", "Chariup", "Dudijup", "Dwalgan","Meribup", "Murtin", "Tone", "Warrup", "Yerramin", "Yeticup")
trap.points = trap.points %>% filter(Study.AreaName %in% site.subset)

trap.points.jsdm = trap.points %>% 
  dplyr::select(UTM_E, UTM_N,LocationID.x, LocationNa, Study.AreaName, District, Round, Start, End) %>%
  mutate(Start = as.POSIXct(Start, format="%d/%m/%y"),
         End = as.POSIXct(End, format="%d/%m/%y")) %>%
  rename(X = "UTM_E", Y="UTM_N", LocationID = "LocationID.x", LocationName = "LocationNa")


#### Step 2: Create Study Area Mask ####
study.region = trap.points %>% st_buffer(dist=10000) %>% group_by() %>% summarise()
mask = raster(extent(study.region), res=100, crs = st_crs(study.region)) # Resolution of 100m x 100m 
study.mask = fasterize(study.region, mask)

writeRaster(study.mask, "Data_Clean/covariate_StudyRegion_StudyMask.tif", overwrite=TRUE)

#### Step 3: DEM covariate ####
# DEM from Geoscience Australia
dem = raster("~/Dropbox/_research/data/spatial/dem/dem-9s.asc")
crs(dem) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
dem = projectRaster(dem, crs=crs(trap.points))
dem.crop = crop(dem, study.mask)
dem.crop = resample(dem.crop, study.mask, method='ngb')

trap.points.jsdm$elev = extract(dem.crop, trap.points.jsdm)

writeRaster(dem.crop, "Data_Clean/covariate_StudyRegion_DEM.tif", overwrite=TRUE)

#### Step 4: Fire Covariates ####
# Fire History from data.wa.gov.au / DBCA
# Process will be to create tsf and fire count rasters and extract by date
# But then save raster at start of survey and end of survey
firehistory = read_sf("~/Dropbox/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Covariates/FireHistoryDec18_WarrenJarrah.shp")
firehistory = st_transform(firehistory, st_crs(trap.points.jsdm))
trap.dates =  na.omit(unique(trap.points.jsdm$Start))

# Create TSF raster for start of period
firehistory$TSF = as.numeric(difftime(trap.dates[1],firehistory$FIH_DATE1, unit="weeks"))/52.25
firehistory_temp = subset(firehistory, TSF > 0)
tsf_raster = fasterize(firehistory_temp, study.mask, field="TSF", fun='min', background=NA)
writeRaster(tsf_raster, "Data_Clean/covariate_StudyRegion_tsf.start.Oct16.tif", overwrite=TRUE)

# Create TSF raster for end of period
firehistory$TSF = as.numeric(difftime(trap.dates[20],firehistory$FIH_DATE1, unit="weeks"))/52.25
firehistory_temp = subset(firehistory, TSF > 0)
tsf_raster = fasterize(firehistory_temp, study.mask, field="TSF", fun='min', background=NA)
writeRaster(tsf_raster, "Data_Clean/covariate_StudyRegion_tsf.end.Oct17.tif")

# Extract survey specific tsf for each site
fire.list = list()
for (i in 1:length(trap.dates)) {
  start.date = trap.dates[i]
  firehistory$TSF = as.numeric(difftime(start.date,firehistory$FIH_DATE1, unit="weeks"))/52.25
  
  #time since fire
  firehistory_temp = subset(firehistory, TSF > 0)
  tsf_raster = fasterize(firehistory_temp, mask, field="TSF", fun='min')
  
  extract_traps = trap.points.jsdm %>% filter(Start == start.date)
  
  tsf.point = extract(tsf_raster, extract_traps)
  fire.data = cbind(extract_traps, tsf.point)
  
  fire.list[[as.character(start.date)]] = fire.data
}


trap.points.jsdm = do.call(rbind, fire.list)

#### Step 5: NV distance to edge ####
# Native Veg from data.wa.gov.au
native.veg = read_sf("~/Dropbox/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Covariates/nativeveg_warrenjarrah.shp")
native.veg1 = native.veg %>% group_by() %>% summarise()
native.veg1 = native.veg1 %>% st_transform(st_crs(trap.points))
# Create raster of non-native vegetation by creating inverse of nv raster
nv.ras = fasterize(native.veg1, study.mask, fun = 'last')
nv.ras[nv.ras == 1] <- 2
nv.ras[is.na(nv.ras)] <- 1
nv.ras[nv.ras == 2] <- NA
# Create raster of distance to native vegetation
nv.dist = distance(nv.ras)
plot(nv.dist)
plot(trap.points, add=TRUE)

trap.points.jsdm$nv.dist = extract(nv.dist, trap.points.jsdm)

writeRaster(nv.dist, "Data_Clean/covariate_StudyRegion_NVdist.tif", overwrite=TRUE)

#### Step XX: Add to detection history df ####
dethist = readRDS("Data_Processing/camtrapR.dethist.UpperWarren.multispp.21112021.RData")
trap.points.jsdm = trap.points.jsdm %>% st_drop_geometry()
dethist$Sites = left_join(dethist$Sites, trap.points.jsdm)

saveRDS(dethist, "Data_Processing/camtrapR.dethist.UpperWarren.multispp.21112021.RData")

