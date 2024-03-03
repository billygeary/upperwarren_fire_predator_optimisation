##+++++++++++++++++++++++++++++++++++++++++++++++++++++##
## 2: Create raster covariates Multi Species Modelling ##
##+++++++++++++++++++++++++++++++++++++++++++++++++++++##

## Billy Geary
## Feb 2022
## Upper Warren Mammal Project

#### Setup ####
library(sf)
library(fasterize)
library(lubridate)
library(raster)
library(exactextractr)
library(tidyverse)

std.crs = st_crs(32750) # Standard CRS to put everything into: https://epsg.io/32750
source(here::here('Scripts', 'covariate_helper_functions.R'))

#### Step 1: Read in Trap Points & Create Study Area Mask 
######## Step 1: Read in Trap Points ####
trap.points = read_sf("~/Library/CloudStorage/Dropbox/Billy/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Eradicat/Camera_Locations.shp")
trap.points = trap.points %>% filter(Treatment == "Ground") 

# Get dates of inital trapping and clean up covariate table
ct.table.site = read.csv("~/Library/CloudStorage/Dropbox/Billy/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Eradicat data/camera.operational.csv")
ct.cam.list = read.csv("~/Library/CloudStorage/Dropbox/Billy/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Eradicat data/Site.List.csv")
trap.points = left_join(trap.points, ct.cam.list, by=c("LocationNa" = "LocationName"))
trap.points = left_join(trap.points, ct.table.site, by=c("Study.AreaName" = "Site"))

# Constrain data to just Upper Warren
site.subset = c("Balban","Boyicup", "Chariup", "Dudijup", "Dwalgan","Meribup", "Murtin", "Tone", "Warrup", "Yerramin", "Yeticup")
trap.points = trap.points %>% filter(Study.AreaName %in% site.subset)

sf::write_sf(trap.points, "Data_Processing/shapefiles/trap.points.shp")

trap.points.jsdm = trap.points %>% 
  st_transform(std.crs) %>%
  dplyr::select(UTM_E, UTM_N,LocationID.x, LocationNa, Study.AreaName, District, Round, Start, End) %>%
  mutate(Start = as_date(Start, format="%d/%m/%y"),
         End = as_date(End, format="%d/%m/%y")) %>%
  rename(X = "UTM_E", Y="UTM_N", LocationID = "LocationID.x", LocationName = "LocationNa")

trap.points.jsdm.buffer400m = trap.points.jsdm %>% st_buffer(dist=400)
trap.points.jsdm.buffer1km = trap.points.jsdm %>% st_buffer(dist=1000)

# Get Dates
dates = unique(trap.points.jsdm$Start)

# Create Mask 
study.region = read_sf("Data_Processing/shapefiles/UpperWarrenMask.shp")
study.region = study.region %>% st_transform(std.crs) %>%
  st_buffer(dist=50) %>% group_by() %>% summarise() %>% st_as_sf()
r = raster(res= 50, extent(st_buffer(study.region, dist=5000)), crs=std.crs$input) # Create raster with same extent and crs as transects

lease = read_sf(here::here("Data_Processing", "shapefiles", "DBCAlease.shp")) %>% st_transform(std.crs) %>% group_by() %>% summarise()
study.region = rbind(study.region, lease)
study.region = study.region %>% group_by() %>% summarise()

base_raster = fasterize(study.region, r, background = 0) # Create base raster, all cells =1
names(base_raster) <- "base_raster"
base_raster_mask = base_raster
base_raster_mask[base_raster_mask==0] <- NA

xm<-matrix(xFromCell(base_raster,c(1:1774680)),nrow=1380,byrow=TRUE)
easting = raster(xm, xmn=417881.2, xmx=482181.2, ymn=6175973, ymx=6244973)
crs(easting) <- crs(base_raster)
names(easting)<-"easting"
ym<-matrix(yFromCell(base_raster,c(1:1774680)),nrow=1380,byrow=TRUE)
northing = raster(ym, xmn=417881.2, xmx=482181.2, ymn=6175973, ymx=6244973)
crs(northing) <- crs(base_raster)
names(northing)<-"northing"


#### Step 2: Create raster buffers for focal area analysis
window_500m = raster::focalWeight(base_raster, d = 500, type='circle')
window_1km = raster::focalWeight(base_raster, d = 1000, type='circle')
window_3km = raster::focalWeight(base_raster, d = 3000, type='circle')
window_5km = raster::focalWeight(base_raster, d = 5000, type='circle')

#### Step 3: Fox Baiting #### 
#  Calculate intensities for each time period
bait_intensity_3yrs_Oct2016 = baiting_lag_rasters(as.Date("1/10/2016", format= "%d/%m/%Y"), lagmonth=38, ras = base_raster)
bait_intensity_3yrs_Oct2017 = baiting_lag_rasters(as.Date("23/10/2017", format= "%d/%m/%Y"), lagmonth=38, ras = base_raster)

#### Step 4: Time Since Fire ####
firehistory = read_sf("~/Library/CloudStorage/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Covariates/FireHistoryDec18_WarrenJarrah.shp")
firehistory = st_transform(firehistory, st_crs(trap.points.jsdm))
trap.dates = dates

# Create TSF raster for start of period
firehistory$TSF = as.numeric(difftime(trap.dates[1],firehistory$FIH_DATE1, unit="weeks"))/52.25
firehistory_temp = subset(firehistory, TSF > 0)
tsf_Oct2016 = fasterize(firehistory_temp, base_raster, field="TSF", fun='min', background=NA)
names(tsf_Oct2016) <- "tsf_Oct2016"
writeRaster(tsf_Oct2016, "Data_Clean/covariate_StudyRegion_tsf.start.Oct16.tif", overwrite=TRUE)

# Create TSF raster for end of period
firehistory$TSF = as.numeric(difftime(trap.dates[11],firehistory$FIH_DATE1, unit="weeks"))/52.25
firehistory_temp = firehistory %>% filter(TSF > 0)
tsf_Oct2017 = fasterize(firehistory_temp, base_raster, field="TSF", fun='min', background=NA)
names(tsf_Oct2017) <- "tsf_Oct2017"
writeRaster(tsf_Oct2017, "Data_Clean/covariate_StudyRegion_tsf.end.Oct17.tif", overwrite=TRUE)

#### Step 5: Fire Severity ####
# Create fire severity rasters for a given date
mga50 <- "+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs"
firehist <- st_read("~/Library/CloudStorage/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Fire Severity prep/For Analysis/firehistory_upperwarren.shp") %>%
  st_transform(std.crs)

dates = firehist %>% dplyr::select(OBJECTID, FIH_DATE1) %>% mutate(OBJECTID = as.character(OBJECTID)) %>% st_drop_geometry()

# Create mask of severity study area
mask = raster("~/Library/CloudStorage/Dropbox/Billy/_research/_PhD/07_UpperWarren_WA/data/Upper Warren Fire Severity/Upper_Warren_GIS/Fire_severity_allClass_sum_Warren_1995to2021_2021-11-08_z50.tif")
mask[mask>0] <- 0

# Fix banding
fireras = raster("~/Library/CloudStorage/Dropbox/Billy/_research/_PhD/07_UpperWarren_WA/data/Upper Warren Fire Severity/CBI_maps/FireSev_Southern Jarrah Forest_675250_qd_L8_dNBRmax.tif")
fire675250 = filter(firehist, OBJECTID==675250) %>% mutate(d=1)
fire675250 = fasterize(fire675250, fireras,field='d', background = NA)
fire675250 = sum(fire675250,fireras,na.rm=TRUE)
fire675250[fire675250==1]<-6
fire675250[fire675250==0]<-NA
fire675250 = fire675250-1


# Max fire severity in the last x years
sev.6yrs.2016 = severity_raster(polyint = study.region, firehist = firehist, reference.date = "1/10/2016", lagyears = 6)
sev.10yrs.2016 = severity_raster(polyint = study.region, firehist = firehist, reference.date = "1/10/2016", lagyears = 10)
sev_6yrs_Oct2016 = resample(sev.6yrs.2016, base_raster, method='ngb'); names(sev_6yrs_Oct2016) <- "sev_6yrs_Oct2016"
sev_10yrs_Oct2016 = resample(sev.10yrs.2016, base_raster, method='ngb'); names(sev_10yrs_Oct2016) <- "sev_10yrs_Oct2016"

sev.6yrs.2017 = severity_raster(polyint = study.region, firehist = firehist, reference.date = "23/10/2017", lagyears = 6)
sev.10yrs.2017 = severity_raster(polyint = study.region, firehist = firehist, reference.date = "23/10/2017", lagyears = 10)
sev_6yrs_Oct2017 = resample(sev.6yrs.2017, base_raster, method='ngb'); names(sev_6yrs_Oct2017) <- "sev_6yrs_Oct2017"
sev_10yrs_Oct2017 = resample(sev.10yrs.2017, base_raster, method='ngb'); names(sev_10yrs_Oct2017) <- "sev_10yrs_Oct2017"

# Proportion of severe fire in the last x years
propsev_6yrs_Oct2016_500m = focal(sev_6yrs_Oct2016, w=window_500m, fun = 'sum'); names(propsev_6yrs_Oct2016_500m) <- "propsev_6yrs_Oct2016_500m"
propsev_10yrs_Oct2016_500m = focal(sev_10yrs_Oct2016, w=window_500m, fun ='sum'); names(propsev_10yrs_Oct2016_500m) <- "propsev_10yrs_Oct2016_500m"
propsev_6yrs_Oct2016_1km = focal(sev_6yrs_Oct2016, w=window_1km, fun ='sum'); names(propsev_6yrs_Oct2016_1km) <- "propsev_6yrs_Oct2016_1km"
propsev_10yrs_Oct2016_1km = focal(sev_10yrs_Oct2016, w=window_1km, fun ='sum'); names(propsev_10yrs_Oct2016_1km) <- "propsev_10yrs_Oct2016_1km"

propsev_6yrs_Oct2017_500m = focal(sev_6yrs_Oct2017, w=window_500m, fun ='sum'); names(propsev_6yrs_Oct2017_500m) <- "propsev_6yrs_Oct2017_500m"
propsev_10yrs_Oct2017_500m = focal(sev_10yrs_Oct2017, w=window_500m, fun ='sum'); names(propsev_10yrs_Oct2017_500m) <- "propsev_10yrs_Oct2017_500m"
propsev_6yrs_Oct2017_1km = focal(sev_6yrs_Oct2017, w=window_1km, fun ='sum'); names(propsev_6yrs_Oct2017_1km) <- "propsev_6yrs_Oct2017_1km"
propsev_10yrs_Oct2017_1km = focal(sev_10yrs_Oct2017, w=window_1km, fun ='sum'); names(propsev_10yrs_Oct2017_1km) <- "propsev_10yrs_Oct2017_1km"

#### Step 6: Severity Count ####
# Count of severe fires in the last x years
sev.count.20yrs.2016 = severe_count_raster(polyint = study.region, firehist = firehist, reference.date = "1/10/2016", lagyears = 20)
sevcount_20yrs_Oct2016 = resample(sev.count.20yrs.2016, base_raster, method='ngb'); names(sevcount_20yrs_Oct2016) <- "sevcount_20yrs_Oct2016"

sev.count.20yrs.2017 = severe_count_raster(polyint = study.region, firehist = firehist, reference.date = "23/10/2017", lagyears = 20)
sevcount_20yrs_Oct2017 = resample(sev.count.20yrs.2017, base_raster, method='ngb'); names(sevcount_20yrs_Oct2017) <- "sevcount_20yrs_Oct2017"

#### Step 7: Native Vegetation ####
nv = read_sf("~/Library/CloudStorage/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Covariates/nativeveg_warrenjarrah.shp")
nv = nv %>% st_transform(std.crs)
nv_raster = fasterize(nv, base_raster, background=0)

prop_nv_500m = focal(nv_raster, w=window_500m, fun ='sum'); names(prop_nv_500m) <- "prop_nv_500m"
prop_nv_1km = focal(nv_raster, w=window_1km, fun ='sum'); names(prop_nv_1km) <- "prop_nv_1km"
prop_nv_3km = focal(nv_raster, w=window_3km, fun ='sum'); names(prop_nv_3km) <- "prop_nv_3km"
prop_nv_5km = focal(nv_raster, w=window_5km, fun ='sum'); names(prop_nv_5km) <- "prop_nv_5km"

#### Step 8: Agricultural Fragmentation ####
managed_tenure = read_sf("~/Library/CloudStorage/Dropbox/Billy/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Covariates/Managed_Tenure_WarrenJarrah.shp")
lease = read_sf(here::here("Data_Processing", "shapefiles", "DBCAlease.shp")) %>% st_transform(std.crs) %>% group_by() %>% summarise()
managed_tenure = managed_tenure %>% st_transform(std.crs) %>% 
  st_buffer(dist=50) %>% group_by() %>% summarise()
managed_tenure = rbind(managed_tenure, lease)
managed_tenure = managed_tenure %>% group_by() %>% summarise()
managed_tenure_raster = fasterize(managed_tenure, base_raster, background = 0)

# Distance to Agricultural Land
region = base_raster
region[region==0] <- 1
ag_land = region - managed_tenure_raster
ag_land[ag_land == 0] <- NA
dist_to_ag = distance(ag_land); names(dist_to_ag) <- "dist_to_ag"

# Proportion of surrounding landscape buffered
ag_land[is.na(ag_land)] <- 0
prop_ag_500m = focal(ag_land, w=window_500m, fun ='sum'); names(prop_ag_500m) <- "prop_ag_500m"
prop_ag_1km = focal(ag_land, w=window_1km, fun ='sum'); names(prop_ag_1km) <- "prop_ag_1km"
prop_ag_3km = focal(ag_land, w=window_3km, fun ='sum'); names(prop_ag_3km) <- "prop_ag_3km"
prop_ag_5km = focal(ag_land, w=window_5km, fun ='sum'); names(prop_ag_5km) <- "prop_ag_5km"

#### Step 9: Distance to Minor Hydo ####
minor_hydro = read_sf(here::here("Data_Processing", "shapefiles",  "UW_Minor_Hydro.shp"))
minor_hydro = minor_hydro %>% st_transform(std.crs) %>% st_buffer(dist=10) %>% group_by() %>% summarise() # Takes a few minutes
minor_hydro_raster = fasterize(minor_hydro, base_raster, background = NA)
dist_to_minorhydro = distance(minor_hydro_raster) # Takes < 1 min
names(dist_to_minorhydro) <- "dist_to_minorhydro"

#### Step 10: Distance to Major Hydro ####
major_hydro = read_sf(here::here("Data_Processing", "shapefiles", "UW_MajorHydro.shp"))
major_hydro = major_hydro %>% st_set_crs("+proj=longlat +datum=WGS84") %>% st_transform(std.crs) %>% st_buffer(dist=10) %>% group_by() %>% summarise() # Takes a few minutes
major_hydro_raster = fasterize(major_hydro, base_raster, background = NA)
dist_to_majorhydro = distance(major_hydro_raster) # Takes < 1 min
names(dist_to_majorhydro) <- "dist_to_majorhydro"

#### Step 11: Distance to All Hydro ###3
all_hydro = rbind(minor_hydro, major_hydro)
all_hydro_raster = fasterize(all_hydro, base_raster, background = NA)
dist_to_allhydro = distance(all_hydro_raster) # Takes < 1 min
names(dist_to_allhydro) <- "dist_to_allhydro"

#### Step 12: Topograpic Wetness Index ####
# Taken from CSIRO : https://data.csiro.au/collection/csiro:5588
twi = raster("~/Library/CloudStorage/Dropbox/Billy/_Research/Data/spatial/twi/twi_3s.tif")
clippoly = study.region %>% st_buffer(dist=5000) %>% st_transform(st_crs(twi)) 
twi = raster::crop(twi, extent(clippoly))
twi = projectRaster(twi, crs=std.crs$input)
twi = resample(twi, base_raster)
names(twi) <- "topographic_wetness_index"

#### Step 13: Road Density ####
roads = read_sf(here::here("Data_Processing", "shapefiles", "UW_all_roads.shp"))
roads_filtered = roads %>% filter(!RDS_RDS_LE %in% c("Track restricted", 
                                                     "Type 4 - track",
                                                     "Unused - to be determined", 
                                                     "6E - track (4x4)",
                                                     "6C - NonDECrdDECmanaged minor road (unsealed) restricted",
                                                     "6B - NonDECrdDECmanaged secondary road (unsealed)"))  %>%
  st_transform(std.crs) %>% st_buffer(dist=50) %>% group_by() %>% summarise() # Takes a few minutes

filtered_roads_raster = fasterize(roads_filtered, base_raster, background = 0)

prop_filtered_roads_1km = focal(filtered_roads_raster, w=window_1km, fun ='sum'); names(prop_filtered_roads_1km) <- "prop_filtered_roads_1km"
prop_filtered_roads_3km = focal(filtered_roads_raster, w=window_3km, fun ='sum'); names(prop_filtered_roads_3km) <- "prop_filtered_roads_3km"
prop_filtered_roads_5km = focal(filtered_roads_raster, w=window_5km, fun ='sum'); names(prop_filtered_roads_5km) <- "prop_filtered_roads_5km"

#### Step 14: Landscape Position and Veg ####
veg = read_sf(here::here("Data_Processing", "shapefiles", "ForestEcosystems_WarrenJarrah.shp"))
veg = veg %>% st_transform(std.crs) %>% mutate(Veg = ifelse(Descript %in% c("Jarrah - North East", 
                                                                            "Jarrah - North West",
                                                                            "Jarrah - Blackwood Plateau",
                                                                            "Jarrah woodland",
                                                                            "Jarrah - Leeuwin Ridge",
                                                                            "Jarrah - Sandy Basins",
                                                                            "Jarrah - South",
                                                                            "Jarrah - Unicup",
                                                                            "Jarrah - Mt Lindesay",
                                                                            "Jarrah/Yellow Tingle",
                                                                            "Jarrah/Rates Tingle",
                                                                            "Jarrah/Red Tingle"), "Jarrah",
                                                            ifelse(Descript %in% c("Western Wandoo forest",
                                                                                   "Western Wandoo woodland"), "Wandoo",
                                                                   ifelse(Descript %in% c("Karri - West Coast", "Karri - South Coast",
                                                                                          "Karri/Yellow Tingle", "Karri/Rates Tingle",
                                                                                          "Karri/Red Tingle", "Karri - Main Belt"), 
                                                                          "Karri",
                                                                          ifelse(Descript %in% c("Rocky outcrops", "Swamps",
                                                                                                 "Bullich and Yate", "Peppermint and coastal heathland"), "Other",
                                                                                 Descript)))))  %>% 
  st_crop(extent(base_raster)) 

veg.lookup = data.frame(Veg = unique(veg$Veg),
                        VegID = c(100, 200, 300, 400, 500))

veg = veg %>% left_join(veg.lookup) %>% as_Spatial() %>% st_as_sf()

veg.ras = fasterize(veg, base_raster, field = "VegID", background=0) # Make raster - values correspond to lookup

## Landscape position
pos = read_sf(here::here("Data_Processing", "shapefiles", "Vegetation_Complexes_WarrenJarrah.shp"))
pos = pos %>% st_transform(std.crs) %>% 
  mutate(SUBCATEGOR = ifelse(Hectares == 38.4371169110, "Uplands", SUBCATEGOR)) %>%
  st_crop(extent(base_raster))%>% group_by(SUBCATEGOR) %>% summarise() %>% 
  mutate(LandscapePos = ifelse(SUBCATEGOR %in% c("(Needs to be reviewed)", "<NA>"), NA, SUBCATEGOR))

pos.lookup = data.frame(SUBCATEGOR = unique(pos$LandscapePos), PosID = c(1, 2, 3, 4, 5))
pos = pos %>% left_join(pos.lookup) %>% as_Spatial() %>% st_as_sf()

landscape.position = fasterize(pos, base_raster, background=0, field = 'PosID')

forest.position = landscape.position + veg.ras

## Reclassify as landscape position but splitting wandoo and jarrah
forest.position[forest.position %in% c(1, 101, 201, 301, 401, 501, 300, 100)] <- 0 # The NA landscape positions
forest.position[forest.position %in% c(2, 102, 202, 302, 402, 502)] <- 1 # Depressions and swamps on uplands
forest.position[forest.position %in% c(3, 103, 203, 303, 403, 503)] <- 2 # Uplands
forest.position[forest.position %in% c(4, 5, 104, 105, 204, 205, 404, 405, 504, 505)] <- 3 # Other valleys, valley floor and swamps
forest.position[forest.position %in% c(304, 305)] <- 4 # Wandoo Valleys and Valley Floor and Swamps 

names(forest.position) <- "forest.position"
names(landscape.position) <- "landscape.position"


#### Step 15: Mean Annual Rainfall (WorldClim Bio12)
mean.rainfall = raster("~/Library/CloudStorage/Dropbox/Billy/_research/data/spatial/WorldClim/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")
clippoly = study.region %>% st_buffer(dist=5000) %>% st_transform(st_crs(mean.rainfall)) 
mean.rainfall = raster::crop(mean.rainfall, extent(clippoly))
mean.rainfall = projectRaster(mean.rainfall, crs=std.crs$input)
mean.rainfall = resample(mean.rainfall, base_raster)
names(mean.rainfall) <- "mean_annual_rainfall"

#### Step 15: Compile into Stack and Save #### 
covariates = stack(base_raster, 
                   # Easting and Northing
                   easting, northing,
                   # Fox Baiting
                   bait_intensity_3yrs_Oct2016, bait_intensity_3yrs_Oct2017,
                   # Time Since Fire
                   tsf_Oct2016, tsf_Oct2017,
                   # Fire Severity
                   propsev_6yrs_Oct2016_500m, propsev_6yrs_Oct2016_1km, propsev_10yrs_Oct2016_500m, propsev_10yrs_Oct2016_1km,
                   propsev_6yrs_Oct2017_500m, propsev_6yrs_Oct2017_1km, propsev_10yrs_Oct2017_500m, propsev_10yrs_Oct2017_1km,
                   # Severe Fires Count
                   sevcount_20yrs_Oct2016, sevcount_20yrs_Oct2017,
                   # Native Veg
                   prop_nv_500m, prop_nv_3km, prop_nv_1km, prop_nv_5km, 
                   # Agriculture
                   dist_to_ag, prop_ag_500m, prop_ag_3km, prop_ag_1km, prop_ag_5km, 
                   # Hydrology
                   dist_to_minorhydro, dist_to_majorhydro, dist_to_allhydro, twi,
                   # Roads
                   prop_filtered_roads_1km, prop_filtered_roads_3km, prop_filtered_roads_5km,
                   # Landscape Position
                   landscape.position, forest.position,
                   # Mean Annual Rainfall
                   mean.rainfall)
# Mask theem
covariates_masked = mask(covariates, base_raster_mask)

# Save the actual raster files as raster stack
covariates_masked = readAll(covariates_masked)
save(covariates_masked, file=here::here("Data_Clean", "covariates_stack.RData"))

#### Step 16: Extract for trap points and dates ####
trap.dates = unique(trap.points.jsdm$Start)
# Create dataframe for covariates
covariates.df = trap.points.jsdm %>% st_drop_geometry()

### Temporally Variable Covariates ###
# Fox Baiting
transect.baiting.out.400 = parallel::mclapply(trap.dates, baiting_lags, polys = trap.points.jsdm.buffer400m, lagmonth=38, ras=base_raster, mc.cores = 3)
baiting.out.400 = do.call('rbind', transect.baiting.out.400)
baiting.out.400$Date = as.Date(baiting.out.400$Date)
names(baiting.out.400)[4]<-"Mean_Intensity_400"

covariates.df = left_join(covariates.df, baiting.out.400, by=c("LocationName", "Start"="Date"))  

# Time Since Fire
# Extract survey specific tsf for each site
fire.list = list()
for (i in 1:length(trap.dates)) {
  start.date = trap.dates[i]
  firehistory$TSF = as.numeric(difftime(start.date,firehistory$FIH_DATE1, unit="weeks"))/52.25
  #time since fire
  firehistory_temp = subset(firehistory, TSF > 0)
  tsf_raster = fasterize(firehistory_temp, mask, field="TSF", fun='min')
  extract_traps = trap.points.jsdm %>% dplyr::filter(Start == start.date)
  extract_traps_buff = trap.points.jsdm.buffer400m %>% dplyr::filter(Start == start.date)
  tsf.point = raster::extract(tsf_raster, extract_traps)
  fire.data = cbind(extract_traps, tsf.point)
  fire.list[[as.character(start.date)]] = fire.data
}

tsf.out = do.call(rbind, fire.list)
covariates.df$tsf.point = ifelse(is.na(tsf.out$tsf.point), max(tsf.out$tsf.point, na.rm=TRUE), tsf.out$tsf.point)

# Fire Severity
sev.6yrs = lapply(trap.dates, prop_severe_raster, polyint = project.area, window_buffer=window_500m, firehist = firehist, lagyears = 6)
sev.6yrs = stack(sev.6yrs)
names(sev.6yrs) = trap.dates
crs(sev.6yrs) <- crs(trap.points.jsdm)

prop_severe_500 = data.frame(raster::extract(sev.6yrs, trap.points))

propsev500.data = trap.points.jsdm %>% st_drop_geometry() %>% 
  dplyr::select(X,Y,LocationID, LocationName) %>% 
  cbind(prop_severe_500) %>%
  pivot_longer(cols = 5:15, names_to="Date",values_to = "PropSevere500") %>% 
  mutate(Date = as.Date((gsub("X", "", Date)), format="%Y.%m.%d"))

covariates.df = left_join(covariates.df, propsev500.data, by=c("X","Y","LocationID","LocationName", "Start"="Date"))  


#### Static Covariates 
covariates.df = covariates.df %>%
  mutate(prop_nv_500m = raster::extract(prop_nv_500m, trap.points),
         prop_nv_1km = raster::extract(prop_nv_1km, trap.points),
         prop_nv_3km = raster::extract(prop_nv_3km, trap.points),
         prop_nv_5km = raster::extract(prop_nv_5km, trap.points),
         dist_to_ag = raster::extract(dist_to_ag, trap.points),
         prop_ag_500m = raster::extract(prop_ag_500m, trap.points),
         prop_ag_1km = raster::extract(prop_ag_1km, trap.points),
         prop_ag_3km = raster::extract(prop_ag_3km, trap.points),
         prop_ag_5km = raster::extract(prop_ag_5km, trap.points),
         dist_to_majorhydro = raster::extract(dist_to_majorhydro, trap.points),
         dist_to_allhydro = raster::extract(dist_to_allhydro, trap.points),
         twi = raster::extract(twi, trap.points),
         prop_filtered_roads_3km = raster::extract(prop_filtered_roads_3km, trap.points), 
         landscape_position = raster::extract(landscape.position, trap.points), 
         forest_position = raster::extract(forest.position, trap.points),
         rainfall = raster::extract(mean.rainfall, trap.points))
    
# Check for colinearity in the covariates
covs = covariates.df %>% 
  dplyr::select(X,Y,11:29)
cor.covs = cor(covs, use='complete.obs', method='pearson')
corrplot::corrplot(cor.covs)
high.cors = cor.covs %>% as.table() %>% as.data.frame() %>% filter(Freq < -0.7 | Freq >0.7) %>% filter(Freq !=1) %>%
  group_by(Var1) %>% summarise(Variables=paste(Var2, collapse=", "))

#### Step 17: Add covariates to detection history df ####
## Occupancy
dethist = readRDS("Data_Processing/camtrapR.dethist.UpperWarren.multispp.RData")
dethist$Sites = left_join(dethist$Sites[,1:9], covariates.df)
saveRDS(dethist, "Data_Processing/camtrapR.dethist.UpperWarren.multispp09122023.RData")
  
## Abundance
counthist = readRDS("Data_Processing/camtrapR.counthist.UpperWarren.multispp.RData")
counthist$Sites = dethist$Sites
saveRDS(counthist, "Data_Processing/camtrapR.counthist.UpperWarren.multispp09122023.RData")




  