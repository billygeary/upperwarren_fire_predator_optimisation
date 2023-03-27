#### Script containing covariate functions ####
#### Fox Baiting Input Data ####
bait.crs = st_crs(20350)
# Read in baiting spatial data
# First do Manjimup Western Shield ground baiting
balban = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Baiting/Balban.shp")
camelar = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Baiting/Camelar.shp")
keninup = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Baiting/Keninup.shp")
kingston = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Baiting/Kingston.shp")
perupcore = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Baiting/Perup Core.shp")
perup_yendicup = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Baiting/Perup_Yendicup_Baiting.shp")

perupcore = perupcore  %>% st_union() %>% st_as_sf()
balban = balban %>% st_transform(bait.crs) %>% st_union() %>% st_buffer(dist=2500) %>% st_as_sf() 
camelar = camelar %>% st_set_crs(bait.crs) %>% st_union() %>% st_buffer(dist=2500) %>% st_as_sf() 
keninup = keninup %>% st_set_crs(bait.crs) %>% st_union() %>% st_buffer(dist=2500) %>% st_as_sf() 
kingston = kingston %>% st_set_crs(bait.crs) %>% st_union() %>% st_buffer(dist=2500) %>% st_as_sf() 
perupcore = perupcore %>% st_set_crs(bait.crs) %>% st_union() %>% st_buffer(dist=2500) %>% st_as_sf() 
yendicup = perup_yendicup %>% filter(Program == "Yendicup Core Baiting") %>% st_union() %>% st_buffer(dist=2500) %>% st_as_sf() 
# Perup whole is an amalgamation of balban, camelar and keninup where we cant distinguish the bait numbers between transects
perup_whole = st_union(balban, camelar) 
perup_whole = st_union(perup_whole, keninup)

# Add in extra Western Shield ground baiting (Denbarker & Lake Muir)
den.muir = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Predator data/BaitingData Request/baiting data from Marika 2021/Baiting History/Denbarker.Lake Muir 15.08.19.shp")
den.muir= den.muir %>% st_set_crs(bait.crs) %>% st_union()  %>% st_buffer(dist=2500) %>% st_as_sf() 

# Add in De Longraft road pre-October 2012
# Discussed with AW and MM on 20/04/2021
# Quarterly baiting at 50 baits per quarter estimate
delandgraaft = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Predator data/BaitingData Request/baiting data from Marika 2021/Baiting History/deLandgrafft.shp")
delandgraaft = delandgraaft %>% st_transform(bait.crs) %>% st_union() %>% st_buffer(dist=2500) %>% st_as_sf() 

# Add in Translocations
translocations = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Predator data/BaitingData Request/baiting data from Marika 2021/Baiting History/Translocations.shp")
translocations = translocations %>% st_transform(bait.crs) %>% st_buffer(dist=2500) %>% st_as_sf() 
warrup_trans = translocations %>% filter(Zone =="Warrup") %>% st_union() %>% st_as_sf()
yendicup_trans = translocations %>% filter(Zone =="Yendicup") %>% st_union() %>% st_as_sf()
walcott_trans = translocations %>% filter(Zone =="Walcott") %>% st_union() %>% st_as_sf()

# Add in 2021 post fire baiting efforts
corbal_pf = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Predator data/BaitingData Request/baiting data from Marika 2021/Baiting History/CorbalPostFireBaiting2021_22.shp")
corbal_pf = corbal_pf %>% st_set_crs(std.crs) %>% st_transform(bait.crs) %>% st_buffer(dist=2500) %>% st_as_sf() 

weinup_pf = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Predator data/BaitingData Request/baiting data from Marika 2021/Baiting History/WeinupPostFireBaiting2021.shp")
weinup_pf =  weinup_pf %>% st_set_crs(std.crs) %>% st_transform(bait.crs) %>% st_buffer(dist=2500) %>% st_as_sf() 

# Add in FPC Baiting
coonan4 = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Predator data/BaitingData Request/FPC/FPC_Coonan4Baiting.shp")
coonan4 =  coonan4 %>% st_transform(bait.crs) %>% st_buffer(dist=2500) %>% st_as_sf() 

coonan5 = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Predator data/BaitingData Request/FPC/FPC_Coonan5Baiting.shp")
coonan5 =  coonan5 %>% st_transform(bait.crs) %>% st_buffer(dist=2500) %>% st_as_sf() 

yeticup = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Predator data/BaitingData Request/FPC/FPC_YeticupBaiting.shp")
yeticup =  yeticup %>% st_transform(bait.crs) %>% st_buffer(dist=2500) %>% st_as_sf() 

# Last do Western Shield Aerial baiting for three cells
aerial = read_sf("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Spatial data/Baiting/WesternShield Flight Cells.shp")
manjimup = aerial %>% st_transform(bait.crs) %>% filter(CELL_NAME == "Manjimup") %>% st_union() %>% st_as_sf() 
shannon = aerial %>% st_transform(bait.crs) %>% filter(CELL_NAME == "Shannon") %>% st_union() %>% st_as_sf() 
denbarker = aerial %>% st_transform(bait.crs) %>% filter(CELL_NAME == "Denbarker") %>% st_union() %>% st_as_sf() 


# Baiting count data across whole study region
baits = read.csv("~/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Predator data/predator.baiting.combined_BG_10052022.csv")
baits = baits %>% filter(Purpose != "FPC Fox Baiting") %>%
  mutate(StartDate = as.Date(StartDate, format= "%d/%m/%y"))

# Check which months baiting happens to make sure we include key baiting runs
bait.summary = baits %>% 
  mutate(StartDate = as.Date(StartDate, format= "%d/%m/%y")) %>%
  mutate(Year = lubridate::year(StartDate), Month = lubridate::month(StartDate), FirstMonth = lubridate::floor_date(StartDate, 'month')) %>%
  filter(Year %in% c(2016:2022)) %>% 
  group_by(Year, FirstMonth, Aerial_Ground) %>% summarise(BaitCount = sum(BaitCount))

aerial.baits = baits %>% filter(Aerial_Ground == 'Aerial') %>% 
  mutate(Year = lubridate::year(StartDate)) %>%
  group_by(BaitingCell, Year) %>%
  summarise(Runs = n(), 
            TotalBaits = sum(BaitCount))

#### Fox Baiting Rasters ####
baiting_lag_rasters <- function(date, lagmonth, ras){
  date.num = as.Date(date)
  # Get baiting events within tme period
  lag_baits = baits %>% filter(StartDate < date.num  & StartDate > (date.num-months(lagmonth)))  # Baiting runs in previous lag months
  
  # Aerial baiting
  manjimup_area = manjimup %>% st_area() %>% as.numeric()
  manjimup_lag = lag_baits %>% filter(Zone == "Manjimup Cell")
  manjimup$lag = sum(manjimup_lag$BaitCount)/(manjimup_area/1000000) # Baits per km in previous 3 months
  manjimup_lag_raster = fasterize(manjimup, ras, field="lag")
  
  # Shannon
  shannon_area = shannon %>% st_area() %>% as.numeric()
  shannon_lag = lag_baits %>% filter(Zone == "Shannon")
  shannon$lag = sum(shannon_lag$BaitCount)/(shannon_area/1000000) # Baits per km in previous 3 months
  shannon_lag_raster = fasterize(shannon, ras, field="lag")
  
  # Denbarker
  denbarker_area = denbarker %>% st_area() %>% as.numeric()
  denbarker_lag = lag_baits %>% filter(Zone == "Denbarker")
  denbarker$lag = sum(denbarker_lag$BaitCount)/(denbarker_area/1000000) # Baits per km in previous 3 months
  denbarker_lag_raster = fasterize(denbarker, ras, field="lag")
  
  # Kingston
  kingston_area = kingston %>% st_area() %>% as.numeric()
  kingston_lag = lag_baits %>% filter(Zone == "Kingston") 
  kingston$lag = sum(kingston_lag$BaitCount)/(kingston_area/1000000)
  kingston_lag_raster = fasterize(kingston, ras, field="lag")
  
  # Perup (whole) -- only pre-2010
  perup_area = perup_whole %>% st_area() %>% as.numeric()
  perup_lag = lag_baits %>% filter(Zone == "Perup ground (Ken, Cam, Bal)")
  perup_whole$lag = sum(perup_lag$BaitCount)/(perup_area/1000000)
  perup_lag_raster = fasterize(perup_whole, ras, field="lag")
  
  # De landgraaft
  delandgraaft_area = delandgraaft %>% st_area() %>% as.numeric()
  delandgraaft_lag = lag_baits %>% filter(Zone == "deLandgraaft") 
  delandgraaft$lag = sum(delandgraaft_lag$BaitCount)/(delandgraaft_area/1000000)
  delandgraaft_lag_raster = fasterize(delandgraaft, ras, field="lag")
  
  # Perup Core -- only post=
  perupcore_area = perupcore %>% st_area() %>% as.numeric()
  perupcore_lag = lag_baits %>% filter(Zone == "PerupCore") 
  perupcore$lag = sum(perupcore_lag$BaitCount)/(perupcore_area/1000000)
  perupcore_lag_raster = fasterize(perupcore, ras, field="lag")
  
  # Balban
  balban_area = balban %>% st_area() %>% as.numeric()
  balban_lag = lag_baits %>% filter(Zone == "Balban") 
  balban$lag= sum(balban_lag$BaitCount)/(balban_area/1000000)
  balban_lag_raster = fasterize(balban, ras, field="lag")
  
  # Camelar
  camelar_area = camelar %>% st_area() %>% as.numeric()
  camelar_lag = lag_baits %>% filter(Zone == "Camelar")  
  camelar$lag = sum(camelar_lag$BaitCount)/(camelar_area/1000000)
  camelar_lag_raster = fasterize(camelar, ras, field="lag")
  
  # Keninuo
  keninup_area = keninup %>% st_area() %>% as.numeric()
  keninup_lag = lag_baits %>% filter(Zone == "Keninup") 
  keninup$lag = sum(keninup_lag$BaitCount)/(keninup_area/1000000)
  keninup_lag_raster = fasterize(keninup, ras, field="lag")
  
  # Denbarker - Lake Muir Ground Baiting
  den.muir_area = den.muir %>% st_area() %>% as.numeric()
  den.muir_lag = lag_baits %>% filter(BaitingCell == "Denbarker" & Purpose == "Western Shield" & Aerial_Ground == "Ground") 
  den.muir$lag = sum(den.muir_lag$BaitCount)/(den.muir_area/1000000)
  den.muir_lag_raster = fasterize(den.muir, ras, field="lag")
  
  # Translocations
  ## Warrup
  warrup_trans_area = warrup_trans %>% st_area() %>% as.numeric()
  warrup_trans_lag = lag_baits %>% filter(Zone == "Warrup" & Purpose == "Translocation") 
  warrup_trans$lag = sum(warrup_trans_lag$BaitCount)/(warrup_trans_area/1000000)
  warrup_trans_lag_raster = fasterize(warrup_trans, ras, field="lag")
  
  ## Yendicup
  yendicup_trans_area = yendicup_trans %>% st_area() %>% as.numeric()
  yendicup_trans_lag = lag_baits %>% filter(Zone == "Yendicup" & Purpose == "Translocation") 
  yendicup_trans$lag = sum(yendicup_trans_lag$BaitCount)/(yendicup_trans_area/1000000)
  yendicup_trans_lag_raster = fasterize(yendicup_trans, ras, field="lag")
  
  ## Walcott
  walcott_trans_area = walcott_trans %>% st_area() %>% as.numeric()
  walcott_trans_lag = lag_baits %>% filter(Zone == "Walcott" & Purpose == "Translocation") 
  walcott_trans$lag = sum(walcott_trans_lag$BaitCount)/(walcott_trans_area/1000000)
  walcott_trans_lag_raster = fasterize(walcott_trans, ras, field="lag")
  
  ## Post Fire Baiting
  corbal_pf_area = corbal_pf %>% st_area() %>% as.numeric()
  corbal_pf_lag = lag_baits %>% filter(Zone == "Corbal_PostFire")
  corbal_pf$lag = sum(corbal_pf_lag$BaitCount)/(corbal_pf_area/1000000)
  corbal_pf_lag_raster = fasterize(corbal_pf, ras, field="lag")
  
  weinup_pf_area = weinup_pf %>% st_area() %>% as.numeric()
  weinup_pf_lag = lag_baits %>% filter(Zone == "Weinup_PostFire")
  weinup_pf$lag = sum(weinup_pf_lag$BaitCount)/(weinup_pf_area/1000000)
  weinup_pf_lag_raster = fasterize(weinup_pf, ras, field="lag")
  
  ## FPC Baiting
  coonan4_area = coonan4 %>% st_area() %>% as.numeric()
  coonan4_lag = lag_baits %>% filter(Zone == "FPC_Coonan4Baiting")
  coonan4$lag = sum(coonan4_lag$BaitCount)/(coonan4_area/1000000)
  coonan4_lag_raster = fasterize(coonan4, ras, field="lag")
  
  coonan5_area = coonan5 %>% st_area() %>% as.numeric()
  coonan5_lag = lag_baits %>% filter(Zone == "FPC_Coonan5Baiting")
  coonan5$lag = sum(coonan5_lag$BaitCount)/(coonan5_area/1000000)
  coonan5_lag_raster = fasterize(coonan5, ras, field="lag")
  
  yeticup_area = yeticup %>% st_area() %>% as.numeric()
  yeticup_lag = lag_baits %>% filter(Zone == "FPC_YeticupBaiting")
  yeticup$lag = sum(yeticup_lag$BaitCount)/(yeticup_area/1000000)
  yeticup_lag_raster = fasterize(yeticup, ras, field="lag")
  
  ## Sum them up across region
  sum_lag = sum(manjimup_lag_raster, shannon_lag_raster, denbarker_lag_raster, delandgraaft_lag_raster,
                perup_lag_raster, perupcore_lag_raster, balban_lag_raster, kingston_lag_raster,
                balban_lag_raster, camelar_lag_raster, keninup_lag_raster, den.muir_lag_raster,
                walcott_trans_lag_raster, yendicup_trans_lag_raster, warrup_trans_lag_raster,
                coonan4_lag_raster, coonan5_lag_raster, yeticup_lag_raster,
                corbal_pf_lag_raster, weinup_pf_lag_raster, na.rm=TRUE)
  
  sum_lag[sum_lag==0]<-NA
  # Raster out
  bait.intensity = sum_lag
  names(bait.intensity) <- paste0("bait_intensity_", date.num, "_", lagmonth,"month_lag")
  
  return(bait.intensity)
}

#### Severity Raster ####
severity_raster = function(polyint,firehist, reference.date, lagyears){
  # Set up analysis inputs
  shp = firehist 
  ply = study.region
  d = as.Date(reference.date, format="%d/%m/%Y")
  y = lagyears
  
  # Read in severity layers
  tlist <- data.frame(name = list.files("~/Dropbox/Billy/_Research/Upper Warren Fire Severity/CBI_maps", pattern = "tif"),
                      path = list.files("~/Dropbox/Billy/_Research/Upper Warren Fire Severity/CBI_maps", pattern = "tif", full=TRUE))
  tlist <- tlist %>% 
    mutate(fnum = str_split_fixed(name, "_", 7)[,3]) %>%
    left_join(dates, by=c('fnum'='OBJECTID'))
  
  # Set up output df
  df.out <- data.frame()
  
  pti <- ply # Select transect
  # Find the fire severity rasters we need
  int <- st_intersection(pti, shp) # All fires within transect  
  tifi <- tlist %>% 
    filter(fnum %in% int$OBJECTID) %>% # Get the relevant fires that intersect
    filter(FIH_DATE1 < d & FIH_DATE1 > d-lubridate::years(y)) # Get the relvant fires within moving window time period
  
  if(length(tifi$path) > 0){
    # Create binary map of severity within five years
    t <- raster(tifi$path[1])
    # Check CRS consistency
    if(as.character(crs(t)) == mga50){t=t} else{t=projectRaster(t, crs = mga50, method='ngb')}
    # Check resolution 
    if(res(t)[1] == res(mask)[1]){t=t} else{t = resample(t, mask, method='ngb')}
    
    t <- mosaic(mask, t, fun=sum)
    
    if(length(tifi$path) > 1){
      # If 2 or more fires, loop through all severity layers and create a mosaic of cells severely burnt
      for (i in 2:length(tifi$path)){
        if(tifi$fnum[i]== "675250"){
          ti=fire675250 # Fill holes in 675250
        }else{ti=raster(tifi$path[i])}
        # Check CRS consistency
        if(as.character(crs(ti)) == mga50){ti=ti} else{ti=projectRaster(ti, crs = mga50, method='ngb')}
        # Check resolution 
        if(res(ti)[1] == res(mask)[1]){ti=ti} else{ti = resample(ti, mask, method='ngb')}
        
        t <- mosaic(t, ti, fun = max)
      }
    }else{t=t}
  }else{t=t}
  
  # Raster Out
  t[t<4]<-0
  t[t>3]<-1
  t
}

#### Proportion Severe Raster ####
prop_severe_raster = function(polyint, window_buffer, firehist, reference.date, lagyears){
  # Set up analysis inputs
  shp = firehist 
  ply = study.region
  d = as.Date(reference.date, format="%d/%m/%Y")
  y = lagyears
  
  # Read in severity layers
  tlist <- data.frame(name = list.files("~/Dropbox/Billy/_Research/Upper Warren Fire Severity/CBI_maps", pattern = "tif"),
                      path = list.files("~/Dropbox/Billy/_Research/Upper Warren Fire Severity/CBI_maps", pattern = "tif", full=TRUE))
  tlist <- tlist %>% 
    mutate(fnum = str_split_fixed(name, "_", 7)[,3]) %>%
    left_join(dates, by=c('fnum'='OBJECTID'))
  
  # Set up output df
  df.out <- data.frame()
  
  pti <- ply # Select transect
  # Find the fire severity rasters we need
  int <- st_intersection(pti, shp) # All fires within transect  
  tifi <- tlist %>% 
    filter(fnum %in% int$OBJECTID) %>% # Get the relevant fires that intersect
    filter(FIH_DATE1 < d & FIH_DATE1 > d-lubridate::years(y)) # Get the relvant fires within moving window time period
  
  if(length(tifi$path) > 0){
    # Create binary map of severity within five years
    t <- raster(tifi$path[1])
    # Check CRS consistency
    if(as.character(crs(t)) == mga50){t=t} else{t=projectRaster(t, crs = mga50, method='ngb')}
    # Check resolution 
    if(res(t)[1] == res(mask)[1]){t=t} else{t = resample(t, mask, method='ngb')}
    
    t <- mosaic(mask, t, fun=sum)
    
    if(length(tifi$path) > 1){
      # If 2 or more fires, loop through all severity layers and create a mosaic of cells severely burnt
      for (i in 2:length(tifi$path)){
        if(tifi$fnum[i]== "675250"){
          ti=fire675250 # Fill holes in 675250
        }else{ti=raster(tifi$path[i])}
        # Check CRS consistency
        if(as.character(crs(ti)) == mga50){ti=ti} else{ti=projectRaster(ti, crs = mga50, method='ngb')}
        # Check resolution 
        if(res(ti)[1] == res(mask)[1]){ti=ti} else{ti = resample(ti, mask, method='ngb')}
        
        t <- mosaic(t, ti, fun = max)
      }
    }else{t=t}
  }else{t=t}
  
  # Raster Out
  t[t<4]<-0
  t[t>3]<-1
  t = focal(t, window_buffer, fun='sum')
}

#### Severe Count Raster ####
severe_count_raster = function(polyint,firehist, reference.date, lagyears){
  # Set up analysis inputs
  shp = firehist
  ply = study.region
  d = as.Date(reference.date, format="%d/%m/%Y")
  y = lagyears
  
  # Read in severity layers
  tlist <- data.frame(name = list.files("~/Dropbox/Billy/_Research/Upper Warren Fire Severity/CBI_maps", pattern = "tif"),
                      path = list.files("~/Dropbox/Billy/_Research/Upper Warren Fire Severity/CBI_maps", pattern = "tif", full=TRUE))
  tlist <- tlist %>% 
    mutate(fnum = str_split_fixed(name, "_", 7)[,3]) %>%
    left_join(dates, by=c('fnum'='OBJECTID'))
  
  # Set up output df
  df.out <- data.frame()
  
  pti <- ply # Select transect
  # Find the fire severity rasters we need
  int <- st_intersection(pti, shp) # All fires within transect  
  tifi <- tlist %>% 
    filter(fnum %in% int$OBJECTID) %>% # Get the relevant fires that intersect
    filter(FIH_DATE1 < d & FIH_DATE1 > d-lubridate::years(y)) # Get the relvant fires within moving window time period
  
  if(length(tifi$path) > 0){
    # Create binary map of severity within five years
    t <- raster(tifi$path[1])
    # Check CRS consistency
    if(as.character(crs(t)) == mga50){t=t} else{t=projectRaster(t, crs = mga50, method='ngb')}
    # Check resolution 
    if(res(t)[1] == res(mask)[1]){t=t} else{t = resample(t, mask, method='ngb')}
    
    t <- mosaic(mask, t, fun=sum)
    t[t<4] <- 0
    t[t>3] <- 1
    
    if(length(tifi$path) > 1){
      # If 2 or more fires, loop through all severity layers and create a mosaic of cells severely burnt
      for (i in 2:length(tifi$path)){
        if(tifi$fnum[i]== "675250"){
          ti=fire675250 # Fill holes in 675250
        }else{ti=raster(tifi$path[i])}
        # Check CRS consistency
        if(as.character(crs(ti)) == mga50){ti=ti} else{ti=projectRaster(ti, crs = mga50, method='ngb')}
        # Check resolution 
        if(res(ti)[1] == res(mask)[1]){ti=ti} else{ti = resample(ti, mask, method='ngb')}
        
        ti[ti<4] <- 0
        ti[ti>3] <- 1 
        t <- mosaic(t, ti, fun = sum)
      }
    }else{t=t}
  }else{t=t}
  
  # Raster Out
  t
}

#### Baiting Lags Extraction ####
baiting_lags <- function(date, polys, lagmonth, ras){
  date.num = date
  
  # Get baiting events within tme period
  lag_baits = baits %>% filter(StartDate < date.num  & StartDate > (date.num-months(lagmonth)))  # Baiting runs in previous lag months
  
  # Aerial baiting
  manjimup_area = manjimup %>% st_area() %>% as.numeric()
  manjimup_lag = lag_baits %>% filter(Zone == "Manjimup Cell")
  
  manjimup$lag = sum(manjimup_lag$BaitCount)/(manjimup_area/1000000) # Baits per km in previous 3 months
  manjimup_lag_raster = fasterize(manjimup, ras, field="lag")
  
  # Shannon
  shannon_area = shannon %>% st_area() %>% as.numeric()
  shannon_lag = lag_baits %>% filter(Zone == "Shannon")
  
  shannon$lag = sum(shannon_lag$BaitCount)/(shannon_area/1000000) # Baits per km in previous 3 months
  shannon_lag_raster = fasterize(shannon, ras, field="lag")
  
  # Denbarker
  denbarker_area = denbarker %>% st_area() %>% as.numeric()
  denbarker_lag = lag_baits %>% filter(Zone == "Denbarker")
  
  denbarker$lag = sum(denbarker_lag$BaitCount)/(denbarker_area/1000000) # Baits per km in previous 3 months
  denbarker_lag_raster = fasterize(denbarker, ras, field="lag")
  
  # Kingston
  kingston_area = kingston %>% st_area() %>% as.numeric()
  kingston_lag = lag_baits %>% filter(Zone == "Kingston") 
  
  kingston$lag = sum(kingston_lag$BaitCount)/(kingston_area/1000000)
  kingston_lag_raster = fasterize(kingston, ras, field="lag")
  
  # Perup (whole) -- only pre-2010
  perup_area = perup_whole %>% st_area() %>% as.numeric()
  perup_lag = lag_baits %>% filter(Zone == "Perup ground (Ken, Cam, Bal)")
  
  perup_whole$lag = sum(perup_lag$BaitCount)/(perup_area/1000000)
  perup_lag_raster = fasterize(perup_whole, ras, field="lag")
  
  # De landgraaft
  delandgraaft_area = delandgraaft %>% st_area() %>% as.numeric()
  delandgraaft_lag = lag_baits %>% filter(Zone == "deLandgraaft") 
  
  delandgraaft$lag = sum(delandgraaft_lag$BaitCount)/(delandgraaft_area/1000000)
  
  delandgraaft_lag_raster = fasterize(delandgraaft, ras, field="lag")
  
  # Perup Core -- only post=
  perupcore_area = perupcore %>% st_area() %>% as.numeric()
  perupcore_lag = lag_baits %>% filter(Zone == "PerupCore") 
  
  perupcore$lag = sum(perupcore_lag$BaitCount)/(perupcore_area/1000000)
  
  perupcore_lag_raster = fasterize(perupcore, ras, field="lag")
  
  # Balban
  balban_area = balban %>% st_area() %>% as.numeric()
  balban_lag = lag_baits %>% filter(Zone == "Balban") 
  
  balban$lag= sum(balban_lag$BaitCount)/(balban_area/1000000)
  balban_lag_raster = fasterize(balban, ras, field="lag")
  
  # Camelar
  camelar_area = camelar %>% st_area() %>% as.numeric()
  camelar_lag = lag_baits %>% filter(Zone == "Camelar")  
  
  camelar$lag = sum(camelar_lag$BaitCount)/(camelar_area/1000000)
  
  camelar_lag_raster = fasterize(camelar, ras, field="lag")
  
  # Keninuo
  keninup_area = keninup %>% st_area() %>% as.numeric()
  keninup_lag = lag_baits %>% filter(Zone == "Keninup") 
  
  keninup$lag = sum(keninup_lag$BaitCount)/(keninup_area/1000000)
  keninup_lag_raster = fasterize(keninup, ras, field="lag")
  
  # Denbarker - Lake Muir Ground Baiting
  den.muir_area = den.muir %>% st_area() %>% as.numeric()
  den.muir_lag = lag_baits %>% filter(BaitingCell == "Denbarker" & Purpose == "Western Shield" & Aerial_Ground == "Ground") 
  
  den.muir$lag = sum(den.muir_lag$BaitCount)/(den.muir_area/1000000)
  den.muir_lag_raster = fasterize(den.muir, ras, field="lag")
  
  # Translocations
  ## Warrup
  warrup_trans_area = warrup_trans %>% st_area() %>% as.numeric()
  warrup_trans_lag = lag_baits %>% filter(Zone == "Warrup" & Purpose == "Translocation") 
  
  warrup_trans$lag = sum(warrup_trans_lag$BaitCount)/(warrup_trans_area/1000000)
  warrup_trans_lag_raster = fasterize(warrup_trans, ras, field="lag")
  
  ## Yendicup
  yendicup_trans_area = yendicup_trans %>% st_area() %>% as.numeric()
  yendicup_trans_lag = lag_baits %>% filter(Zone == "Yendicup" & Purpose == "Translocation") 
  
  yendicup_trans$lag = sum(yendicup_trans_lag$BaitCount)/(yendicup_trans_area/1000000)
  yendicup_trans_lag_raster = fasterize(yendicup_trans, ras, field="lag")
  
  ## Walcott
  walcott_trans_area = walcott_trans %>% st_area() %>% as.numeric()
  walcott_trans_lag = lag_baits %>% filter(Zone == "Walcott" & Purpose == "Translocation") 
  
  walcott_trans$lag = sum(walcott_trans_lag$BaitCount)/(walcott_trans_area/1000000)
  
  walcott_trans_lag_raster = fasterize(walcott_trans, ras, field="lag")
  
  ## Post Fire Baiting
  corbal_pf_area = corbal_pf %>% st_area() %>% as.numeric()
  corbal_pf_lag = lag_baits %>% filter(Zone == "Corbal_PostFire")
  corbal_pf$lag = sum(corbal_pf_lag$BaitCount)/(corbal_pf_area/1000000)
  corbal_pf_lag_raster = fasterize(corbal_pf, ras, field="lag")
  
  weinup_pf_area = weinup_pf %>% st_area() %>% as.numeric()
  weinup_pf_lag = lag_baits %>% filter(Zone == "Weinup_PostFire")
  weinup_pf$lag = sum(weinup_pf_lag$BaitCount)/(weinup_pf_area/1000000)
  weinup_pf_lag_raster = fasterize(weinup_pf, ras, field="lag")
  
  ## FPC Baiting
  coonan4_area = coonan4 %>% st_area() %>% as.numeric()
  coonan4_lag = lag_baits %>% filter(Zone == "FPC_Coonan4Baiting")
  coonan4$lag = sum(coonan4_lag$BaitCount)/(coonan4_area/1000000)
  coonan4_lag_raster = fasterize(coonan4, ras, field="lag")
  
  coonan5_area = coonan5 %>% st_area() %>% as.numeric()
  coonan5_lag = lag_baits %>% filter(Zone == "FPC_Coonan5Baiting")
  coonan5$lag = sum(coonan5_lag$BaitCount)/(coonan5_area/1000000)
  coonan5_lag_raster = fasterize(coonan5, ras, field="lag")
  
  yeticup_area = yeticup %>% st_area() %>% as.numeric()
  yeticup_lag = lag_baits %>% filter(Zone == "FPC_YeticupBaiting")
  yeticup$lag = sum(yeticup_lag$BaitCount)/(yeticup_area/1000000)
  yeticup_lag_raster = fasterize(yeticup, ras, field="lag")
  
  ## Sum them up across region
  sum_lag = sum(manjimup_lag_raster, shannon_lag_raster, denbarker_lag_raster, delandgraaft_lag_raster,
                perup_lag_raster, perupcore_lag_raster, balban_lag_raster, kingston_lag_raster,
                balban_lag_raster, camelar_lag_raster, keninup_lag_raster, den.muir_lag_raster,
                walcott_trans_lag_raster, yendicup_trans_lag_raster, warrup_trans_lag_raster,
                coonan4_lag_raster, coonan5_lag_raster, yeticup_lag_raster,
                corbal_pf_lag_raster, weinup_pf_lag_raster, na.rm=TRUE)
  
  # Extract mean baiting intensity for each polygon
  lag_extract = data.frame(Date = paste(date.num),
                           Lag = paste(lagmonth), 
                           LocationName = st_drop_geometry(polys)$LocationName,
                           Mean_Intensity = raster::extract(sum_lag, polys, fun=mean))
  
  baiting.intensity.list = lag_extract
  baiting.intensity.list
}
