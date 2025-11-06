##++++++++++++++++++++++++++++++++++++++++++++++++++##
## 1: Create data frame for Multispp NMixture Model ##
##++++++++++++++++++++++++++++++++++++++++++++++++++##

#### Step 1: Set Up ####
library(ggplot2)
library(dplyr)
library(lubridate)
library(camtrapR)
library(unmarked)

# Read in raw data
eradicat = read.csv("~/Library/CloudStorage/Dropbox/Billy/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Eradicat data/eradicat.camera.data.csv")
ct.table.site = read.csv("~/Library/CloudStorage/Dropbox/Billy/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Eradicat data/camera.operational.csv")
ct.cam.list = read.csv("~/Library/CloudStorage/Dropbox/Billy/_research/_PhD/07_UpperWarren_WA/data/from Adrian Feb19/Eradicat data/Site.List.csv")

# Constrain data to just Upper Warren
site.subset = c("Balban","Boyicup", "Chariup", "Dudijup", "Dwalgan","Meribup", "Murtin", "Tone", "Warrup", "Yerramin", "Yeticup")

#### Step 2: Clean up Data #### 
cams = eradicat %>% rename(Site=Study.AreaName) %>% filter(Site %in% site.subset) %>% group_by(Site, LocationName, LocationID) %>% summarise() 
ct.table.site = ct.table.site %>% filter(Site %in% site.subset)
ct.table = left_join(ct.table.site, cams, by="Site")
ct.table$Start = as.Date(ct.table$Start, format="%d/%m/%y")
ct.table$End = as.Date(ct.table$End, format="%d/%m/%y")
ct.table = ct.table %>% filter(Treatment == 'Ground')

# Create table of when cameras were active
cam.op = cameraOperation(ct.table, stationCol="LocationName", setupCol = "Start", retrievalCol = "End", dateFormat="%Y-%m-%d")
eradicat$DateTimeOriginal = as.POSIXlt(eradicat$Independent.Record.start.date...time, format="%d/%m/%Y %H:%M")
eradicat$Species = eradicat$Species.Name

#### Step 3: Select which sites we want to focus on ####
# We want to filter to just the road transect sites as they are 'more independent' 
eradicat = left_join(eradicat, ct.table.site, by = c('Study.AreaName' = 'Site'))

eradicat.ground = eradicat %>% filter(Treatment == 'Ground') %>% filter(Study.AreaName %in% site.subset)

#### Step 4: Create detection histories ####
focal.species.nights = ct.table 

# Plot some stats to see # of detections
det.stats = eradicat.ground %>% 
  group_by(Species) %>% summarise(DetCount = n(), Sites = length(unique(LocationName)))

det.stats = det.stats %>% filter(DetCount > 50) %>%
  filter(!Species %in% c('Bird sp', 'Small Mammal', 'Unknown'))

ggplot(det.stats) + geom_bar(aes(x=Species, y = DetCount), stat='identity') + 
  geom_text(aes(x=Species, y=DetCount, label=DetCount),vjust=-1) +
  theme_bw() + ylim(0,8500) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(det.stats) + geom_bar(aes(x=Species, y = Sites/550), stat='identity') + 
  geom_text(aes(x=Species, y=Sites/550, label=Sites),vjust=-1) +
  theme_bw() + ylim(0,1) + ylab("Proportion of Sites Detected") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



# Which species do we want to look at?
focal.species = det.stats$Species

# Fix midnight issues by adding one minute
eradicat.ground$DateTimeOriginalFix = eradicat.ground$DateTimeOriginal + lubridate::seconds(1)

eradicat.ground = tidyr::drop_na(eradicat.ground, DateTimeOriginalFix)

# Loop through and create detection histories for each species
# Loop populates a list of detection histories that is ready for JAGS
det.cam = list()
multi.list.cam = list()
for (f in focal.species){
  det = detectionHistory(eradicat.ground, 
                         species = paste(f),
                         camOp = cam.op, 
                         stationCol = "LocationName", 
                         recordDateTimeCol = "DateTimeOriginalFix", 
                         occasionLength = 1, 
                         includeEffort = FALSE,
                         timeZone = 'Australia/Perth',
                         day1="station")
  
  nights = data.frame(rowSums(det$detection_history, na.rm=TRUE))
  nights$Site = row.names(nights)
  colnames(nights) = c(paste(f),"LocationName")
  
  focal.species.nights = merge(focal.species.nights, nights, by="LocationName")
  
  det.cam[[paste(f)]] = det
  
  det.hist = as.matrix(det$detection_history[,2:ncol(det$detection_history)])
  
  multi.list.cam[[paste(f)]] = det.hist
  
}

nsurvey = rowSums(!is.na(multi.list.cam$`Black Rat`), na.rm=TRUE)

UpperWarrenMultiSpp = list(Species = multi.list.cam, 
                           Sites = ct.table,
                           nsurvey = rowSums(!is.na(multi.list.cam[[1]]), na.rm=TRUE))

saveRDS(UpperWarrenMultiSpp, "Data_Processing/camtrapR.dethist.UpperWarren.multispp.RData")

#### Step 5: N Mixture Count Dethist ####
# Loop through and create detection histories for each species
# Loop populates a list of detection histories that is ready for JAGS
det.cam = list()
multi.list.cam = list()
for (f in focal.species){
  det = detectionHistory(eradicat.ground, 
                         species = paste(f),
                         camOp = cam.op, 
                         stationCol = "LocationName", 
                         recordDateTimeCol = "DateTimeOriginalFix", 
                         occasionLength = 1, 
                         output = "count",
                         includeEffort = FALSE,
                         timeZone = 'Australia/Perth',
                         day1="station")
  
  nights = data.frame(rowSums(det$detection_history, na.rm=TRUE))
  nights$Site = row.names(nights)
  colnames(nights) = c(paste(f),"LocationName")
  
  focal.species.nights = merge(focal.species.nights, nights, by="LocationName")
  
  det.cam[[paste(f)]] = det
  
  det.hist = as.matrix(det$detection_history[,2:ncol(det$detection_history)])
  
  multi.list.cam[[paste(f)]] = det.hist
  
}

UpperWarrenMultiSpp = list(Species = multi.list.cam, 
                           Sites = ct.table,
                           nsurvey = rowSums(!is.na(multi.list.cam[[1]]), na.rm=TRUE))

saveRDS(UpperWarrenMultiSpp, "Data_Processing/camtrapR.counthist.UpperWarren.multispp.RData")



