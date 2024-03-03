
# Read in data
dethist = readRDS("Data_Processing/camtrapR.counthist.UpperWarren.multispp.RData")

# Select species we want to include in the Multi Spp Occ model
names(dethist$Species)
dethist$Species =  dethist$Species[c("Woylie", "Chuditch", "Koomal", "Quenda", "Roo",
                                     "Vulpes", "Tammar", "Numbat", "Western Brush Wallaby", "Dunnart")]
spp = names(dethist$Species)

sites = dethist$Sites
species = dethist$Species
# Use apply to sum along the rows
detection_table <- data.frame(
  Site = rownames(species$Woylie),
  sapply(species, function(x) apply(x, 1, sum, na.rm=TRUE))
)

joined = dplyr::left_join(sites, detection_table, by =c("LocationName"= "Site"))

write.csv(joined, "Data_Processing/reporting_rates.csv")
