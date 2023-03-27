
### Make a pretty map for inclusion in publication

library(raster)
library(ggplot2)
library(sf)
load(here::here("Data_Clean", "covariates_stack.RData"))
sites = read_sf("Data_Processing/shapefiles/trap.points.shp")

study.location = covariates_masked$base_raster  %>% as.data.frame(xy=TRUE) %>% na.omit()
base_map = ggplot() + 
  geom_raster(data = study.location, aes(x=x,y=y,fill=base_raster)) + scale_fill_manual(palette = "grey") +
  geom_sf(data=sites) +
  theme_void() +
  theme(legend.position = "bottom")

base_map

tsf = covariates_masked$tsf_Oct2016 %>% as.data.frame(xy=TRUE) %>% na.omit()
bait = covariates_masked$bait_intensity_2016.10.01_38month_lag %>% as.data.frame(xy=TRUE) %>% na.omit()
sev = covariates_masked$propsev_6yrs_Oct2016_500m %>% as.data.frame(xy=TRUE) %>% na.omit()
rain = covariates_masked$mean_annual_rainfall %>% as.data.frame(xy=TRUE) %>% na.omit()
twi = covariates_masked$topographic_wetness_index %>% as.data.frame(xy=TRUE) %>% na.omit()

tsf_map = ggplot(tsf) + geom_raster(aes(x=x,y=y, fill=tsf_Oct2016)) + 
            scale_fill_viridis_c(name="Time Since Fire") +
            theme_void() +
            theme(legend.position = "bottom")

bait_map = ggplot(bait) + geom_raster(aes(x=x,y=y, fill=bait_intensity_2016.10.01_38month_lag)) + 
  scale_fill_viridis_c(name="Bait Intensity") +
  theme_void() +
  theme(legend.position = "bottom")

sev_map = ggplot(sev) + geom_raster(aes(x=x,y=y, fill=propsev_6yrs_Oct2016_500m)) + 
  scale_fill_viridis_c(name="Proportion Severely Burnt") +
  theme_void() +
  theme(legend.position = "bottom")

rain_map = ggplot(rain) + geom_raster(aes(x=x,y=y, fill=mean_annual_rainfall)) + 
  scale_fill_viridis_c(name="Mean Annual Rainfall") +
  theme_void() +
  theme(legend.position = "bottom")

twi_map = ggplot(twi) + geom_raster(aes(x=x,y=y, fill=topographic_wetness_index)) + 
  scale_fill_viridis_c(name="Topographic Wetness Index") +
  theme_void() +
  theme(legend.position = "bottom")

cowplot::plot_grid(tsf_map, sev_map, bait_map,rain_map,twi_map,
                   labels = c("a)", "b)", "c)", "d)", "e)"), ncol=3)

### Make coefficient figure
#### Step 5: Extract Model Coefs ####
library(tidybayes)
library(ggplot2)
library(AHMbook)
library(cowplot)

# Read in data
dethist = readRDS("Data_Processing/camtrapR.counthist.UpperWarren.multispp.RData")
model_output = readRDS(here::here("Data_Processing","nmix_modeloutput_baitfire_plussev.RData"))
model_output = readRDS(here::here("Data_Processing","nmix_modeloutput_16112022.RData"))
source("Scripts/3b_ModelFunctions_Final.R")

model.samples=model_output
nspec = length(model.samples$mean$mean.lambda)
spp = c("Woylie", "Chuditch", "Koomal", "Quenda", "Roo",
        "Vulpes", "Tammar", "Numbat", "Western Brush Wallaby", "Dunnart")


c.b1 = clean_coefs_90ci("beta1", "Mean Annual Rainfall", spp, model.samples)
c.b2 = clean_coefs_90ci("beta2", "TWI", spp, model.samples)
c.b3 = clean_coefs_90ci("beta3", "Baiting", spp, model.samples)
c.b4 = clean_coefs_90ci("beta4", "Time Since Fire", spp, model.samples)
c.b5 = clean_coefs_90ci("beta5", "Time Since Fire^2", spp, model.samples)
c.b6 = clean_coefs_90ci("beta6", "Time Since Fire * Baiting", spp, model.samples)
c.b7 = clean_coefs_90ci("beta7", "Fire Severity", spp, model.samples)


coef.out = rbind(c.b1, c.b2, c.b3, c.b4,c.b5, c.b6, c.b7)
spplabels= data.frame(Species = spp,
                      SpeciesLab = c("Woylie", "Chuditch", "Koomal", "Quenda", "Kangaroo",
                                     "Red Fox", "Tammar Wallaby", "Numbat", "Western Brush Wallaby", "Dunnart"))
coef.out = left_join(coef.out, spplabels)
coef.out$SpeciesLab = as.factor(coef.out$SpeciesLab)
coef.out$SpeciesLab = factor(coef.out$SpeciesLab, levels = rev(levels(coef.out$SpeciesLab)))

coef.plot = ggplot(coef.out) + 
  geom_pointrange(aes(x = SpeciesLab, y = Mean, ymin = Lower, ymax = Upper),
                  shape = ifelse(coef.out$Lower > 0 | coef.out$Upper < 0, 21, 21),
                  fill = ifelse(coef.out$Lower > 0 | coef.out$Upper < 0, "black", "white")) + 
  geom_hline(yintercept = 0, linetype='dashed') + ylab("Coeffcient and 90% CIs") + xlab("") +
  coord_flip() + facet_wrap(~Covariate, ncol=4, scales='free_x') + theme_cowplot()
coef.plot

#### Step 6: Species Co-occurrence ####
sppA.lookup = data.frame(No = 1:nspec, SpeciesA = spp)
sppB.lookup = data.frame(No = 1:nspec, SpeciesB = spp)
sppA.lookup = left_join(sppA.lookup, spplabels, by=c("SpeciesA"="Species"))
sppB.lookup = left_join(sppB.lookup, spplabels, by=c("SpeciesB"="Species"))

rho = 
  model.samples %>%
  tidy_draws() %>%
  gather_variables() %>%
  filter(grepl("rho",.variable)) %>%
  group_by(.variable) %>% 
  summarise(Mean = mean(.value), Lower_2.5 = quantile(.value, probs = 0.025), Upper_97.5 = quantile(.value, probs = 0.975)) %>%
  mutate(SppA = as.numeric(parse_number(sub("\\,.*", "", .variable))),
         SppB = as.numeric(parse_number(sub(".*,", "", .variable)))) %>%
  left_join(sppA.lookup, by = c("SppA"="No")) %>% left_join(sppB.lookup, by = c("SppB"="No"))

spp.corplot = rho %>% filter(Mean <1) %>% ggplot( ) + 
  geom_tile(aes(x=SpeciesLab.x, y=SpeciesLab.y, fill=Mean)) + 
  scale_fill_gradient2(midpoint = 0, low = "red", mid = "white",
                       high = "blue", space = "Lab", name = "Correlation") + 
  theme_cowplot() +xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
spp.corplot

#### Step 6: Make some pretty maps of species distributions ####
library(raster)
load(here::here("Data_Clean", "covariates_stack.RData"))
covariates_masked = raster::stack(covariates_masked)
cov.stack = subset(covariates_masked, c("mean_annual_rainfall", "topographic_wetness_index", "tsf_Oct2017",
                                        "bait_intensity_2017.10.23_38month_lag", "propsev_6yrs_Oct2017_500m"))

cov.stack$tsf_Oct2017[cov.stack$tsf_Oct2017 >33] <- 33 # Truncate
#cov.stack$bait_intensity_2017.10.23_38month_lag[cov.stack$bait_intensity_2017.10.23_38month_lag > 168] <- NA # Truncate

covras.df = data.frame(rasterToPoints(cov.stack))
#covras.df$easting = standardize2match(covras.df$easting, dethist$Sites$X)
#covras.df$northing = standardize2match(covras.df$northing, dethist$Sites$Y)
covras.df$topographic_wetness_index = standardize2match(covras.df$topographic_wetness_index, dethist$Sites$twi)
covras.df$mean_annual_rainfall = standardize2match(covras.df$mean_annual_rainfall, dethist$Sites$rainfall)
covras.df$tsf_Oct2017 = standardize2match(covras.df$tsf_Oct2017, dethist$Sites$tsf.point)
covras.df$propsev_6yrs_Oct2017_500m = standardize2match(covras.df$propsev_6yrs_Oct2017_500m, dethist$Sites$PropSevere500)
covras.df$bait_intensity_2017.10.23_38month_lag = standardize2match(covras.df$bait_intensity_2017.10.23_38month_lag, dethist$Sites$Mean_Intensity_400)
#covras.df$dist_to_ag = standardize2match(covras.df$dist_to_ag, dethist$Sites$dist_to_ag)
#covras.df$dist_to_majorhydro = standardize2match(covras.df$dist_to_majorhydro, dethist$Sites$dist_to_majorhydro)
#covras.df$prop_filtered_roads_3km = standardize2match(covras.df$prop_filtered_roads_3km, dethist$Sites$prop_filtered_roads_3km)

predicted.abundances = data.frame(x=covras.df$x, y=covras.df$y)
for (s in 1:length(spp)){
  species = spp[s]
  prediction = data.frame(exp(model_output$mean$beta0[s] + 
                                model_output$mean$beta1[s] * covras.df$mean_annual_rainfall + 
                                model_output$mean$beta2[s] * covras.df$topographic_wetness_index +
                                model_output$mean$beta3[s] * covras.df$tsf_Oct2017 +
                                model_output$mean$beta4[s] * covras.df$tsf_Oct2017^2 +
                                model_output$mean$beta5[s] * covras.df$bait_intensity_2017.10.23_38month_lag +
                                model_output$mean$beta6[s] * covras.df$tsf_Oct2017 * covras.df$bait_intensity_2017.10.23_38month_lag +
                                model_output$mean$beta7[s] * covras.df$propsev_6yrs_Oct2017_500m))
  names(prediction) = species
  predicted.abundances = cbind(predicted.abundances, prediction)
}


woylie.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,Woylie)), res=res(cov.stack), crs=crs(cov.stack))
chuditch.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,Chuditch)), res=res(cov.stack), crs=crs(cov.stack))
koomal.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,Koomal)), res=res(cov.stack), crs=crs(cov.stack))
quenda.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,Quenda)), res=res(cov.stack), crs=crs(cov.stack))
roo.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,Roo)), res=res(cov.stack), crs=crs(cov.stack))
vulpes.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,Vulpes)), res=res(cov.stack), crs=crs(cov.stack))
tammar.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,Tammar)), res=res(cov.stack), crs=crs(cov.stack))
numbat.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,Numbat)), res=res(cov.stack), crs=crs(cov.stack))
wallaby.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,`Western Brush Wallaby`)), res=res(cov.stack), crs=crs(cov.stack))
dunnart.sdm = rasterFromXYZ(dplyr::select(predicted.abundances, c(x,y,Dunnart)), res=res(cov.stack), crs=crs(cov.stack))

abundance.stack = stack(woylie.sdm, chuditch.sdm, koomal.sdm, quenda.sdm, roo.sdm, vulpes.sdm, tammar.sdm, numbat.sdm, wallaby.sdm, dunnart.sdm)

library(tmap)

tm_shape(abundance.stack) +
  tm_raster(title="") +
  tm_facets(free.scales = TRUE) +
  tm_shape(sites) +
  tm_dots() + 
  tm_layout(panel.labels = names(abundance.stack))



