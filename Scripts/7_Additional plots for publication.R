##+++++++++++++++++++++++++++++++++++++##
## 6: Additional Plots for Publication ##
##+++++++++++++++++++++++++++++++++++++##

### Covariate maps for Supplementary Information

library(raster)
library(ggplot2)
library(cowplot)
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
            theme(legend.position = "bottom") +
  coord_equal()

bait_map = ggplot(bait) + geom_raster(aes(x=x,y=y, fill=bait_intensity_2016.10.01_38month_lag)) + 
  scale_fill_viridis_c(name="Bait Intensity") +
  theme_void() +
  theme(legend.position = "bottom")+
  coord_equal()

sev_map = ggplot(sev) + geom_raster(aes(x=x,y=y, fill=propsev_6yrs_Oct2016_500m)) + 
  scale_fill_viridis_c(name="Proportion Severely Burnt") +
  theme_void() +
  theme(legend.position = "bottom")+
  coord_equal()

rain_map = ggplot(rain) + geom_raster(aes(x=x,y=y, fill=mean_annual_rainfall)) + 
  scale_fill_viridis_c(name="Mean Annual Rainfall") +
  theme_void() +
  theme(legend.position = "bottom")+
  coord_equal()

twi_map = ggplot(twi) + 
  geom_raster(aes(x=x,y=y, fill=topographic_wetness_index)) + 
  scale_fill_viridis_c(name="Topographic Wetness Index") +
  theme_void() +
  theme(legend.position = "bottom")+
  coord_equal()

map.plot = cowplot::plot_grid(tsf_map, sev_map, bait_map,rain_map,twi_map,
                   labels = c("a)", "b)", "c)", "d)", "e)"), ncol=3)

ggsave(plot = map.plot, filename = "Outputs/covariate_map.pdf", dpi = 300, width = 10, height=6)


### Make coefficient figure
#### Step 5: Extract Model Coefs ####
library(tidybayes)
library(ggplot2)
library(AHMbook)
library(cowplot)
library(tidyverse)
# Read in data
dethist = readRDS("Data_Processing/camtrapR.counthist.UpperWarren.multispp.RData")
model_output = readRDS("Data_Clean/nmix_nimblemodel_final.RDS")
source("Scripts/model_helper_functions.R")

model.samples=model_output
spp = c("Chuditch", "Koomal", "Quenda", "Woylie", 
        "Vulpes", "Numbat")
nspec = length(spp)

c.b1 = clean_coefs_90ci("beta1", "Mean Annual Rainfall", spp, model.samples)
c.b2 = clean_coefs_90ci("beta2", "Topographic Wetness Index", spp, model.samples)
c.b3 = clean_coefs_90ci("beta3", "Baiting", spp, model.samples)
c.b4 = clean_coefs_90ci("beta4", "Time Since Fire", spp, model.samples)
c.b5 = clean_coefs_90ci("beta5", "Time Since Fire^2", spp, model.samples)
c.b6 = clean_coefs_90ci("beta6", "Time Since Fire * Baiting", spp, model.samples)
c.b7 = clean_coefs_90ci("beta7", "Fire Severity", spp, model.samples)

coef.out = rbind(c.b1, c.b2, c.b3, c.b4,c.b5, c.b6, c.b7)
spplabels= data.frame(Species = spp,
                      SpeciesLab = c("Chuditch", "Koomal", "Quenda", "Woylie", 
                                     "Red Fox",  "Numbat"))
coef.out = left_join(coef.out, spplabels)
coef.out$SpeciesLab = as.factor(coef.out$SpeciesLab)
coef.out$SpeciesLab = factor(coef.out$SpeciesLab, levels = rev(levels(coef.out$SpeciesLab)))

coef.plot = ggplot(coef.out) + 
  geom_pointrange(aes(x = Covariate, y = Mean, ymin = Lower, ymax = Upper),
                  shape = ifelse(coef.out$Lower > 0 | coef.out$Upper < 0, 21, 21),
                  fill = ifelse(coef.out$Lower > 0 | coef.out$Upper < 0, "black", "white")) + 
  geom_hline(yintercept = 0, linetype='dashed') + ylab("Coeffcient and 90% CIs") + xlab("") +
  coord_flip() + facet_wrap(~SpeciesLab, ncol=3, scales='free_x') + theme_cowplot()
coef.plot

ggsave(plot = coef.plot, filename = "Outputs/abundance_coefficients.pdf", dpi = 300, width = 8, height=4, scale = 1.3)


#### Make Species Co-occurrence Figure ####
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

ggsave(plot = spp.corplot, filename = "Outputs/species_corrplot.pdf", dpi = 300, width = 6, height=5, scale = 1)


