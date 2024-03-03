##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
## 4: Make Abundance Predicitons from Multispp NMixture Model ##
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##

library(tidyverse)
library(tidybayes)
library(AHMbook)
library(parallel)
library(rjags)
library(nimble)
library(cowplot)
source("Scripts/model_helper_functions.R")

#### Step 1: Load data ####
# Read in data
dethist = readRDS("Data_Processing/camtrapR.counthist.UpperWarren.multispp09122023.RData")
model_output = readRDS("Data_Clean/nmix_nimblemodel_final.RDS")
model.summary = MCMCvis::MCMCsummary(model_output)
model.summary$.variable = rownames(model.summary)
model.eta.lam = model.summary %>% filter(grepl("eta.lam", .variable)) 
model.eta.lam = model.eta.lam%>%
  mutate(species = as.numeric(str_match(model.eta.lam$.variable, "eta.lam\\[\\s*(\\d+),\\s*(\\d+)\\]")[, 3]))

model.betas = model_output %>% tidy_draws() %>% gather_variables() %>% filter(grepl("beta", .variable))
model.betas = separate(model.betas, .variable, into=c("beta", "species"))

species.predict = c("Chuditch", "Koomal", "Quenda", "Woylie", 
                    "Vulpes", "Numbat")

nspec=length(species.predict)
spp = species.predict

#### Step 2: Predict Relationships with Baiting and Time Since Fire ####
spp.lookup = data.frame(No = as.character(1:length(species.predict)), Species = species.predict)

### Pred Alone
model.coefs = model_output %>%
  tidy_draws() %>%
  gather_variables() %>%
  filter(grepl("beta",.variable)) %>%
  group_by(.variable) %>% 
  summarise(Mean = mean(.value), 
            Lower_10 = quantile(.value, probs = 0.1),
            Upper_90 = quantile(.value, probs = 0.9),
            Lower_5 = quantile(.value, probs = 0.05),
            Upper_95 = quantile(.value, probs = 0.95),
            Lower_2.5 = quantile(.value, probs = 0.025), 
            Upper_97.5 = quantile(.value, probs = 0.975)) %>%
  tidyr::separate(.variable, into=c("coef", "Spp")) %>%
  left_join(spp.lookup, by = c("Spp"="No")) 

bait.pred = seq(from = min(dethist$Sites$Mean_Intensity_400), to = max(dethist$Sites$Mean_Intensity_400), length.out = 100)
bait.pred.std = standardize2match(bait.pred, dethist$Sites$Mean_Intensity_400)
names(bait.pred.std) <- "Bait"

spp.lookup = data.frame(No = as.character(1:length(species.predict)), Species = species.predict)

pred.out = list()
for(i in 1:nspec){ # Loop over each observed species
  coefs = filter(model.coefs, Spp == as.character(i))
  pred.df = data.frame(
    Species = spp[i],
    Bait = bait.pred, 
    Mean = exp(filter(coefs, coef=="beta0")$Mean + 
                 filter(coefs, coef=="beta3")$Mean * bait.pred.std), # Mean
    LCI = exp(filter(coefs, coef=="beta0")$Lower_10 + 
                filter(coefs, coef=="beta3")$Lower_10 * bait.pred.std), # LCI
    UCI = exp(filter(coefs, coef=="beta0")$Upper_90 + 
                filter(coefs, coef=="beta3")$Upper_90 * bait.pred.std) # UCI
  )
  pred.out[[i]]<- pred.df
}

bait.pred.out = do.call('rbind', pred.out)

spplabels= data.frame(Species = species.predict,
                      SpeciesLab = c("Chuditch", "Koomal", "Quenda", "Woylie", 
                                     "Red Fox", "Numbat"))
bait.pred.out = left_join(bait.pred.out, spplabels)

bait.abundance.plot = ggplot(bait.pred.out) + 
  geom_line(aes(x = Bait, y = Mean)) + 
  geom_ribbon(aes(x = Bait, ymin=LCI, ymax=UCI), alpha=0.4) +
  ylab("Predicted Abundance") + xlab("Baits / km^2") + theme_cowplot() + 
  facet_wrap(~SpeciesLab, scales="free_y", nrow=2)
bait.abundance.plot

### TSF Alone
tsf.pred = seq(from = min(dethist$Sites$tsf.point), to = max(dethist$Sites$tsf.point), length.out = 100)
tsf.pred.std = standardize2match(tsf.pred, dethist$Sites$tsf.point)
spp.lookup = data.frame(No = as.character(1:length(species.predict)), Species = species.predict)
names(tsf.pred.std) <- "TSF"

pred.out = list()
for(i in 1:nspec){ # Loop over each observed species
  coefs = filter(model.coefs, Spp == as.character(i))
  pred.df = data.frame(
    Species = spp[i],
    TSF = tsf.pred,
    Mean = exp(filter(coefs, coef=="beta0")$Mean + 
                 filter(coefs, coef=="beta4")$Mean * tsf.pred.std + 
                 filter(coefs, coef=="beta5")$Mean * tsf.pred.std^2), # Mean
    LCI = exp(filter(coefs, coef=="beta0")$Lower_10 + 
                filter(coefs, coef=="beta4")$Lower_10 * tsf.pred.std + 
                filter(coefs, coef=="beta5")$Lower_10 * tsf.pred.std^2), # LCI
    UCI = exp(filter(coefs, coef=="beta0")$Upper_90 + 
                filter(coefs, coef=="beta4")$Upper_90 * tsf.pred.std + 
                filter(coefs, coef=="beta5")$Upper_90 * tsf.pred.std^2) # UCI
  )
  pred.out[[i]]<- pred.df
}

tsf.pred.out = do.call('rbind', pred.out)

tsf.pred.out= left_join(tsf.pred.out, spplabels)

tsf.abundance.plot = ggplot(tsf.pred.out) + 
  geom_line(aes(x = TSF, y = Mean)) + 
  geom_ribbon(aes(x = TSF, ymin=LCI, ymax=UCI), alpha=0.4) +
  ylab("Predicted Abundance") + xlab("Time Since Fire") + theme_cowplot() + 
  facet_wrap(~SpeciesLab, scales="free_y", nrow=2)
tsf.abundance.plot

### Pred time since fire
tsf.pred = seq(from = min(dethist$Sites$tsf.point), to = max(dethist$Sites$tsf.point), length.out = 100)
tsf.pred.std = standardize2match(tsf.pred, dethist$Sites$tsf.point)
bait.pred = seq(from = min(dethist$Sites$Mean_Intensity_400), to = max(dethist$Sites$Mean_Intensity_400), length.out = 100)
bait.pred.std = standardize2match(bait.pred, dethist$Sites$Mean_Intensity_400)
pred.grid = expand.grid(tsf.pred, bait.pred); names(pred.grid) <-c("TSF", "Bait")
pred.grid.std = expand.grid(tsf.pred.std, bait.pred.std); names(pred.grid.std) <-c("TSF", "Bait")

spp.lookup = data.frame(No = as.character(1:length(species.predict)), Species = species.predict)

pred.out = list()
for(i in 1:nspec){ # Loop over each observed species
  coefs = filter(model.coefs, Spp == as.character(i))
  pred.df = data.frame(
    Species = spp[i],
    TSF = pred.grid$TSF,
    Bait = pred.grid$Bait, 
    Mean = exp(filter(coefs, coef=="beta0")$Mean + 
                 filter(coefs, coef=="beta4")$Mean * pred.grid.std$TSF + 
                 filter(coefs, coef=="beta5")$Mean * pred.grid.std$TSF^2 + 
                 filter(coefs, coef=="beta3")$Mean * pred.grid.std$Bait + 
                 filter(coefs, coef=="beta6")$Mean * pred.grid.std$Bait * pred.grid.std$TSF), # Mean
    LCI = exp(filter(coefs, coef=="beta0")$Lower_10 + 
                filter(coefs, coef=="beta4")$Lower_10 * pred.grid.std$TSF + 
                filter(coefs, coef=="beta5")$Lower_10 * pred.grid.std$TSF^2 + 
                filter(coefs, coef=="beta3")$Lower_10 * pred.grid.std$Bait + 
                filter(coefs, coef=="beta6")$Lower_10 * pred.grid.std$Bait * pred.grid.std$TSF), # LCI
    UCI = exp(filter(coefs, coef=="beta0")$Upper_90 + 
                filter(coefs, coef=="beta4")$Upper_90 * pred.grid.std$TSF + 
                filter(coefs, coef=="beta5")$Upper_90 * pred.grid.std$TSF^2 + 
                filter(coefs, coef=="beta3")$Upper_90 * pred.grid.std$Bait + 
                filter(coefs, coef=="beta6")$Upper_90 * pred.grid.std$Bait * pred.grid.std$TSF) # UCI
  )
  pred.out[[i]]<- pred.df
}

pred.out = do.call('rbind', pred.out)
pred.out$BaitF = as.factor(round(pred.out$Bait, 2))
pred.out.sub = filter(pred.out, pred.out$BaitF %in% c("4.91","80.07", "120.16"))

pred.out.sub = left_join(pred.out.sub, spplabels)

abundance.plot = ggplot(pred.out.sub) + 
  geom_line(aes(x = TSF, y = Mean, colour = BaitF)) + 
  geom_ribbon(aes(x = TSF, ymin=LCI, ymax=UCI, fill = BaitF), alpha=0.4) +
  ylab("Predicted Abundance") + xlab("Time Since Fire") + labs(colour = "Baits / km^2", fill= "Baits / km^2") +
  scale_fill_viridis_d(labels=c("4","80","120")) + scale_colour_viridis_d(labels=c("4","80","120")) + theme_cowplot() + 
  theme(legend.position="bottom") +
  facet_wrap(~SpeciesLab, scales="free_y", nrow=2)

abundance.plot

abundance.plots = plot_grid(bait.abundance.plot, tsf.abundance.plot, abundance.plot, 
                            labels = c("a)", "b)", "c)"), nrow=3)


ggsave(plot = abundance.plots, filename="predictedabundance_plot_final.pdf", path="Outputs", device="pdf", width=4,height=6,units="in",scale=3)

ggsave(plot = abundance.plot, filename="predictedabundance_interaction_plot_final.pdf", path="Outputs", device="pdf", width=4,height=2,units="in",scale=3)


#### Step 4: Create management scenarios to predict species abundances into ####
## High, aerial_groundmonthlyium and low baiting

## Base Scenario (Current Situation)
consistent.covs = dethist$Sites %>% 
  dplyr::select(LocationName, Site, "east" = X, "north" = Y, rainfall, twi) 

# Add in variables for baseline scenario 
base = consistent.covs %>% mutate(bait = dethist$Sites$Mean_Intensity_400)
base[,3:7] = sapply(base[,3:7], FUN = function(x) {as.numeric(standardize(x))}) # Standardise covariates

## Expand to tsf scenarios
make_tsf_scenarios = function(scen_base, bait_val, sev_val, scenario_name){
  scenarios = scen_base
  scenarios = scenarios %>% mutate(bait = standardize2match(bait_val, dethist$Sites$Mean_Intensity_400))
  tsf.scenarios = list()
  for (t in 1:6){
    scenarios = scenarios %>% mutate(tsf = standardize2match(t, dethist$Sites$tsf.point),
                                     tsf_actual = t,
                                     propsev = standardize2match(sev_val, dethist$Sites$PropSevere500))
    tsf.scenarios[[paste0(scenario_name, "_tsf",t)]] <- scenarios
  }
  for (t in 7:33){
    scenarios = scenarios %>% mutate(tsf = standardize2match(t, dethist$Sites$tsf.point),
                                     tsf_actual = t,
                                     propsev = standardize2match(0, dethist$Sites$PropSevere500))
    tsf.scenarios[[paste0(scenario_name, "_tsf",t)]] <- scenarios
  }
  tsf.scenarios = do.call("rbind", tsf.scenarios)
  return(tsf.scenarios)
}

s_b0_sL = make_tsf_scenarios(scen_base =base, bait_val = dethist$Sites$Mean_Intensity_400, sev_val = 0.2, scenario_name = "s_b0_sL")
s_b0_sM = make_tsf_scenarios(scen_base =base, bait_val = dethist$Sites$Mean_Intensity_400, sev_val = 0.5, scenario_name = "s_b0_sM")
s_b0_sH = make_tsf_scenarios(scen_base =base, bait_val = dethist$Sites$Mean_Intensity_400, sev_val = 0.8, scenario_name = "s_b0_sH")

s_bL_sL = make_tsf_scenarios(scen_base =base, bait_val = 5, sev_val = 0.2, scenario_name = "s_bL_sL")
s_bL_sM = make_tsf_scenarios(scen_base =base, bait_val = 5, sev_val = 0.5, scenario_name = "s_bL_sM")
s_bL_sH = make_tsf_scenarios(scen_base =base, bait_val = 5, sev_val = 0.8, scenario_name = "s_bL_sH")

s_bA_sL = make_tsf_scenarios(scen_base =base, bait_val = 59, sev_val = 0.2, scenario_name = "s_bA_sL")
s_bA_sM = make_tsf_scenarios(scen_base =base, bait_val = 59, sev_val = 0.5, scenario_name = "s_bA_sM")
s_bA_sH = make_tsf_scenarios(scen_base =base, bait_val = 59, sev_val = 0.8, scenario_name = "s_bA_sH")

s_bAG_sL = make_tsf_scenarios(scen_base =base, bait_val = 150, sev_val = 0.2, scenario_name = "s_bA_sL")
s_bAG_sM = make_tsf_scenarios(scen_base =base, bait_val = 150, sev_val = 0.5, scenario_name = "s_bA_sM")
s_bAG_sH = make_tsf_scenarios(scen_base =base, bait_val = 150, sev_val = 0.8, scenario_name = "s_bA_sH")

#### Step 3: Predictions of species abundances, sampling from posterior ####
species.predict = c("Chuditch", "Koomal", "Quenda", "Woylie", 
                    "Red Fox", "Numbat")

nsamples = 1000 # how many samples from the posterior

start = Sys.time()
s_b0_sL_abundance = predict.abundance.posterior(s_b0_sL, species.list = species.predict, nsamp=nsamples)
saveRDS(s_b0_sL_abundance, "Data_Processing/UpperWarrenAbundance_s_b0_sL_abundance_scenarios.RData"); rm(s_b0_sL_abundance) 
s_b0_sM_abundance = predict.abundance.posterior(s_b0_sM, species.list = species.predict, nsamp=nsamples)
saveRDS(s_b0_sM_abundance, "Data_Processing/UpperWarrenAbundance_s_b0_sM_abundance_scenarios.RData"); rm(s_b0_sM_abundance) 
s_b0_sH_abundance = predict.abundance.posterior(s_b0_sH, species.list = species.predict, nsamp=nsamples)
saveRDS(s_b0_sH_abundance, "Data_Processing/UpperWarrenAbundance_s_b0_sH_abundance_scenarios.RData"); rm(s_b0_sH_abundance) 
timetaken = Sys.time() - start
timetaken
s_bL_sL_abundance = predict.abundance.posterior(s_bL_sL, species.list = species.predict, nsamp=nsamples)
saveRDS(s_bL_sL_abundance, "Data_Processing/UpperWarrenAbundance_s_bL_sL_abundance_scenarios.RData"); rm(s_bL_sL_abundance) 
s_bL_sM_abundance = predict.abundance.posterior(s_bL_sM, species.list = species.predict, nsamp=nsamples)
saveRDS(s_bL_sM_abundance, "Data_Processing/UpperWarrenAbundance_s_bL_sM_abundance_scenarios.RData"); rm(s_bL_sM_abundance) 
s_bL_sH_abundance = predict.abundance.posterior(s_bL_sH, species.list = species.predict, nsamp=nsamples)
saveRDS(s_bL_sH_abundance, "Data_Processing/UpperWarrenAbundance_s_bL_sH_abundance_scenarios.RData"); rm(s_bL_sH_abundance) 
timetaken = Sys.time() - start
timetaken
s_bA_sL_abundance = predict.abundance.posterior(s_bA_sL, species.list = species.predict, nsamp=nsamples)
saveRDS(s_bA_sL_abundance, "Data_Processing/UpperWarrenAbundance_s_bA_sL_abundance_scenarios.RData"); rm(s_bA_sL_abundance) 
s_bA_sM_abundance = predict.abundance.posterior(s_bA_sM, species.list = species.predict, nsamp=nsamples)
saveRDS(s_bA_sM_abundance, "Data_Processing/UpperWarrenAbundance_s_bA_sM_abundance_scenarios.RData"); rm(s_bA_sM_abundance) 
s_bA_sH_abundance = predict.abundance.posterior(s_bA_sH, species.list = species.predict, nsamp=nsamples)
saveRDS(s_bA_sH_abundance, "Data_Processing/UpperWarrenAbundance_s_bA_sH_abundance_scenarios.RData"); rm(s_bA_sH_abundance) 
timetaken = Sys.time() - start
timetaken
s_bAG_sL_abundance = predict.abundance.posterior(s_bAG_sL, species.list = species.predict, nsamp=nsamples)
saveRDS(s_bAG_sL_abundance, "Data_Processing/UpperWarrenAbundance_s_bAG_sL_abundance_scenarios.RData"); rm(s_bAG_sL_abundance) 
s_bAG_sM_abundance = predict.abundance.posterior(s_bAG_sM, species.list = species.predict, nsamp=nsamples)
saveRDS(s_bAG_sM_abundance, "Data_Processing/UpperWarrenAbundance_s_bAG_sM_abundance_scenarios.RData"); rm(s_bAG_sM_abundance) 
s_bAG_sH_abundance = predict.abundance.posterior(s_bAG_sH, species.list = species.predict, nsamp=nsamples)
saveRDS(s_bAG_sH_abundance, "Data_Processing/UpperWarrenAbundance_s_bAG_sH_abundance_scenarios.RData"); rm(s_bAG_sH_abundance) 

timetaken = Sys.time() - start
timetaken
