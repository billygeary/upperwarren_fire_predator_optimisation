##++++++++++++++++++++++++++++++++++++++++++++++++++##
## 6: Plot the Optimisation Results for Publication ##
##++++++++++++++++++++++++++++++++++++++++++++++++++##

library(prioritizr)
library(tidyverse)
library(ggplot2)
library(cowplot)

###############################################
#### Step 1: Read in the results summaries ####
###############################################

results.paths = list.files(path = "Data_Clean", pattern = "baitfireoptimisationresults_transectscale_summaries_s", full = TRUE)

results = data.frame()
for (i in seq_along(results.paths)){
  t = read.csv(results.paths[i])
  results = rbind(results, t)
}
summaries.out = results

##################################################
#### Step 2: Check the data for all secnarios ####
##################################################
head(summaries.out)

plot = summaries.out %>% drop_na() %>% 
  group_by(.id, Scenario, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years"))) %>%
  separate(Scenario, into=c('Scen', 'Baiting', 'Severity')) %>%
  filter(Baiting != "b0") %>%
  #mutate(Scenario = factor(Scenario, levels = c("Low", "Aerial", "Aerial + Ground (Monthly)", "Aerial + Ground (Fortnightly)", "Current"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=TSFCAT, y=cost/11, fill = Severity)) + facet_wrap(~Baiting, ncol = 1) +
  theme_cowplot() + scale_fill_viridis_d() +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
plot

plot = summaries.out %>% drop_na() %>% 
  group_by(.id, Scenario, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years"))) %>%
  #separate(Scenario, into=c('Scen', 'Baiting', 'Severity')) %>%
  #filter(Baiting != "b0") %>%
  #mutate(Scenario = factor(Scenario, levels = c("Low", "Aerial", "Aerial + Ground (Monthly)", "Aerial + Ground (Fortnightly)", "Current"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=TSFCAT, y=cost/11, fill = TSFCAT)) + facet_wrap(~Scenario, ncol = 3) +
  theme_cowplot() + scale_fill_viridis_d() +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
plot

############################################
#### Step 3: Optimisation Results Plots ####
############################################

#### Step 3a: High severity plots only
results.paths = intersect(list.files(path = "Data_Clean", pattern = "baitfireoptimisationresults_transectscale_summaries_s", full = TRUE),
                          list.files(path = "Data_Clean", pattern = "_sH", full = TRUE))
results.paths = results.paths[2:4]
results.paths

highsev.summaries.out = data.frame()
for (i in seq_along(results.paths)){
  t = read.csv(results.paths[i])
  highsev.summaries.out = rbind(highsev.summaries.out, t)
}

plot = highsev.summaries.out %>% drop_na() %>% 
  group_by(Scenario,.id, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=TSFCAT, y=cost/11, fill = Scenario)) + #facet_wrap(~Baiting, ncol = 1) +
  theme_cowplot() + scale_fill_viridis_d() +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
plot

plot.quantiles = highsev.summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years")),
         Scenario = as.factor(Scenario)) %>%
  group_by(Scenario,TSFCAT) %>% summarise(Median = quantile(cost, probs=0.5),
                                          Mean = mean(cost),
                                          LCI = quantile(cost, probs=0.1),
                                          UCI = quantile(cost, probs=0.90))

facet.labs =  c(s_bL_sH="Low", s_bA_sH="Medium", s_bAG_sH="High")

highsev.plot = highsev.summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years")),
         Scenario = factor(Scenario, levels = c("s_bL_sH", "s_bA_sH", "s_bAG_sH"))) %>% 
  ggplot() + 
  geom_jitter(aes(x=TSFCAT,y=cost/11, colour= TSFCAT), alpha=0.4) + scale_color_viridis_d() + 
  geom_pointrange(data = plot.quantiles, 
                  aes(x=TSFCAT, 
                      y=Mean/11, 
                      ymin=LCI/11, 
                      ymax=UCI/11)) +
  theme_cowplot() + facet_wrap(~Scenario, nrow=3, labeller=labeller(Scenario=facet.labs)) +
  theme(legend.position="none") +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
highsev.plot  

ggsave(plot = highsev.plot, filename="outcome_plot_highsev_final.pdf", path="Outputs", device="pdf", width=4,height=6,units="in",scale=1.5)

#### Step 3b: Low Severity Plots ####
results.paths = intersect(list.files(path = "Data_Clean", pattern = "baitfireoptimisationresults_transectscale_summaries_s", full = TRUE),
                          list.files(path = "Data_Clean", pattern = "_sL", full = TRUE))
results.paths = results.paths[2:4]
results.paths

lowsev.summaries.out = data.frame()
for (i in seq_along(results.paths)){
  t = read.csv(results.paths[i])
  lowsev.summaries.out = rbind(lowsev.summaries.out, t)
}

plot = lowsev.summaries.out %>% drop_na() %>% 
  group_by(Scenario,.id, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=TSFCAT, y=cost/11, fill = Scenario)) + #facet_wrap(~Baiting, ncol = 1) +
  theme_cowplot() + scale_fill_viridis_d() +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
plot

plot.quantiles = lowsev.summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years")),
         Scenario = as.factor(Scenario)) %>%
  group_by(Scenario,TSFCAT) %>% summarise(Median = quantile(cost, probs=0.5),
                                          Mean = mean(cost),
                                          LCI = quantile(cost, probs=0.1),
                                          UCI = quantile(cost, probs=0.90))

facet.labs =  c(s_bL_sL="Low", s_bA_sL="Medium", s_bAG_sL="High")

lowsev.plot = lowsev.summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years")),
         Scenario = factor(Scenario, levels = c("s_bL_sL", "s_bA_sL", "s_bAG_sL"))) %>% 
  ggplot() + 
  geom_jitter(aes(x=TSFCAT,y=cost/11, colour= TSFCAT), alpha=0.4) + scale_color_viridis_d() + 
  geom_pointrange(data = plot.quantiles, 
                  aes(x=TSFCAT, 
                      y=Mean/11, 
                      ymin=LCI/11, 
                      ymax=UCI/11)) +
  theme_cowplot() + facet_wrap(~Scenario, nrow=3, labeller=labeller(Scenario=facet.labs)) +
  theme(legend.position="none") +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
lowsev.plot  

ggsave(plot = lowsev.plot, filename="outcome_plot_lowsev_final.pdf", path="Outputs", device="pdf", width=4,height=6,units="in",scale=1.5)

#### Step 3c: Medium Severity Plots ####
results.paths = intersect(list.files(path = "Data_Clean", pattern = "baitfireoptimisationresults_transectscale_summaries_s", full = TRUE),
                          list.files(path = "Data_Clean", pattern = "_sM", full = TRUE))
results.paths = results.paths[2:4]
results.paths

medsev.summaries.out = data.frame()
for (i in seq_along(results.paths)){
  t = read.csv(results.paths[i])
  medsev.summaries.out = rbind(medsev.summaries.out, t)
}

plot = medsev.summaries.out %>% drop_na() %>% 
  group_by(Scenario,.id, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=TSFCAT, y=cost/11, fill = Scenario)) + #facet_wrap(~Baiting, ncol = 1) +
  theme_cowplot() + scale_fill_viridis_d() +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
plot

plot.quantiles = medsev.summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years")),
         Scenario = as.factor(Scenario)) %>%
  group_by(Scenario,TSFCAT) %>% summarise(Median = quantile(cost, probs=0.5),
                                          Mean = mean(cost),
                                          LCI = quantile(cost, probs=0.1),
                                          UCI = quantile(cost, probs=0.90))

facet.labs =  c(s_bL_sM="Low", s_bA_sM="Medium", s_bAG_sM="High")

medsev.plot = medsev.summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id, TSFCAT) %>% summarise(cost = sum(n)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years")),
         Scenario = factor(Scenario, levels = c("s_bL_sM", "s_bA_sM", "s_bAG_sM"))) %>% 
  ggplot() + 
  geom_jitter(aes(x=TSFCAT,y=cost/11, colour= TSFCAT), alpha=0.4) + scale_color_viridis_d() + 
  geom_pointrange(data = plot.quantiles, 
                  aes(x=TSFCAT, 
                      y=Mean/11, 
                      ymin=LCI/11, 
                      ymax=UCI/11)) +
  theme_cowplot() + facet_wrap(~Scenario, nrow=3, labeller=labeller(Scenario=facet.labs)) +
  theme(legend.position="none") +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
medsev.plot  

ggsave(plot = plot, filename="outcome_plot_medsev_final.pdf", path="Outputs", device="pdf", width=4,height=6,units="in",scale=1.5)


allplot = cowplot::plot_grid(lowsev.plot, medsev.plot, highsev.plot, labels = c("a)", "b)", "c)"), nrow=1)

ggsave(plot = allplot, filename="outcome_plot_allsev_final.pdf", path="Outputs", device="pdf", width=12,height=6,units="in",scale=1.5)

##############################################
#### Step 4: Optimisation Abundance Plots ####
##############################################

#### Step 4a: High severity abundance plots only 
results.paths = intersect(list.files(path = "Data_Clean", pattern = "baitfireoptimisationresults_transectscale_abundances_s", full = TRUE),
                          list.files(path = "Data_Clean", pattern = "_sH", full = TRUE))
results.paths = results.paths[2:4]
results.paths

highsev.abundances.out = data.frame()
for (i in seq_along(results.paths)){
  t = read.csv(results.paths[i])
  highsev.abundances.out = rbind(highsev.abundances.out, t)
}

facet.labs =  data.frame(Scenario = c("s_bL_sH", "s_bA_sH","s_bAG_sH"),
                         ScenarioLab = c("Low", "Medium", "High"))

abundance.plot.quantiles = highsev.abundances.out %>% 
  filter(summary=="overall") %>%
  left_join(facet.labs, by="Scenario") %>%
  group_by(ScenarioLab, Species) %>% summarise(Median = quantile(absolute_held, probs=0.5),
                                               Mean = mean(absolute_held),
                                               LCI = quantile(absolute_held, probs=0.1),
                                               UCI = quantile(absolute_held, probs=0.9))

overall.plot = highsev.abundances.out %>% filter(summary=="overall") %>%
  left_join(facet.labs, by="Scenario") %>%
  mutate(ScenarioLab = factor(ScenarioLab, levels = c("Low", "Medium", "High"))) %>%
  ggplot() +
  geom_jitter(aes(x=ScenarioLab, y = absolute_held, color=ScenarioLab),alpha=0.4) + 
  #ylim(c(0,80000)) +
  geom_pointrange(data = abundance.plot.quantiles, 
                  aes(x=ScenarioLab, 
                      y=Mean, 
                      ymin=LCI, 
                      ymax=UCI)) +
  facet_wrap(~Species, nrow=2, scales="free_y") +
  theme_cowplot() + scale_color_viridis_d(begin=0.2) + 
  ylab("Abundance") + xlab("") + theme(legend.position="none")
overall.plot

ggsave(plot = overall.plot, filename="abundance_plot_highsev_final.pdf", path="Outputs", device="pdf", width=4,height=3,units="in",scale=2.2)

#### Step 4b: Medium severity abundance plots only 
results.paths = intersect(list.files(path = "Data_Clean", pattern = "baitfireoptimisationresults_transectscale_abundances_s", full = TRUE),
                          list.files(path = "Data_Clean", pattern = "_sM", full = TRUE))
results.paths = results.paths[2:4]
results.paths

medsev.abundances.out = data.frame()
for (i in seq_along(results.paths)){
  t = read.csv(results.paths[i])
  medsev.abundances.out = rbind(medsev.abundances.out, t)
}

facet.labs =  data.frame(Scenario = c("s_bL_sM", "s_bA_sM","s_bAG_sM"),
                         ScenarioLab = c("Low", "Medium", "High"))

abundance.plot.quantiles = medsev.abundances.out %>% 
  filter(summary=="overall") %>%
  left_join(facet.labs, by="Scenario") %>%
  group_by(ScenarioLab, Species) %>% summarise(Median = quantile(absolute_held, probs=0.5),
                                               Mean = mean(absolute_held),
                                               LCI = quantile(absolute_held, probs=0.1),
                                               UCI = quantile(absolute_held, probs=0.9))

overall.plot = medsev.abundances.out %>% filter(summary=="overall") %>%
  left_join(facet.labs, by="Scenario") %>%
  mutate(ScenarioLab = factor(ScenarioLab, levels = c("Low", "Medium", "High"))) %>%
  ggplot() +
  geom_jitter(aes(x=ScenarioLab, y = absolute_held, color=ScenarioLab),alpha=0.4) + 
  #ylim(c(0,100000)) +
  geom_pointrange(data = abundance.plot.quantiles, 
                  aes(x=ScenarioLab, 
                      y=Mean, 
                      ymin=LCI, 
                      ymax=UCI)) +
  facet_wrap(~Species, nrow=2, scales="free_y") +
  theme_cowplot() + scale_color_viridis_d(begin=0.2) + 
  ylab("Abundance") + xlab("") + theme(legend.position="none")
overall.plot

ggsave(plot = overall.plot, filename="abundance_plot_medsev_final.pdf", path="Outputs", device="pdf", width=4,height=3,units="in",scale=2.2)

#### Step 4c: Low severity abundance plots only 
results.paths = intersect(list.files(path = "Data_Clean", pattern = "baitfireoptimisationresults_transectscale_abundances_s", full = TRUE),
                          list.files(path = "Data_Clean", pattern = "_sL", full = TRUE))
results.paths = results.paths[2:4]
results.paths

lowsev.abundances.out = data.frame()
for (i in seq_along(results.paths)){
  t = read.csv(results.paths[i])
  lowsev.abundances.out = rbind(lowsev.abundances.out, t)
}

facet.labs =  data.frame(Scenario = c("s_bL_sL", "s_bA_sL","s_bAG_sL"),
                         ScenarioLab = c("Low", "Medium", "High"))

abundance.plot.quantiles = lowsev.abundances.out %>% 
  filter(summary=="overall") %>%
  left_join(facet.labs, by="Scenario") %>%
  group_by(ScenarioLab, Species) %>% summarise(Median = quantile(absolute_held, probs=0.5),
                                               Mean = mean(absolute_held),
                                               LCI = quantile(absolute_held, probs=0.1),
                                               UCI = quantile(absolute_held, probs=0.9))

overall.plot = lowsev.abundances.out %>% filter(summary=="overall") %>%
  left_join(facet.labs, by="Scenario") %>%
  mutate(ScenarioLab = factor(ScenarioLab, levels = c("Low", "Medium", "High"))) %>%
  ggplot() +
  geom_jitter(aes(x=ScenarioLab, y = absolute_held, color=ScenarioLab),alpha=0.4) + 
  #ylim(c(0,100000)) +
  geom_pointrange(data = abundance.plot.quantiles, 
                  aes(x=ScenarioLab, 
                      y=Mean, 
                      ymin=LCI, 
                      ymax=UCI)) +
  facet_wrap(~Species, nrow=2, scales="free_y") +
  theme_cowplot() + scale_color_viridis_d(begin=0.2) + 
  ylab("Abundance") + xlab("") + theme(legend.position="none")
overall.plot

ggsave(plot = overall.plot, filename="abundance_plot_lowsev_final.pdf", path="Outputs", device="pdf", width=4,height=3,units="in",scale=2.2)

#############################################
#### Step 5: Optimisation Transect Plots ####
#############################################

#### Step 5a: Medium Severity Transect Plots #### 
transect.data = intersect(list.files(path = "Data_Clean", pattern = "baitfireoptimisationresults_transectscale_outcomes_", full = TRUE),
                          list.files(path = "Data_Clean", pattern = "_sM", full = TRUE))
transect.data  = transect.data[2:4]
transect.data 

medsev.transects.out = data.frame()
for (i in seq_along(transect.data)){
  t = read.csv(transect.data[i])
  medsev.transects.out = rbind(medsev.transects.out, t)
}


facet.labs =  data.frame(Scenario = c("s_bL_sM", "s_bA_sM","s_bAG_sM"),
                         ScenarioLab = c("Low", "Medium", "High"))

# First we need to summarise to TSFCAT level by summing props across tsf vals, per iteration
transects.summary = medsev.transects.out %>% group_by(.id, Site, TSFCAT, Scenario) %>% summarise(prop_landscape = sum(prop_landscape)) 

# Then we want the mean prop landscape across iterations - this gives the mean proportion of the landscape chosen and CIs
plot.data = transects.summary %>% group_by(Site, TSFCAT, Scenario) %>% summarise(mean_prop_landscape=mean(prop_landscape),
                                                                                         LCI = quantile(prop_landscape, probs = 0.1),
                                                                                         UCI = quantile(prop_landscape, probs = 0.9)) %>%
  # Then add some plotting labels to make the plots nice
  left_join(facet.labs, by="Scenario") %>%
  mutate(ScenarioLab = factor(ScenarioLab, levels = c("Low", "Medium", "High"))) %>%
  mutate(TSFCATlab = gsub(" years", "", TSFCAT)) %>%
  mutate(TSFCATlab = factor(TSFCATlab, levels = c("0-5", "6-10", "11-15", "16-24", ">24")))

big.plot = plot.data %>% ggplot() + 
  geom_pointrange(aes(x=TSFCATlab, y=mean_prop_landscape, ymin =LCI, ymax = UCI, colour = ScenarioLab), position = position_dodge(width = 0.5)) +
  theme_cowplot() + scale_color_viridis_d() + 
  ylab("Proportion of Landscape") + xlab("Time Since Fire") + facet_wrap(~Site, ncol = 3) + labs(color = "Bait Intensity")
big.plot

ggsave(plot = big.plot, filename="transect_outputs_medsev_plot_final.pdf", path="Outputs", device="pdf", width=6.5,height=5,units="in",scale=1.55)

### Per transect
for (t in 1:length(unique(plot.data$Site))){
  transect = unique(plot.data$Site)[t]
  dat = plot.data %>% filter(Site == transect) %>% filter(ScenarioLab=="High")
  
  plot.out = ggplot(dat) + 
    geom_col(aes(x=TSFCATlab, y=mean_prop_landscape, fill=TSFCATlab)) + ylim(c(0,1)) + 
    theme_cowplot() + scale_fill_viridis_d() + theme(legend.position="none") + 
    ylab("Proportion of Landscape") + xlab("Time Since Fire")
  plot.out
  
  ggsave(plot.out,device="png", 
         path=here::here("Outputs"),
         width = 2.5,height =1.25, units = "in", dpi = 320, scale=2.25,
         filename = paste0("transectoptimisation_",transect,".png"))
}


#### Step 5b: High Severity Transect Plots #### 
transect.data = intersect(list.files(path = "Data_Clean", pattern = "baitfireoptimisationresults_transectscale_outcomes_", full = TRUE),
                          list.files(path = "Data_Clean", pattern = "_sH", full = TRUE))
transect.data  = transect.data[2:4]
transect.data 

highsev.transects.out = data.frame()
for (i in seq_along(transect.data)){
  t = read.csv(transect.data[i])
  highsev.transects.out = rbind(highsev.transects.out, t)
}

facet.labs =  data.frame(Scenario = c("s_bL_sH", "s_bA_sH","s_bAG_sH"),
                         ScenarioLab = c("Low", "Medium", "High"))

# First we need to summarise to TSFCAT level by summing props across tsf vals, per iteration
transects.summary = highsev.transects.out %>% group_by(.id, Site, TSFCAT, Scenario) %>% summarise(prop_landscape = sum(prop_landscape)) 

# Then we want the mean prop landscape across iterations - this gives the mean proportion of the landscape chosen and CIs
plot.data = transects.summary %>% group_by(Site, TSFCAT, Scenario) %>% summarise(mean_prop_landscape=mean(prop_landscape),
                                                                                 LCI = quantile(prop_landscape, probs = 0.1),
                                                                                 UCI = quantile(prop_landscape, probs = 0.9)) %>%
  # Then add some plotting labels to make the plots nice
  left_join(facet.labs, by="Scenario") %>%
  mutate(ScenarioLab = factor(ScenarioLab, levels = c("Low", "Medium", "High"))) %>%
  mutate(TSFCATlab = gsub(" years", "", TSFCAT)) %>%
  mutate(TSFCATlab = factor(TSFCATlab, levels = c("0-5", "6-10", "11-15", "16-24", ">24")))

big.plot = plot.data %>% ggplot() + 
  geom_pointrange(aes(x=TSFCATlab, y=mean_prop_landscape, ymin =LCI, ymax = UCI, colour = ScenarioLab), position = position_dodge(width = 0.5)) +
  theme_cowplot() + scale_color_viridis_d() + 
  ylab("Proportion of Landscape") + xlab("Time Since Fire") + facet_wrap(~Site, ncol = 3) + labs(color = "Bait Intensity")
big.plot

ggsave(plot = big.plot, filename="transect_outputs_highsev_plot_final.pdf", path="Outputs", device="pdf", width=6.5,height=5,units="in",scale=1.55)

### Per transect
for (t in 1:length(unique(plot.data$Site))){
  transect = unique(plot.data$Site)[t]
  dat = plot.data %>% filter(Site == transect) %>% filter(ScenarioLab=="High")
  
  plot.out = ggplot(dat) + 
    geom_col(aes(x=TSFCATlab, y=mean_prop_landscape, fill=TSFCATlab)) + ylim(c(0,1)) + 
    theme_cowplot() + scale_fill_viridis_d() + theme(legend.position="none") + 
    ylab("Proportion of Landscape") + xlab("Time Since Fire")
  plot.out
  
  ggsave(plot.out,device="png", 
         path=here::here("Outputs"),
         width = 2.5,height =1.25, units = "in", dpi = 320, scale=2.25,
         filename = paste0("transectoptimisation_highsev_",transect,".png"))
}
