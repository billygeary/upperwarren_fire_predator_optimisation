library(prioritizr)
library(tidyverse)
library(cowplot)

#### Plot Results of Scenarios for Paper ####
scenarios = intersect(list.files(path = "Data_Processing", pattern = "abundance_scenarios.RData", full = TRUE),
                    list.files(path = "Data_Processing", pattern = "Results_NoStd", full = TRUE))
scenarios
s.nos = c(10, 4, 7)

#### Regional Fire Age Distributions ####
extract_optimal_firehist = function(solution.out){
  solution.summary = eval_n_summary(solution.out$problem, solution.out$solved)
  solution.summary$tsf = parse_number(solution.summary$summary)
  return(solution.summary)
}

results = data.frame()
for (s in s.nos){
  df = readRDS(scenarios[s])
  name = gsub("Data_Processing/Results_NoStd_UpperWarrenAbundance_","", scenarios[s])
  name = gsub("_abundance_scenarios.RData","", name)
  scen.sum = lapply(df, function(x){extract_optimal_firehist(x)})
  scen.sum.dt = data.table::rbindlist(scen.sum, idcol=TRUE)
  summaries.out = scen.sum.dt
  summaries.out$TSFCAT[summaries.out$tsf < 6] = "0-5 years"
  summaries.out$TSFCAT[summaries.out$tsf > 5 & summaries.out$tsf < 11] = "6-10 years"
  summaries.out$TSFCAT[summaries.out$tsf > 10 & summaries.out$tsf < 16] = "11-15 years"
  summaries.out$TSFCAT[summaries.out$tsf > 15 & summaries.out$tsf < 26] = "16-24 years"
  summaries.out$TSFCAT[summaries.out$tsf > 24] = ">24 years"
  summaries.out$Scenario = name
  results = rbind(results, summaries.out)
}

write.csv(results, "Data_Clean/tsfprops_results_highsev.csv")

summaries.out = read.csv("Data_Clean/tsfprops_results_highsev.csv")

plot = summaries.out %>% drop_na() %>% 
  group_by(Scenario,.id, TSFCAT) %>% summarise(cost = sum(cost)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=TSFCAT, y=cost/11, fill = Scenario)) + #facet_wrap(~Baiting, ncol = 1) +
  theme_cowplot() + scale_fill_viridis_d() +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
plot

plot.quantiles = summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id, TSFCAT) %>% summarise(cost = sum(cost)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years")),
         Scenario = as.factor(Scenario)) %>%
  group_by(Scenario,TSFCAT) %>% summarise(Median = quantile(cost, probs=0.5),
                                 Mean = mean(cost),
                                          LCI = quantile(cost, probs=0.1),
                                          UCI = quantile(cost, probs=0.90))

facet.labs =  c(s_bL_sH="Low", s_bA_sH="Medium", s_bAG_sH="High")
  
plot = summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id, TSFCAT) %>% summarise(cost = sum(cost)) %>%
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
plot  

ggsave(plot = plot, filename="outcome_plot_highsev.pdf", path="Outputs", device="pdf", width=4,height=6,units="in",scale=1.5)


#### Regional Abundance Outcomes ####
extract_abundances = function(solution.out){
  rep_summary = eval_feature_representation_summary(solution.out$problem, solution.out$solved)
  rep_summary$tsf = parse_number(rep_summary$summary)
  rep_summary = left_join(rep_summary, spp.lookup)
  return(rep_summary)
}

spp.lookup = data.frame(Species = c("Woylie", "Chuditch", "Numbat", "Quenda", "Koomal"),
                        feature = as.factor(1:5))


abundance.results = data.frame()
for (s in s.nos){
  print(paste(s,"started at", Sys.time()))
  df = readRDS(scenarios[s])
  name = gsub("Data_Processing/Results_NoStd_UpperWarrenAbundance_","", scenarios[s])
  name = gsub("_abundance_scenarios.RData","", name)
  abundance.sum = lapply(df, function(x){extract_abundances(x)})
  abundance.sum.dt = data.table::rbindlist(abundance.sum, idcol=TRUE)
  summaries.out = abundance.sum.dt
  summaries.out$TSFCAT[summaries.out$tsf < 6] = "0-5 years"
  summaries.out$TSFCAT[summaries.out$tsf > 5 & summaries.out$tsf < 11] = "6-10 years"
  summaries.out$TSFCAT[summaries.out$tsf > 10 & summaries.out$tsf < 16] = "11-15 years"
  summaries.out$TSFCAT[summaries.out$tsf > 15 & summaries.out$tsf < 26] = "16-24 years"
  summaries.out$TSFCAT[summaries.out$tsf > 24] = ">24 years"
  summaries.out$Scenario = name
  abundance.results = rbind(abundance.results, summaries.out)
  print(paste(s,"done at", Sys.time()))
}

write.csv(abundance.results, "Data_Clean/abundance_results_highsev.csv")
ab.summaries.out = read.csv("Data_Clean/abundance_results_highsev.csv")

### Overall outcomes
facet.labs =  data.frame(Scenario = c("s_bL_sH", "s_bA_sH","s_bAG_sH"),
                ScenarioLab = c("Low", "Medium", "High"))

abundance.plot.quantiles = ab.summaries.out %>% 
  filter(summary=="overall") %>%
  left_join(facet.labs, by="Scenario") %>%
  group_by(ScenarioLab, Species) %>% summarise(Median = quantile(absolute_held, probs=0.5),
                                                   Mean = mean(absolute_held),
                                                   LCI = quantile(absolute_held, probs=0.1),
                                                   UCI = quantile(absolute_held, probs=0.9))

overall.plot = ab.summaries.out %>% filter(summary=="overall") %>%
  left_join(facet.labs, by="Scenario") %>%
  mutate(ScenarioLab = factor(ScenarioLab, levels = c("Low", "Medium", "High"))) %>%
  ggplot() +
  geom_jitter(aes(x=ScenarioLab, y = absolute_held, color=ScenarioLab),alpha=0.4) + ylim(c(0,9000)) +
  geom_pointrange(data = abundance.plot.quantiles, 
                  aes(x=ScenarioLab, 
                      y=Mean, 
                      ymin=LCI, 
                      ymax=UCI)) +
  facet_wrap(~Species, nrow=2) +
  theme_cowplot() + scale_color_viridis_d(begin=0.2) + 
  ylab("Abundance") + xlab("") + theme(legend.position="none")
overall.plot

ggsave(plot = overall.plot, filename="abundance_plot_highsev.pdf", path="Outputs", device="pdf", width=4,height=3,units="in",scale=2.2)


#### Transect Specific Fire Age Classes ####
library(ggplot2)
library(tidyverse)
plot = summaries.out %>% drop_na() %>%
  ggplot() + 
  geom_point(aes(y=cost/11,x=tsf)) + xlab("Proportion of Landscape)") +
  geom_smooth(aes(y=cost/11,x=tsf), 
              method="gam",
              method.args=list(family="binomial"))+
  facet_wrap(~Scenario)
plot









plot = summaries.out %>% drop_na() %>% 
  group_by(Scenario,.id, TSFCAT, Species) %>% summarise(absolute_held = sum(absolute_held)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=TSFCAT, y=absolute_held, fill = Scenario)) + facet_wrap(~Species, ncol = 3) +
  theme_cowplot() + scale_fill_viridis_d() +
  ylab("Abundance") + xlab("Time Since Fire Category")
plot

plot.quantiles = summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id,Species, TSFCAT) %>% summarise(absolute_held = sum(absolute_held)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years")),
         Scenario = as.factor(Scenario)) %>%
  group_by(Scenario,TSFCAT, Species) %>% summarise(Median = quantile(absolute_held, probs=0.5),
                                          Mean = mean(absolute_held),
                                          LCI = quantile(absolute_held, probs=0.05),
                                          UCI = quantile(absolute_held, probs=0.95))

facet.labs =  c(s_bL_sM="Low", s_bA_sM="Aerial", s_bAG_sM="Aerial+Ground")

plot = summaries.out %>% drop_na() %>% 
  group_by(Scenario, .id, TSFCAT, Species) %>% summarise(absolute_held = sum(absolute_held)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years")),
         Scenario = factor(Scenario, levels = c("s_bL_sM", "s_bA_sM", "s_bAG_sM"))) %>% 
  ggplot() + 
  #geom_jitter(aes(x=TSFCAT,y=absolute_held, colour= TSFCAT), alpha=0.4) + scale_color_viridis_d() + 
  geom_pointrange(data = plot.quantiles, 
                  aes(x=TSFCAT, 
                      y=Mean, 
                      ymin=LCI, 
                      ymax=UCI)) +
  theme_cowplot() + facet_grid(Species~Scenario,labeller=labeller(Scenario=facet.labs),scales="free_y") +
  theme(legend.position="none") +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
plot  





plot = abundance.sum.dt %>% na.omit() %>%
  ggplot() + geom_line(aes(x=tsf, y = absolute_held, colour = .id)) + 
  facet_wrap(~Species, scale = "free_y")
plot

plot = abundance.sum.dt %>% na.omit() %>%
  group_by(.id, TSFCAT, Species) %>% summarise(abundance = sum(absolute_held)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years"))) %>%
  group_by(TSFCAT, Species) %>% summarise(Median = quantile(abundance, probs=0.5),
                                          LCI = quantile(abundance, probs=0.1),
                                          UCI = quantile(abundance, probs=0.90)) %>%
  #separate(Scenario, into=c('Scen', 'Baiting', 'Severity')) %>%
  #filter(Baiting != "b0") %>%
  #mutate(Scenario = factor(Scenario, levels = c("Low", "Aerial", "Aerial + Ground (Monthly)", "Aerial + Ground (Fortnightly)", "Current"))) %>%
  ggplot() + 
  geom_pointrange(aes(x=TSFCAT, y=Median, 
                      ymin=LCI, 
                      ymax=UCI)) + 
  facet_wrap(~Species, ncol = 1, scales="free_y")
plot

