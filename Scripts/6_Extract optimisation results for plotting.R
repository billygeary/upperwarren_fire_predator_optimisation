library(prioritizr)
library(tidyverse)
library(ggplot2)
library(cowplot)

#### Step 1: Load Custom Functions
results.paths = list.files(path = "Data_Processing", pattern = "Results", full = TRUE)

#### Step 4: Interrogate Results ####
# Plot the number of times a tsf category was chosem
extract_optimal_firehist = function(solution.out){
  solution.summary = eval_n_summary(solution.out$problem, solution.out$solved)
  solution.summary$tsf = parse_number(solution.summary$summary)
  solution.summary$TSFCAT[solution.summary$tsf < 6] = "0-5 years"
  solution.summary$TSFCAT[solution.summary$tsf > 5 & solution.summary$tsf < 11] = "6-10 years"
  solution.summary$TSFCAT[solution.summary$tsf > 10 & solution.summary$tsf < 16] = "11-15 years"
  solution.summary$TSFCAT[solution.summary$tsf > 15 & solution.summary$tsf < 26] = "16-24 years"
  solution.summary$TSFCAT[solution.summary$tsf > 24] = ">24 years"
  solution.summary$Scenario = name
  return(solution.summary)
}

# This takes a long time
summaries = list()
for (i in 1:length(results.paths)){
  name = gsub("Data_Processing/Results_UpperWarrenAbundance_","", results.paths[i])
  name = gsub("_abundance_scenarios.RData","", name)
  results = readRDS(results.paths[i])
  scen.sum = lapply(results, function(x){extract_optimal_firehist(x)})
  scen.sum.dt = data.table::rbindlist(scen.sum, idcol=TRUE);scen.sum.dt$Scenario = paste0(i)
  scen.sum.dt$Scenario = name
  summaries[[i]] <- scen.sum.dt
  print(i)
}

summaries.out = data.table::rbindlist(summaries, idcol=FALSE)

write.csv(summaries.out, "Data_Clean/baitfireoptimisationresults_transectscale_summaries.csv")
summaries.out = read.csv("Data_Clean/baitfireoptimisationresults_transectscale_summaries.csv")

head(summaries.out)

plot = summaries.out %>% drop_na() %>% 
  group_by(.id, Scenario, TSFCAT) %>% summarise(cost = sum(cost)) %>%
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
  group_by(.id, Scenario, TSFCAT) %>% summarise(cost = sum(cost)) %>%
  mutate(TSFCAT = factor(TSFCAT, levels = c("0-5 years", "6-10 years", "11-15 years", "16-24 years", ">24 years"))) %>%
  #separate(Scenario, into=c('Scen', 'Baiting', 'Severity')) %>%
  #filter(Baiting != "b0") %>%
  #mutate(Scenario = factor(Scenario, levels = c("Low", "Aerial", "Aerial + Ground (Monthly)", "Aerial + Ground (Fortnightly)", "Current"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=TSFCAT, y=cost/11, fill = TSFCAT)) + facet_wrap(~Scenario, ncol = 3) +
  theme_cowplot() + scale_fill_viridis_d() +
  ylab("Proportion of Landscape") + xlab("Time Since Fire Category")
plot
