##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
## 5: Find the Optimal Fire Age Distributions per Scenario ##
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++##

library(prioritizr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(gurobi)
library(terra)
library(pbmcapply)
#######################################
#### Step 1: Load Custom Functions ####
#######################################
# Clean scenarios 
clean_scenarios = function(scen, spp.list){
  scenarios = data.table::rbindlist(scen, idcol=FALSE)
  scenarios = scenarios %>% dplyr::select(tsf, Scenario, Sample, Site, LocationName, all_of(spp.list)) %>%
  mutate(scen = row_number())  # Assign row numbers directly
  
  scenarios.sample = scenarios %>% 
    separate_wider_delim(cols=Scenario, names=c("S", "Baiting","Severity", "TSF"), delim="_")
  
  pu.lookup = data.frame(LocationName = unique(scenarios.sample$LocationName))
  pu.lookup$pu = parse_number(pu.lookup$LocationName)
  scenarios.sample = left_join(scenarios.sample, pu.lookup, by = "LocationName")
  return(scenarios.sample)
}

# Convert to prioritizR
scenarios_to_priortizr = function(scenario.samples) {
  samples.list = list()
  for (n in 1:length(unique(scenario.samples$Sample))){ # Loop through the samples (n=1000)
    scenarios = scenario.samples %>% filter(Sample == n)
    tsf.val = unique(scenarios$tsf)
    scenarios.list = list()
    for(i in 1:length(tsf.val)){ # Loop through the tsf vals (33)
      df = subset(scenarios,tsf==tsf.val[i])
      stack.list = list()
      for (s in 9:(8+length(spp.to.optimise))){ # Loop through the species
        mat = df%>% pivot_wider(id_cols=pu, names_from=Site, values_from=all_of(s)) %>% dplyr::select(2:12) %>% as.matrix()
        mat = as.matrix(colSums(mat, na.rm=TRUE)) # Do we want to do at site or transect level
        ras = terra::rast(mat)
        stack.list[[colnames(df)[s]]] <- ras
      }
      stack = terra::rast(stack.list)
      scenarios.list[[paste0("tsf",i)]] <- stack
    }
    samples.list[[n]] <- scenarios.list
  }
  return(samples.list)
}

# Build prioritizR Zones
build.zones = function(x){lapply(x, function(y){do.call('zones',y)})}

# Create targets
calculate_custom_targets = function(scenarios){
  targets = scenarios %>% 
    group_by(LocationName) %>% 
    summarise(Woylie = quantile(Woylie, probs = 0.9), 
              Chuditch = quantile(Chuditch, probs=0.9), 
              Numbat = quantile(Numbat,probs = 0.9), 
              Quenda = quantile(Quenda, probs=0.9)) %>%
    dplyr::select(Chuditch, Quenda, Woylie, Numbat) %>% colSums()
  targets_list <- tibble::tibble(
    feature = c(paste0(1:4)), 
    zone = list(paste0("tsf",1:33), paste0("tsf",1:33), paste0("tsf",1:33), paste0("tsf",1:33)),
    sense = c(">=", ">=", ">=", ">="), 
    type = c("absolute", "absolute", "absolute", "absolute"),
    target = targets)
  
  return(targets_list)
}

# Solve the prioritizR problem
solve_prioritizr_problem = function(z, planning.units, targets){
  p <- problem(planning.units, z) %>%
    add_min_largest_shortfall_objective(budget = 10000000) %>%
    add_manual_targets(targets) %>% 
    add_proportion_decisions() %>% 
    add_mandatory_allocation_constraints() %>% 
    add_gurobi_solver(verbose = FALSE, 
                      numeric_focus = TRUE)
  solved.p = solve(p, 
                   force=TRUE, 
                   run_checks = FALSE)
  p.out = list(problem = p, solved = terra::wrap(solved.p)) 
  return(p.out)
}


# Big function for running through each of the steps above per scenario
run_optimisation_procedure = function(scenario.file, spp.to.optimise){
  scen = readRDS(scenario.file)
  cleaned.scenarios = clean_scenarios(scen, spp.to.optimise)
  print(paste("Scenarios cleaned"))
  # Convert the scenarios to raster grids for input into Prioritizr
  scen.p = scenarios_to_priortizr(cleaned.scenarios)
  print(paste("PrioritizR inputs built"))
  # Build Zones Raster Stack
  # Zones are management actions
  z.scen = build.zones(scen.p)
  print(paste("Zones built"))
  # Build planning unit raster (this assigns costs)
  pu = terra::rast(resolution = terra::res(scen.p[[1]]$tsf1$Woylie), 
                   ext = terra::ext(scen.p[[1]]$tsf1$Woylie),
                   vals=1, crs=NA)
  n <- 33  # number of zones
  pu_stack <- c(lapply(1:n, function(i) pu)) 
  pu_stack <- terra::rast(pu_stack) # Convert to SpatRaster stack
  #Calculate custom targets for scenario 
  targets_scenarios = calculate_custom_targets(cleaned.scenarios)
  print(paste("Targets created"))
  rm(cleaned.scenarios)
  gc()
  #Solve the optimisation
  pb = progressBar(min = 0, max = 1000, initial = 0, style = "ETA") 
  p.scenario = list()
  for (s in 1:1000){
    p.solved = solve_prioritizr_problem(z = z.scen[[s]], planning.units=pu_stack, targets = targets_scenarios)
    p.scenario[[s]] <- p.solved
    setTxtProgressBar(pb,s)
  }
  close(pb)
  print(paste("Optimisation done"))
  gc()
  return(p.scenario)
}

# Extract optimal fire history results
extract_optimal_firehist = function(solution.out){
  solution.summary = eval_n_summary(solution.out$problem, terra::unwrap(solution.out$solved))
  solution.summary$tsf = parse_number(solution.summary$summary)
  solution.summary$TSFCAT = NA
  solution.summary$TSFCAT[solution.summary$tsf < 6] = "0-5 years"
  solution.summary$TSFCAT[solution.summary$tsf > 5 & solution.summary$tsf < 11] = "6-10 years"
  solution.summary$TSFCAT[solution.summary$tsf > 10 & solution.summary$tsf < 16] = "11-15 years"
  solution.summary$TSFCAT[solution.summary$tsf > 15 & solution.summary$tsf < 26] = "16-24 years"
  solution.summary$TSFCAT[solution.summary$tsf > 24] = ">24 years"
  solution.summary$Scenario = name
  return(solution.summary)
}

# Extract the abundance outcome results
extract_abundances = function(solution.out){
  rep_summary = eval_feature_representation_summary(solution.out$problem, terra::unwrap(solution.out$solved))
  rep_summary$tsf = parse_number(rep_summary$summary)
  spp.lookup = data.frame(Species = c("Chuditch","Quenda", "Woylie", "Numbat"),
                          feature = as.factor(1:4))
  rep_summary = left_join(rep_summary, spp.lookup, by="feature")
  return(rep_summary)
}

# Extract the transect scale results
extract_transect_summaries = function(solution) {
  solution.out = terra::unwrap(solution$solved)
  summary.out = solution.out %>%
    as.matrix() %>% 
    as.data.frame() %>%
    mutate(Site = sites) %>%
    pivot_longer(cols = 1:33, names_to="tsf", values_to="prop_landscape") %>%
    mutate(tsf = as.numeric(gsub("tsf","",tsf)))
  
  summary.out$TSFCAT[summary.out$tsf < 6] = "0-5 years"
  summary.out$TSFCAT[summary.out$tsf > 5 & summary.out$tsf < 11] = "6-10 years"
  summary.out$TSFCAT[summary.out$tsf > 10 & summary.out$tsf < 16] = "11-15 years"
  summary.out$TSFCAT[summary.out$tsf > 15 & summary.out$tsf < 26] = "16-24 years"
  summary.out$TSFCAT[summary.out$tsf > 24] = ">24 years"
  
  return(summary.out)
}

################################
#### Step 2: Load Scenarios ####
################################
# The data is set out such that each column represents a species and
# each row represents a time since fire value for a given site. Each matrix is the abundance predictions at different
# sites.  
# the.data[i,j,k] then is the abundance at each time since fire i, of the j-th species in the k-th site. 

# TSF by species matrix
scenarios = setdiff(list.files(path = "Data_Processing", pattern = "abundance_scenarios.RData", full = TRUE),
                    list.files(path = "Data_Processing", pattern = "Results", full = TRUE))


scenarios.names = setdiff(list.files(path = "Data_Processing", pattern = "abundance_scenarios.RData", full = FALSE),
                          list.files(path = "Data_Processing", pattern = "Results", full = FALSE))

# Which species' abundance are we trying to maximise?
spp.to.optimise = c("Chuditch","Quenda", "Woylie", "Numbat")
sites = unique(readRDS(scenarios[1])[[1]]$Site)

######################################################################################
#### Step 3: Build, Run Optimisation Algorithm and Save Outputs for each Scenario ####
######################################################################################
# Takes about 0.75 hours to run each scenario through this loop

for (i in 1:length(scenarios)){
  # Do the optimisation - this takes 1.5-2 hours per iteration
  t1 = Sys.time()
  print(paste("Starting new scenario at", t1))
  scen.name = scenarios.names[i]
  results = run_optimisation_procedure(scenarios[i], spp.to.optimise)
  print(paste0("Scenario ", i, " optimisation complete after ", (Sys.time()-t1) ))
  
  # Make the name for saving 
  name = gsub(".RData", "", scen.name)
  name = gsub("UpperWarrenAbundance_", "",name)
  name = gsub("_abundance_scenarios", "", name)
  
  # Extract the optimisation results and save
  scen.sum =  suppressWarnings(lapply(results, function(x){extract_optimal_firehist(x)}))
  scen.sum.dt = data.table::rbindlist(scen.sum, idcol=TRUE)
  scen.sum.dt$Scenario = name
  write.csv(scen.sum.dt, paste0("Data_Clean/baitfireoptimisationresults_transectscale_summaries_", name, ".csv"))
  
  # Extract the abundance results and save
  abundance.sum = suppressWarnings(lapply(results, function(x){extract_abundances(x)}))
  abundance.sum.dt = data.table::rbindlist(abundance.sum, idcol=TRUE)
  abundance.sum.dt$Scenario = name
  write.csv(abundance.sum.dt, paste0("Data_Clean/baitfireoptimisationresults_transectscale_abundances_", name, ".csv"))
  
  # Extract the transect-scale results and save
  sites = sites
  transect.results = lapply(results, function(x){extract_transect_summaries(x)})
  transect.results.dt = data.table::rbindlist(transect.results, idcol=TRUE)
  transect.results.dt$Scenario = name
  
  write.csv(transect.results.dt, paste0("Data_Clean/baitfireoptimisationresults_transectscale_outcomes_", name, ".csv"))
  print(paste0("Scenario ", i, " saved after ", (Sys.time()-t1) ))
}

# Finished