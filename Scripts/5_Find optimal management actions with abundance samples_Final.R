##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
## 5: Find the Optimal Fire Age Distributions per Scenario ##
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++##

library(prioritizr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(gurobi)
library(terra)

#### Step 1: Load Custom Functions
# Clean scenarios 
clean_scenarios = function(scen, spp.list){
  scenarios.list = list()
  scenarios = data.table::rbindlist(scen, idcol=FALSE)
  scenarios = scenarios %>% dplyr::select(tsf, Scenario, Sample, Site, LocationName, all_of(spp.list))  # Remove foxes 
  scenarios$scen = rownames(scenarios)
  
  scenarios.sample = scenarios %>% 
    separate(Scenario, into=c("S", "Baiting","Severity", "TSF"), sep="_") #%>%
    # mutate(Woylie = Woylie/max(Woylie),
    #        Chuditch = Chuditch/max(Chuditch),
    #        Numbat = Numbat/max(Numbat),
    #        Quenda = Quenda/max(Quenda),
    #        Koomal = Koomal/max(Koomal)) # Standardize the abundances so they are comparable between species
  
  pu.lookup = data.frame(LocationName = unique(scenarios.sample$LocationName))
  pu.lookup$pu = tidyr::extract_numeric(pu.lookup$LocationName)
  scenarios.sample = left_join(scenarios.sample, pu.lookup)
  return(scenarios.sample)
}

# Convert to prioritizR
scenarios_to_priortizr = function(scenario.samples) {
  samples.list = list()
  for (n in 1:length(unique(scenario.samples$Sample))){
    scenarios = scenario.samples %>% filter(Sample == n)
    
    scenarios.list = list()
    for(i in 1:length(unique(scenarios$tsf))){
      tsf.val = unique(scenarios$tsf)[i]
      df = subset(scenarios,tsf==tsf.val)
      stack.list = list()
      for (s in 9:(8+length(spp.to.optimise))){
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
    summarise(Woylie = max(Woylie), Chuditch = max(Chuditch), 
              Numbat = max(Numbat), Quenda = max(Quenda), Koomal = max(Koomal)) %>%
    dplyr::select(Woylie, Chuditch, Numbat, Quenda, Koomal) %>% colSums()
  
  targets_list <- tibble::tibble(
    feature = c(paste0(1:5)), 
    zone = list(paste0("tsf",1:33), paste0("tsf",1:33), paste0("tsf",1:33), paste0("tsf",1:33), paste0("tsf",1:33)),
    sense = c(">=", ">=", ">=", ">=", ">="), 
    type = c("absolute", "absolute", "absolute", "absolute", "absolute"),
    target = targets)
  
  return(targets_list)
}

# Solve the prioritizR problem
solve_prioritizr_problem = function(z){
  p <- problem(pu_stack, z) %>%
    add_min_largest_shortfall_objective(budget = 10000000) %>%
    add_manual_targets(targets_scenarios) %>% 
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

#### Step 2: Load Scenarios ####
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
spp.to.optimise = c("Woylie", "Chuditch", "Numbat", "Quenda", "Koomal")
spp.list = spp.to.optimise

#### Step 3: Build and Run Optimisation Algorithm for each Scenario ####
# Takes about two hours to run each scenario through this loop 

s = 12
  t1 = Sys.time()
  print(paste("Starting new scenario at", t1))
  scen = readRDS(scenarios[s])
  scen.name = scenarios.names[s]
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
  p.scenario = pbapply::pblapply(z.scen, solve_prioritizr_problem)
  print(paste("Optimisation done"))
  # Save to folder
  #filename = gsub(".RData", ".rds", scen.name)
  #saveRDS(p.scenario, paste0("Data_Processing/Results_NoStd_", scen.name)) 
  gc()
  print(paste0("Scenario ", s, " complete after ", (Sys.time()-t1) ))

  
results = p.scenario
name = gsub(".RData", "", scen.name)
name = gsub("UpperWarrenAbundance_", "",name)
name = gsub("_abundance_scenarios", "", name)
name


# Extract the optimisation results and save
scen.sum = lapply(results, function(x){extract_optimal_firehist(x)})
scen.sum.dt = data.table::rbindlist(scen.sum, idcol=TRUE)
scen.sum.dt$Scenario = name
write.csv(scen.sum.dt, paste0("Data_Clean/baitfireoptimisationresults_transectscale_summaries_", name, ".csv"))

# Extract the abundance outcome results and save
extract_abundances = function(solution.out){
  rep_summary = eval_feature_representation_summary(solution.out$problem, terra::unwrap(solution.out$solved))
  rep_summary$tsf = parse_number(rep_summary$summary)
  spp.lookup = data.frame(Species = c("Woylie", "Chuditch", "Numbat", "Quenda", "Koomal"),
                          feature = as.factor(1:5))
  rep_summary = left_join(rep_summary, spp.lookup)
  return(rep_summary)
}

abundance.sum = lapply(results, function(x){extract_abundances(x)})
abundance.sum.dt = data.table::rbindlist(abundance.sum, idcol=TRUE)
abundance.sum.dt$Scenario = name
write.csv(abundance.sum.dt, paste0("Data_Clean/baitfireoptimisationresults_transectscale_abundances_", name, ".csv"))

# Extract the transect scale results
sites = unique(scen[[1]]$Site)
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

transect.results = lapply(results, function(x){extract_transect_summaries(x)})
transect.results.dt = data.table::rbindlist(transect.results, idcol=TRUE)
transect.results.dt$Scenario = name

write.csv(transect.results.dt, paste0("Data_Clean/baitfireoptimisationresults_transectscale_outcomes_", name, ".csv"))


rm(cleaned.scenarios, p.scenario, scen.p, pu_stack, targets_scenarios, z.scen)
gc()
print(paste0("Scenario ", s, " saved after ", (Sys.time()-t1) ))

# Finished
