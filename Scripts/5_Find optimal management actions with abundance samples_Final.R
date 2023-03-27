library(prioritizr)
library(tidyverse)
library(ggplot2)
library(cowplot)

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
        ras = raster::raster(mat)
        stack.list[[colnames(df)[s]]] <- ras
      }
      stack = raster::stack(stack.list)
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
  p.out = list(problem = p, solved = solved.p) 
  return(p.out)
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

for (s in 1:length(scenarios)){
  print(paste("Starting new scenario"))
  scen = readRDS(scenarios[s])
  scen.name = scenarios.names[s]
  cleaned.scenarios = clean_scenarios(scen, spp.to.optimise)
  # Convert the scenarios to raster grids for input into Prioritizr
  scen.p = scenarios_to_priortizr(cleaned.scenarios)
  print(paste("PrioritizR inputs built"))
  # Build Zones Raster Stack
  # Zones are management actions
  z.scen = build.zones(scen.p)
  print(paste("Zones built"))
  # Build planning unit raster (this assigns costs)
  pu = raster(resolution = res(scen.p[[1]]$tsf1$Woylie), 
              ext = extent(scen.p[[1]]$tsf1$Woylie),
              vals=1, crs=NA)
  n <- 33  # number of zones
  pu_stack <- stack(lapply(1:n, function(i) pu)) 
  
  #Calculate custom targets for scenario 
  targets_scenarios = calculate_custom_targets(cleaned.scenarios)
  print(paste("Targets created"))
  rm(cleaned.scenarios)
  gc()
  #Solve the optimisation
  p.scenario = pbapply::pblapply(z.scen, solve_prioritizr_problem)
  print(paste("Optimisation done"))
  # Save to folder
  saveRDS(p.scenario, paste0("Data_Processing/Results_NoStd_", scen.name)) 
  rm(p.scenario)
  rm(scen.p)
  rm(z.scen)
  gc()
  print(paste0("Scenario ", s, " complete"))
}


