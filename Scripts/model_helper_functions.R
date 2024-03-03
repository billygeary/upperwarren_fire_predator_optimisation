# Function to extract covariates
clean_coefs_jags = function(b.name, coef.name,spp.list, model.list){
  
  b.lookup = data.frame(`.variable` = paste0(b.name,"[",1:nspec,"]"),
                        Species = spp.list,
                        Covariate = coef.name) 
  
  coefs = data.frame(Mean = model.list$mean[b.name],
                     Lower_2.5 = model.list$q2.5[b.name],
                     Upper_97.5 = model.list$q97.5[b.name])
  names(coefs) = c("Mean", "Lower_2.5", "Upper_97.5")
  
  var.out = cbind(b.lookup, coefs)
  
  var.out
}

clean_coefs_forest_jags = function(b.name,forest.type, coef.name,spp.list, model.list){
  
  b.lookup = data.frame(`.variable` = paste0(b.name,"[", forest.type,1:nspec,"]"),
                        Species = spp.list,
                        Covariate = coef.name) 
  
  coefs = data.frame(Mean = model.list$mean[[b.name]][,forest.type],
                     Lower_2.5 = model.list$q2.5[[b.name]][,forest.type],
                     Upper_97.5 = model.list$q97.5[[b.name]][,forest.type])
  names(coefs) = c("Mean", "Lower_2.5", "Upper_97.5")
  
  var.out = cbind(b.lookup, coefs)
  
  var.out
}

# Function to extract covariates
clean_coefs_90ci = function(b.name, coef.name,spp.list, model.list){
  b.lookup = data.frame(`.variable` = paste0(b.name,"[",1:nspec,"]"),
                        Species = spp.list,
                        Covariate = coef.name) 
  
  var.out = model.list %>%
    tidy_draws() %>%
    gather_variables() %>%
    filter(.variable %in% b.lookup$.variable) %>%
    left_join(b.lookup, by='.variable') %>%
    group_by(Species, Covariate) %>% summarise(Mean = mean(.value), 
                                               Lower = quantile(.value, probs = 0.055), 
                                               Upper= quantile(.value, probs = 0.945))
  var.out
}


#### Function to extract predictions from posterior coefficient estimates
predict.abundance.posterior = function(strategy, species.list, nsamp){
  predictions = list()
  for(n in 1:nsamp){
    pred.out = data.frame()
    for(i in 1:length(species.list)){ # Loop over each observed species
      b0 = model.betas %>% filter(beta == "beta0" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b1 = model.betas %>% filter(beta == "beta1" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b2 = model.betas %>% filter(beta == "beta2" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b3 = model.betas %>% filter(beta == "beta3" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b4 = model.betas %>% filter(beta == "beta4" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b5 = model.betas %>% filter(beta == "beta5" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b6 = model.betas %>% filter(beta == "beta6" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b7 = model.betas %>% filter(beta == "beta7" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      
      eta.lam.val = model.eta.lam %>% filter(species == i) %>% dplyr::select(mean) %>% unlist() %>% as.numeric()
      
      strat = strategy
      pred.df = data.frame(
        Sample = n,
        Scenario = rownames(strat),
        LocationName = strat$LocationName,
        Site = strat$Site,
        tsf_actual = strat$tsf_actual,
        Species = species.list[i],
        Abundance = exp(sample(b0, size=1) + 
                          sample(b1, size=1) * strat$rainfall +
                          sample(b2, size=1) * strat$twi +
                          sample(b3, size=1) * strat$bait +
                          sample(b4, size=1) * strat$tsf +
                          sample(b5, size=1) * strat$tsf^2 +
                          sample(b6, size=1) * strat$tsf * strat$bait +
                          sample(b7, size=1) * strat$propsev +
                          eta.lam.val
        )
      )
      pred.out = rbind(pred.out, pred.df)
    }
    pred.out = pivot_wider(pred.out, id_cols = c(Sample,Scenario, Site, LocationName, tsf_actual), names_from = Species, values_from = Abundance)
    covs = data.frame(tsf = strat$tsf)
    pred.out = cbind(covs, pred.out)
    predictions[[n]] <- pred.out
  }
  return(predictions)
}

#### Predict median abundance
predict.abundance.mean= function(strategy, species.list){
  pred.out = data.frame()
  for(i in 1:length(species.list)){ # Loop over each observed species
    strat = strategy
    pred.df = data.frame(
      LocationName = strat$LocationName,
      Site = strat$Site,
      Species = species.list[i],
      Abundance = exp(model_output$mean$beta0[i] + 
                        model_output$mean$beta1[i]  * strat$rainfall +
                        model_output$mean$beta2[i]  * strat$twi +
                        model_output$mean$beta3[i]  * strat$bait +
                        model_output$mean$beta4[i]  * strat$tsf +
                        model_output$mean$beta5[i]  * strat$tsf^2 +
                        model_output$mean$beta6[i]  * strat$tsf * strat$bait +
                        model_output$mean$beta7[i]  * strat$propsev +
                        model_output$mean$eta.lam[,i]
      )
    )
    pred.out = rbind(pred.out, pred.df)
  }
  pred.out = pivot_wider(pred.out, id_cols = c(Site, LocationName), names_from = Species, values_from = Abundance)
  covs = data.frame(tsf = strat$tsf, tsf_actual = strat$tsf_actual)
  pred.out = cbind(covs, pred.out)
  
  return(pred.out)
}

#### Function to extract predictions from posterior coefficient estimates
predict.occupancy.posterior = function(strategy, species.list, nsamp){
  predictions = list()
  for(n in 1:nsamp){
    pred.out = data.frame()
    for(i in 1:length(species.list)){ # Loop over each observed species
      b0 = model.betas %>% filter(beta == "beta0" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b1 = model.betas %>% filter(beta == "beta1" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b2 = model.betas %>% filter(beta == "beta2" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b3 = model.betas %>% filter(beta == "beta3" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b4 = model.betas %>% filter(beta == "beta4" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b5 = model.betas %>% filter(beta == "beta5" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b6 = model.betas %>% filter(beta == "beta6" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b7 = model.betas %>% filter(beta == "beta7" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b8 = model.betas %>% filter(beta == "beta8" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b9 = model.betas %>% filter(beta == "beta9" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      b10 = model.betas %>% filter(beta == "beta10" & species == i) %>% dplyr::select(.value) %>% unlist() %>% as.numeric()
      
      strat = strategy
      pred.df = data.frame(
        Sample = n,
        LocationName = strat$LocationName,
        Site = strat$Site,
        Species = species.list[i],
        Abundance = 1 / (1 + exp(-(sample(b0, size=1) + 
                                     sample(b1, size=1) * strat$east +
                                     sample(b2, size=1) * strat$north +
                                     sample(b3, size=1) * strat$twi +
                                     sample(b4, size=1) * strat$tsf +
                                     sample(b5, size=1) * strat$tsf^2 +
                                     sample(b6, size=1) * strat$propsev +
                                     sample(b7, size=1) * strat$bait + 
                                     sample(b8, size=1) * strat$tsf * strat$bait +
                                     sample(b9, size=1) * strat$propnv +
                                     sample(b10, size=1) * strat$hydro)))
      )
      pred.out = rbind(pred.out, pred.df)
    }
    pred.out = pivot_wider(pred.out, id_cols = c(Sample, Site, LocationName), names_from = Species, values_from = Abundance)
    covs = data.frame(tsf = strat$tsf)
    pred.out = cbind(covs, pred.out)
    predictions[[n]] <- pred.out
  }
  predictions <- do.call("rbind", predictions)
  return(predictions)
}


# Calculate GMA

calculate_gma = function(x) {
  site.abundances = x[,2:length(species.predict)] # A matrix where the rows are sites and the columns are species
  abundance.sums = colSums(site.abundances)
  gma = exp(mean(log(abundance.sums)))
  return(gma)
}

