#### JAGS Fire N Mixture Model ####
run_nmix_fire_jagscode = function(seed, data, nt, ni, nb, nc, na){
  y <- data$C
  nsites <- dim(y)[1]
  nreps <- dim(y)[2]
  nspec <- dim(y)[3]
  maxC <- apply(y, c(1,3), max, na.rm = TRUE)
  maxC[maxC == -Inf] <- NA
  
  cat(file = "nmix_fire_jags.txt", "
   model {
    # Priors
    # Intercepts and coefficients all fixed effects
    for(k in 1:nspec){
      mean.lambda[k] <- exp(beta0[k])
      beta0[k] ~ dnorm(0, 0.1)
      alpha0[k] <- logit(mean.p[k])
      alpha1[k] ~ dnorm(0, 0.1)
      mean.p[k] ~ dunif(0,1)
      beta1[k] ~ dnorm(0, 0.1)
      beta2[k] ~ dnorm(0, 0.1)
      beta3[k] ~ dnorm(0, 0.1)
      beta4[k] ~ dnorm(0, 0.1)
      beta5[k] ~ dnorm(0, 0.1)
      beta6[k] ~ dnorm(0, 0.1)      
      beta7[k] ~ dnorm(0, 0.1)
    }
    # Specify MVN prior for random site effects in lambda for each species
    for (i in 1:nsites){
      eta.lam[i,1:nspec] ~ dmnorm(mu.eta[], Omega[,])
    }
    for (k in 1:nspec){
      mu.eta[k] <- 0
    }
    # Vague inverse Wishart prior for variance-covariance matrix
    Omega[1:nspec,1:nspec] ~ dwish(R[,], df)
    Sigma2[1:nspec,1:nspec] <- inverse(Omega[,])
    
    # Scale var/covar matrix to become the correlation matrix
    for (i in 1:nspec){
      for (k in 1:nspec){
        rho[i,k] <- Sigma2[i,k] / (sqrt(Sigma2[i,i]) * sqrt(Sigma2[k,k]))
      }
    }
    # Likelihood
    # Ecological model for true abundance
    for (i in 1:nsites){
      for(k in 1:nspec){
        N[i,k] ~ dpois(lambda[i,k])
        log(lambda[i,k]) <- beta0[k] + 
          beta1[k] * east[i] + 
          beta2[k] * north[i] +
          beta3[k] * bait[i] + # Bait Intensity
          beta4[k] * tsf[i] + beta5[k] * pow(tsf[i],2) + # Time since fire with quadratic term
          beta6[k] * tsf[i] * bait[i] + # Interaction between time since fire and baiting
          beta7[k] * propsev[i] + # Proportion burnt severely
          eta.lam[i,k]
          
        # Observation model for replicated counts
        for (j in 1:nreps){
          C[i,j,k] ~ dbin(p[i,j,k], N[i,k])
          logit(p[i,j,k]) <- alpha0[k] + alpha1[k]*date[i]
        }
      }
    }
  }")
  
  
  # Initial values
  Nst <- maxC + 1
  Nst[is.na(Nst)] <- 1
  modelInits <- function(){list(N = Nst, 
                                mean.p = rep(0.2, nspec), 
                                Omega = diag(nspec))}
  
  # Parameters monitored
  params <- c('mean.lambda', 'mean.p', 'alpha0','alpha1',
              'beta0', 'beta1', 'beta2','beta3','beta4','beta5','beta6','beta7',
              'eta.lam', 'Sigma2', 'rho', 'N')
  
  
  # Run in Jags
  out <- jagsUI::jags(data, 
                      modelInits, 
                      params, 
                      "nmix_fire_jags.txt", 
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
  
  return(out)
}

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
                                               Lower = quantile(.value, probs = 0.05), 
                                               Upper= quantile(.value, probs = 0.95))
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
                          model_output$mean$eta.lam[,i]
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

