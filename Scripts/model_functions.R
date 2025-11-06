# Abundance models

library(nimble)


# Model One - Linear
run_linear_model = function(input_data){
  # Convert to radians
  deg2rad <- function(deg) {(deg * pi) / (180)}
  input_data$date_rad = deg2rad(input_data$date*360)
  
  linear_model_code = nimbleCode({ 
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
    }
    # Specify MVN prior for random site effects in lambda for each species
    for (i in 1:nsites){
      eta.lam[i,1:nspec] ~ dmnorm(mu.eta[1:nspec], Omega[,])
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
          beta1[k] * rainfall[i] + 
          beta2[k] * twi[i] +
          beta3[k] * bait[i] + # Bait Intensity
          
          
          # Option 1 - Linear
          beta4[k] * tsf[i] +
          beta5[k] * tsf[i] * bait[i] + # Interaction between time since fire and baiting
          beta6[k] * propsev[i] + # Proportion burnt severely
          eta.lam[i,k]
        
        # Chi-squared for N
        N_exp[i, k] <- lambda[i, k]
        chi2_N_obs[i, k] <- pow(N[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        N_sim[i, k] ~ dpois(lambda[i, k])
        chi2_N_sim[i, k] <- pow(N_sim[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        
        # Bernoilli RN model - observation model for replicated presence-absence data
        for (j in 1:nreps){
          C[i,j,k] ~ dbern(Pstar[i,j,k])
          Pstar[i,j,k] <- 1-pow((1 - p[i,j,k]), N[i,k])
          logit(p[i,j,k]) <- alpha0[k] + alpha1[k]*date_rad[i]
        }
      }
    }
    # Aggregate chi-squared statistics
    chi2_N_obs_total <- sum(chi2_N_obs[1:nsites, 1:nspec])
    chi2_N_sim_total <- sum(chi2_N_sim[1:nsites, 1:nspec])
  })
  
  # Parameters monitored
  params <- c('mean.lambda', 'mean.p', 'alpha0','alpha1',
              'beta0', 'beta1', 'beta2','beta3',
              'beta4', 'beta5','beta6',
              'eta.lam', 'Sigma2', 'rho', 'N', 'lambda',
              'chi2_N_obs_total', 'chi2_N_sim_total')
  ni <- 500000 ; nb <- 100000 ; nt <- 400 ; na = 10000; nc = 3
  #ni <- 400 ; nb <- 100 ; nt <- 3 ; na = 1; nc = 3
  
  # Initial values
  Nst <- maxC
  Nst[is.na(Nst)] <- 0
  Nst = Nst + 1
  modelInits <- function(){list(N = Nst, 
                                N_sim = Nst,
                                mean.lambda = rep(1, nspec), 
                                beta0 = rep(0, nspec),
                                beta1 = rep(0, nspec),
                                beta2 = rep(0, nspec),
                                beta3 = rep(0, nspec),
                                beta4 = rep(0, nspec),
                                beta5 = rep(0, nspec),
                                beta6 = rep(0, nspec),
                                alpha1 = rep(0, nspec),
                                mean.p = rep(0.5, nspec), 
                                eta.lam = array(1, dim = c(548, nspec)),
                                lambda = array(1, dim = c(548, nspec)),
                                Omega = diag(nspec))}
  
  linear_model_out <- nimbleMCMC(
    code = linear_model_code,
    constants = input_data, ## provide the combined data & constants as constants
    inits = modelInits,
    monitors = params,
    niter = ni,
    nburnin = nb,
    nchains = nc,
    thin = nt,
    samples=TRUE,
    samplesAsCodaMCMC = TRUE,
    WAIC = TRUE)
  
  return(linear_model_out)
}

# Model One - Linear
run_linear_noint_model = function(input_data){
  # Convert to radians
  deg2rad <- function(deg) {(deg * pi) / (180)}
  input_data$date_rad = deg2rad(input_data$date*360)
  
  linear_model_code = nimbleCode({ 
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
    }
    # Specify MVN prior for random site effects in lambda for each species
    for (i in 1:nsites){
      eta.lam[i,1:nspec] ~ dmnorm(mu.eta[1:nspec], Omega[,])
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
          beta1[k] * rainfall[i] + 
          beta2[k] * twi[i] +
          beta3[k] * bait[i] + # Bait Intensity
          
          
          # Option 1 - Linear - no interaction
          beta4[k] * tsf[i] +
          beta5[k] * propsev[i] + # Proportion burnt severely
          eta.lam[i,k]
        
        # Chi-squared for N
        N_exp[i, k] <- lambda[i, k]
        chi2_N_obs[i, k] <- pow(N[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        N_sim[i, k] ~ dpois(lambda[i, k])
        chi2_N_sim[i, k] <- pow(N_sim[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        
        # Bernoilli RN model - observation model for replicated presence-absence data
        for (j in 1:nreps){
          C[i,j,k] ~ dbern(Pstar[i,j,k])
          Pstar[i,j,k] <- 1-pow((1 - p[i,j,k]), N[i,k])
          logit(p[i,j,k]) <- alpha0[k] + alpha1[k]*date_rad[i]
        }
      }
    }
    # Aggregate chi-squared statistics
    chi2_N_obs_total <- sum(chi2_N_obs[1:nsites, 1:nspec])
    chi2_N_sim_total <- sum(chi2_N_sim[1:nsites, 1:nspec])
  })
  
  # Parameters monitored
  params <- c('mean.lambda', 'mean.p', 'alpha0','alpha1',
              'beta0', 'beta1', 'beta2','beta3',
              'beta4', 'beta5',
              'eta.lam', 'Sigma2', 'rho', 'N', 'lambda',
              'chi2_N_obs_total', 'chi2_N_sim_total')
  ni <- 500000 ; nb <- 100000 ; nt <- 400 ; na = 10000; nc = 3
  #ni <- 400 ; nb <- 100 ; nt <- 3 ; na = 1; nc = 3
  
  # Initial values
  Nst <- maxC
  Nst[is.na(Nst)] <- 0
  Nst = Nst + 1
  modelInits <- function(){list(N = Nst, 
                                N_sim = Nst,
                                mean.lambda = rep(1, nspec), 
                                beta0 = rep(0, nspec),
                                beta1 = rep(0, nspec),
                                beta2 = rep(0, nspec),
                                beta3 = rep(0, nspec),
                                beta4 = rep(0, nspec),
                                beta5 = rep(0, nspec),
                                alpha1 = rep(0, nspec),
                                mean.p = rep(0.5, nspec), 
                                eta.lam = array(1, dim = c(548, nspec)),
                                lambda = array(1, dim = c(548, nspec)),
                                Omega = diag(nspec))}
  
  linear_noint_model_out <- nimbleMCMC(
    code = linear_model_code,
    constants = input_data, ## provide the combined data & constants as constants
    inits = modelInits,
    monitors = params,
    niter = ni,
    nburnin = nb,
    nchains = nc,
    thin = nt,
    samples=TRUE,
    samplesAsCodaMCMC = TRUE,
    WAIC = TRUE)
  
  return(linear_noint_model_out)
}

# Model Two - Square Root Transformation
run_sqrt_model = function(input_data){
  deg2rad <- function(deg) {(deg * pi) / (180)}
  
  input_data$date_rad = deg2rad(input_data$date*360)
  
  sqrt_model_code = nimbleCode({ 
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
    }
    # Specify MVN prior for random site effects in lambda for each species
    for (i in 1:nsites){
      eta.lam[i,1:nspec] ~ dmnorm(mu.eta[1:nspec], Omega[,])
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
          beta1[k] * rainfall[i] + 
          beta2[k] * twi[i] +
          beta3[k] * bait_log[i] + # Bait Intensity
          
          # Option 2 - Sqrt Transformation
          beta4[k] * tsf_sqrt[i] + 
          
          beta5[k] * tsf_sqrt[i] * bait[i] + # Interaction between time since fire and baiting
          beta6[k] * propsev[i] + # Proportion burnt severely
          eta.lam[i,k]
        
        # Chi-squared for N
        N_exp[i, k] <- lambda[i, k]
        chi2_N_obs[i, k] <- pow(N[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        N_sim[i, k] ~ dpois(lambda[i, k])
        chi2_N_sim[i, k] <- pow(N_sim[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        
        # Bernoilli RN model - observation model for replicated presence-absence data
        for (j in 1:nreps){
          C[i,j,k] ~ dbern(Pstar[i,j,k])
          Pstar[i,j,k] <- 1-pow((1 - p[i,j,k]), N[i,k])
          logit(p[i,j,k]) <- alpha0[k] + alpha1[k]*date_rad[i]
        }
      }
    }
    # Aggregate chi-squared statistics
    chi2_N_obs_total <- sum(chi2_N_obs[1:nsites, 1:nspec])
    chi2_N_sim_total <- sum(chi2_N_sim[1:nsites, 1:nspec])
  })
  
  # Parameters monitored
  params <- c('mean.lambda', 'mean.p', 'alpha0','alpha1',
              'beta0', 'beta1', 'beta2','beta3',
              'beta4', 'beta5','beta6',
              'eta.lam', 'Sigma2', 'rho', 'N', 'lambda',
              'chi2_N_obs_total', 'chi2_N_sim_total')
  ni <- 500000 ; nb <- 100000 ; nt <- 400 ; na = 10000; nc = 3
  #ni <- 400 ; nb <- 100 ; nt <- 3 ; na = 1; nc = 3
  
  # Initial values
  Nst <- maxC
  Nst[is.na(Nst)] <- 0
  Nst = Nst + 1
  modelInits <- function(){list(N = Nst, 
                                N_sim = Nst,
                                mean.lambda = rep(1, nspec), 
                                beta0 = rep(0, nspec),
                                beta1 = rep(0, nspec),
                                beta2 = rep(0, nspec),
                                beta3 = rep(0, nspec),
                                beta4 = rep(0, nspec),
                                beta5 = rep(0, nspec),
                                beta6 = rep(0, nspec),
                                alpha1 = rep(0, nspec),
                                mean.p = rep(0.5, nspec), 
                                eta.lam = array(1, dim = c(548, nspec)),
                                lambda = array(1, dim = c(548, nspec)),
                                Omega = diag(nspec))}
  
  sqrt_model_out <- nimbleMCMC(
    code = sqrt_model_code,
    constants = input_data, ## provide the combined data & constants as constants
    inits = modelInits,
    monitors = params,
    niter = ni,
    nburnin = nb,
    nchains = nc,
    thin = nt,
    samples=TRUE,
    samplesAsCodaMCMC = TRUE,
    WAIC = TRUE)
  
  return(sqrt_model_out)
}


# Model Three - Polynomial
run_poly_model = function(input_data){
  deg2rad <- function(deg) {(deg * pi) / (180)}
  
  input_data$date_rad = deg2rad(input_data$date*360)
  
  poly_model_code = nimbleCode({ 
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
      eta.lam[i,1:nspec] ~ dmnorm(mu.eta[1:nspec], Omega[,])
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
          beta1[k] * rainfall[i] + 
          beta2[k] * twi[i] +
          beta3[k] * bait[i] + # Bait Intensity
          # Option 3 - Quadtractic term
          beta4[k] * tsf[i] + beta5[k] * pow(tsf[i],2) +
          beta6[k] * tsf_sqrt[i] * bait[i] + # Interaction between time since fire and baiting
          beta7[k] * propsev[i] + # Proportion burnt severely
          eta.lam[i,k]
        
        # Chi-squared for N
        N_exp[i, k] <- lambda[i, k]
        chi2_N_obs[i, k] <- pow(N[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        N_sim[i, k] ~ dpois(lambda[i, k])
        chi2_N_sim[i, k] <- pow(N_sim[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        
        # Bernoilli RN model - observation model for replicated presence-absence data
        for (j in 1:nreps){
          C[i,j,k] ~ dbern(Pstar[i,j,k])
          Pstar[i,j,k] <- 1-pow((1 - p[i,j,k]), N[i,k])
          logit(p[i,j,k]) <- alpha0[k] + alpha1[k]*date_rad[i]
        }
      }
    }
    # Aggregate chi-squared statistics
    chi2_N_obs_total <- sum(chi2_N_obs[1:nsites, 1:nspec])
    chi2_N_sim_total <- sum(chi2_N_sim[1:nsites, 1:nspec])
  })
  
  # Parameters monitored
  params <- c('mean.lambda', 'mean.p', 'alpha0','alpha1',
              'beta0', 'beta1', 'beta2','beta3',
              'beta4','beta5','beta6', 'beta7',
              'eta.lam', 'Sigma2', 'rho', 'N', 'lambda',
              'chi2_N_obs_total', 'chi2_N_sim_total')
  ni <- 500000 ; nb <- 100000 ; nt <- 400 ; na = 10000; nc = 3
  #ni <- 400 ; nb <- 100 ; nt <- 3 ; na = 1; nc = 3
  
  # Initial values
  Nst <- maxC
  Nst[is.na(Nst)] <- 0
  Nst = Nst + 1
  modelInits <- function(){list(N = Nst, 
                                N_sim = Nst,
                                mean.lambda = rep(1, nspec), 
                                beta0 = rep(0, nspec),
                                beta1 = rep(0, nspec),
                                beta2 = rep(0, nspec),
                                beta3 = rep(0, nspec),
                                beta4 = rep(0, nspec),
                                beta5 = rep(0, nspec),
                                beta6 = rep(0, nspec),
                                beta7 = rep(0, nspec),
                                alpha1 = rep(0, nspec),
                                mean.p = rep(0.5, nspec), 
                                eta.lam = array(1, dim = c(548, nspec)),
                                lambda = array(1, dim = c(548, nspec)),
                                Omega = diag(nspec))}
  
  poly_model_out <- nimbleMCMC(
    code = poly_model_code,
    constants = input_data, ## provide the combined data & constants as constants
    inits = modelInits,
    monitors = params,
    niter = ni,
    nburnin = nb,
    nchains = nc,
    thin = nt,
    samples=TRUE,
    samplesAsCodaMCMC = TRUE,
    WAIC = TRUE)
  
  return(poly_model_out)
}


# Model Four - Splines
# TSF Splines


run_spline_model = function(input_data){
  deg2rad <- function(deg) {(deg * pi) / (180)}
  
  input_data$tsf_spline = splines::ns(covs_scaled$tsf.point, df=2)
  input_data$date_rad = deg2rad(input_data$date*360)
  
  
  
 spline_model_code = nimbleCode({ 
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
      beta4_1[k] ~ dnorm(0, 0.1)
      beta4_2[k] ~ dnorm(0, 0.1)
      beta5_1[k] ~ dnorm(0, 0.1)
      beta5_2[k] ~ dnorm(0, 0.1)
      beta6[k] ~ dnorm(0, 0.1)
    }
    # Specify MVN prior for random site effects in lambda for each species
    for (i in 1:nsites){
      eta.lam[i,1:nspec] ~ dmnorm(mu.eta[1:nspec], Omega[,])
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
          beta1[k] * rainfall[i] + 
          beta2[k] * twi[i] +
          beta3[k] * bait[i] + # Bait Intensity
          
          # Option 4 - Splines
          beta4_1[k] * tsf_spline[i, 1] + beta4_2[k] * tsf_spline[i, 2] + 
          beta5_1[k] * tsf_spline[i,1] * bait[i] + beta5_2[k] * tsf_spline[i,2] * bait[i] + # Interaction between time since fire and baiting
          beta6[k] * propsev[i] + # Proportion burnt severely
          eta.lam[i,k]
        
        # Chi-squared for N
        N_exp[i, k] <- lambda[i, k]
        chi2_N_obs[i, k] <- pow(N[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        N_sim[i, k] ~ dpois(lambda[i, k])
        chi2_N_sim[i, k] <- pow(N_sim[i, k] - N_exp[i, k], 2) / N_exp[i, k]
        
        # Bernoilli RN model - observation model for replicated presence-absence data
        for (j in 1:nreps){
          C[i,j,k] ~ dbern(Pstar[i,j,k])
          Pstar[i,j,k] <- 1-pow((1 - p[i,j,k]), N[i,k])
          logit(p[i,j,k]) <- alpha0[k] + alpha1[k]*date_rad[i]
        }
      }
    }
    # Aggregate chi-squared statistics
    chi2_N_obs_total <- sum(chi2_N_obs[1:nsites, 1:nspec])
    chi2_N_sim_total <- sum(chi2_N_sim[1:nsites, 1:nspec])
  })
  
  # Parameters monitored
  params <- c('mean.lambda', 'mean.p', 'alpha0','alpha1',
              #'alpha1_cos', 'alpha1_sin',
              'beta0', 'beta1', 'beta2','beta3',
              'beta4_1', 
              'beta4_2', 
              'beta5_1', 
              'beta5_2', 
              'beta6',
              'eta.lam', 'Sigma2', 'rho', 'N', 'lambda',
              'chi2_N_obs_total', 'chi2_N_sim_total')
  ni <- 500000 ; nb <- 100000 ; nt <- 400 ; na = 10000; nc = 3
  #ni <- 400 ; nb <- 100 ; nt <- 3 ; na = 1; nc = 3
  
  # Initial values
  Nst <- maxC
  Nst[is.na(Nst)] <- 0
  Nst = Nst + 1
  modelInits <- function(){list(N = Nst, 
                                N_sim = Nst,
                                mean.lambda = rep(1, nspec), 
                                beta0 = rep(0, nspec),
                                beta1 = rep(0, nspec),
                                beta2 = rep(0, nspec),
                                beta3 = rep(0, nspec),
                                beta4_1 = rep(0, nspec),
                                beta4_2 = rep(0, nspec),
                                beta5_1 = rep(0, nspec),
                                beta5_2 = rep(0, nspec),
                                beta6 = rep(0, nspec),
                                alpha1 = rep(0, nspec),
                                mean.p = rep(0.5, nspec), 
                                eta.lam = array(1, dim = c(548, nspec)),
                                lambda = array(1, dim = c(548, nspec)),
                                Omega = diag(nspec))}
  
  spline_model_out <- nimbleMCMC(
    code = spline_model_code,
    constants = input_data, ## provide the combined data & constants as constants
    inits = modelInits,
    monitors = params,
    niter = ni,
    nburnin = nb,
    nchains = nc,
    thin = nt,
    samples=TRUE,
    samplesAsCodaMCMC = TRUE,
    WAIC = TRUE)
  return(spline_model_out)
}