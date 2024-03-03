# Reproducible example for multi-species binomial N-mixture model issue. 
# Model based off example in 8.5.1 in AHM Vol 2
library(nimble)
dethist = readRDS(here::here("example_data.RData"))

# Get info about the data
nsite = nrow(dethist$Sites) # number of sites
nrep = max(sapply(dethist$Species, ncol)) # maximum number of replicate surveys per season
nspec = length(dethist$Species) # number of species

# Convert det histories to a 3d array
multispp.data = dethist$Species
spp.names = names(multispp.data)
multispp.data = array(unlist(multispp.data), 
                      dim = c(dim(multispp.data[[1]]), length(multispp.data)),
                      dimnames = list(site = rownames(multispp.data[[1]]), rep = colnames(multispp.data[[1]]), sps = names(multispp.data)))

dim(multispp.data) == c(nsite, nrep, nspec) # Check the dimensions match

#### Step 2: Setup Model Inputs ####
y <- multispp.data
class(y) <- 'numeric'
nsites <- dim(y)[1]
nreps <- dim(y)[2]
nspec <- dim(y)[3]
maxC <- apply(y, c(1,3), max, na.rm = TRUE)

# Standardise coefficients
covs_scaled = covs %>% dplyr::select(-c(landscape_position, forest_position)) %>%
  mutate(date = julian(dethist$Sites$Start, origin = as.Date("2016-10-17"))) %>%
  sapply(FUN = function(x) {as.numeric(scale(x))}) %>% as.data.frame()

# Compile into dataframe
Rmat <- diag(nspec)      # Identity matrix
df <- nspec + 1

# Bundle and summarize data set
bdata <- list(C = y, 
              nsites = nsites, 
              nspec = nspec,
              nreps = nreps,
              rainfall = covs_scaled$rainfall,
              date=covs_scaled$date,
              R = Rmat, 
              df = df)


model_code = nimbleCode({ 
  # Priors
  # Intercepts and coefficients all fixed effects
  for(k in 1:nspec){
    mean.lambda[k] <- exp(beta0[k])
    beta0[k] ~ dnorm(0, 0.1)
    alpha0[k] <- logit(mean.p[k])
    alpha1[k] ~ dnorm(0, 0.1)
    mean.p[k] ~ dunif(0,1)
    beta1[k] ~ dnorm(0, 0.1)
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
        eta.lam[i,k]
      
      # Observation model for replicated counts
      for (j in 1:nreps){
        C[i,j,k] ~ dbin(p[i,j,k], N[i,k])
        logit(p[i,j,k]) <- alpha0[k] + alpha1[k]*date[i]
      }
    }
  }})

# Parameters to monitor
params <- c('mean.lambda', 'mean.p', 'alpha0','alpha1',
            'beta0', 'beta1',  'eta.lam', 'Sigma2', 'rho', 'N')

# MCMC settings
ni <- 1000 ; nb <- 500; nt <- 1 ; nc=3 

# Initial values
Nst <- maxC
Nst[is.na(Nst)] <- 0
Nst = Nst + 1
modelInits <- function(){list(N = Nst, 
                              mean.lambda = rep(1, 6), 
                              Omega = diag(6))}

model_out <- nimbleMCMC(
  code = model_code,
  constants = bdata, ## provide the combined data & constants as constants
  inits = modelInits,
  monitors = params,
  niter = ni,
  nburnin = nb,
  nchains = nc,
  thin = nt,
  samplesAsCodaMCMC=TRUE)

# Check the Rhat vals so we can quickly see which parameters have been sampled or not (NA vals)
library(MCMCvis)
sums = MCMCsummary(model_out)
summary(is.na(sums$Rhat))


