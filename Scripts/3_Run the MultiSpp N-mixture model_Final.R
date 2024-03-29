##+++++++++++++++++++++++++++++++++++++##
## 3: Fit the Multispp N-Mixture Model ##
##+++++++++++++++++++++++++++++++++++++##

library(nimble)
library(tidyverse)
# Read in data
dethist = readRDS("Data_Processing/camtrapR.counthist.UpperWarren.multispp09122023.RData")

# Select species we want to include in the N-Mixture Model
names(dethist$Species)
# dethist$Species =  dethist$Species[c("Chuditch", "Koomal", "Quenda", "Woylie", "Vulpes", "Numbat",
#                                      #"Roo","Tammar",  "Western Brush Wallaby", "Dunnart"
#                                      )]
dethist$Species =  dethist$Species[c("Chuditch","Quenda", "Woylie", "Vulpes", "Numbat")]
spp = names(dethist$Species)

dethist$Species = lapply(dethist$Species, function(x) x[,1:28]) # Drop last two survey occassions

# Check for colinearity in the covariates
covs = dethist$Sites %>% 
  dplyr::select(X,Y,10:32) %>% dplyr::select(-c(Study.AreaName, Lag))
cor.covs = cor(covs, use='complete.obs', method='pearson')
corrplot::corrplot(cor.covs)

high.cors = cor.covs %>% as.table() %>% as.data.frame() %>% filter(Freq < -0.7 | Freq >0.7) %>% filter(Freq !=1) %>%
  group_by(Var1) %>% summarise(Variables=paste(Var2, collapse=", "))

# Get info about the data
nsite = nrow(dethist$Sites) # number of sites
nrep = max(sapply(dethist$Species, ncol)) # maximum number of replicate surveys per season
nspec = length(dethist$Species)

# Convert det histories to a 3d array
multispp.data = dethist$Species
spp.names = names(multispp.data)
multispp.data = array(unlist(multispp.data), 
                      dim = c(dim(multispp.data[[1]]), length(multispp.data)),
                      dimnames = list(site = rownames(multispp.data[[1]]), rep = colnames(multispp.data[[1]]), sps = names(multispp.data)))

dim(multispp.data) == c(nsite, nrep, nspec) # Check the dimensions match

hist(multispp.data) # Check the counts

#### Step 2: Setup Model Inputs ####
y <- multispp.data
class(y) <- 'numeric'
nsites <- dim(y)[1]
nreps <- dim(y)[2]
nspec <- dim(y)[3]
maxC <- apply(y, c(1,3), max, na.rm = TRUE)
maxC[maxC == -Inf] <- NA

# Standardise coefficients
covs_scaled = covs %>% dplyr::select(-c(landscape_position, forest_position)) %>%
  mutate(date = julian(dethist$Sites$Start, origin = as.Date("2016-10-17"))) %>%
  sapply(FUN = function(x) {as.numeric(scale(x))}) %>% as.data.frame()

covs_scaled$landscape_position = covs$landscape_position
covs_scaled$forest_position = covs$forest_position
# Compile into dataframe
Rmat <- diag(nspec)      # Identity matrix
df <- nspec + 1

# Lookup
transects= data.frame(Site = unique(dethist$Sites$Site),
                      SiteNo = 1:length(unique(dethist$Sites$Site)))

dethist$Sites = left_join(dethist$Sites, transects, by="Site")

# Bundle and summarize data set
bdata <- list(C = y, 
              nsites = nsites, 
              nspec = nspec,
              nreps = nreps,
              #transect = dethist$Sites$SiteNo,
              #east = covs_scaled$X,
              #north = covs_scaled$Y,
              propsev = covs_scaled$PropSevere500,
              bait = covs_scaled$Mean_Intensity_400,
              tsf=covs_scaled$tsf.point,
              #ag = covs_scaled$prop_ag_3km,
              #propnv = covs_scaled$prop_nv_3km,
              twi = covs_scaled$twi,
              #hydro = covs_scaled$dist_to_majorhydro,
              #road = covs_scaled$prop_filtered_roads_3km,
              rainfall = covs_scaled$rainfall,
              date=covs_scaled$date,
              #landscape = covs_scaled$forest_position,
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
  }})

# Parameters monitored
params <- c('mean.lambda', 'mean.p', 'alpha0','alpha1',
            'beta0', 'beta1', 'beta2','beta3','beta4','beta5','beta6','beta7',
            'eta.lam', 'Sigma2', 'rho', 'N')
ni <- 400000 ; nb <- 100000 ; nt <- 300 ; na = 10000; nc = 3

# Initial values
Nst <- maxC
Nst[is.na(Nst)] <- 0
Nst = Nst + 1
modelInits <- function(){list(N = Nst, 
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
                              mean.p = rep(0.2, nspec), 
                              eta.lam = array(1, dim = c(548, nspec)),
                              lambda = array(1, dim = c(548, nspec)),
                              Omega = diag(nspec))}

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

 # Check the model convergence
library(MCMCvis)
sums = MCMCsummary(model_out)
summary(is.na(sums$Rhat))
check = filter(sums, Rhat>1.1) # which actual values have high Rhats
check

params.to.check = rownames(check)
MCMCvis::MCMCtrace(model_out,params=params.to.check, ISB=FALSE, Rhat = TRUE)
# ESS and Rhats look good

saveRDS(model_out, "Data_Clean/nmix_nimblemodel_final.RDS")
