##++++++++++++++++++++++++++++++++++++##
## 2: Run the JAGS multispp occ model ##
##++++++++++++++++++++++++++++++++++++##

## Billy Geary
## March 2021

## Upper Warren Mammal Project

#### Step 1: Set Up ####
library(jagsUI)
library(dplyr)
library(AHMbook)
library(corrplot)
library(ggplot2)

# Read in data
dethist = readRDS("Data_Processing/camtrapR.dethist.multispp.12032021.RData")

# Select species
dethist$Species =  dethist$Species[c("Woylie", "Chuditch", "Koomal", "Quenda", "Feral Cat", "Vulpes")]

# Get info about the data
nsite = nrow(dethist$Sites) # number of sites
nrep = max(sapply(dethist$Species, ncol)) # maximum number of replicate surveys per season
nspec = length(dethist$Species)

# Convert det histories to a 3d array
multispp.data = dethist$Species
multispp.data = array(unlist(multispp.data), 
                      dim = c(dim(multispp.data[[1]]), length(multispp.data)),
                      dimnames = list(site = rownames(multispp.data[[1]]), rep = colnames(multispp.data[[1]]), sps = names(multispp.data)))

dim(multispp.data) == c(nsite, nrep, nspec) # Check the dimensions match

# Create the detection/nondetection (1/0) array
y <- multispp.data
class(y) <- 'numeric'

y[y > 1] <- 1  ## 'Y' replaced with 'y'
str(y)

# Check data for one species, and pull them out from 3D array
(tmp <- y[, , "Woylie"])

# Frequency distribution of number of surveys actually carried out per site
table(nsurveys <- apply(!is.na(y[,,1]), 1, sum))

# Observed number of occupied sites
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
# For the 'all NA' site, max returns -Inf with a warning
tmp[tmp == -Inf] <- NA         # Change -Inf to NA
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))

# Plot species 'occurrence frequency' distribution 
plot(sort(obs.occ), xlab = "Species number", ylab = "Number of sites with detections")

# Get observed number of species per site
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum))     # Compute and print sorted species counts

plot(table(C), xlab = "Observed number of species",
     ylab = "Number of Sites", frame = FALSE)
abline(v = mean(C, na.rm = TRUE), col = "blue", lwd = 3)


#### Step 2: Specify a Basic JAGS Model ####

# 11.6.1 Simplest community occupancy model: n-fold single species
#        occupancy model with species treated as fixed effects
# Collapse 3D detection/nondetection data to 2D detection frequencies
ysum <- apply(y, c(1,3), sum, na.rm = TRUE) # Collapse to detection frequency
ysum[NAsites,] <- NA                     # Have to NA out sites with NA data

# Bundle and summarize data set
str(win.data <- list(ysum = ysum, 
                     M = nrow(ysum), 
                     J = dethist$nsurvey,
                     nspec = dim(ysum)[2]) )

# Specify model in BUGS language
sink("Scripts/model5.txt")
cat("
model {
  # Priors
  for(k in 1:nspec){          # Loop over species
    psi[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  # Ecological model for latent occurrence z (process model)
  for(k in 1:nspec){          # Loop over species
    for (i in 1:M) {         # Loop over sites
      z[i,k] ~ dbern(psi[k])
    }
  }
  # Observation model for observed data y
  for(k in 1:nspec){          # Loop over species
    for (i in 1:M) {
      mup[i,k] <- z[i,k] * p[k]
      ysum[i,k] ~ dbin(mup[i,k], J[i])
    }
  }
  # Derived quantities
  for(k in 1:nspec){          # Loop over species
    Nocc.fs[k] <- sum(z[,k]) # Add up number of occupied sites among the 267
  }
  for (i in 1:M) {            # Loop over sites
    Nsite[i] <- sum(z[i,])   # Add up number of occurring species at each site
  }
}
",fill = TRUE)
sink()



#### Step 3: Run the Basic JAGs Models ####
# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as inits for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, psi = rep(0.4, nspec), p = rep(0.4, nspec))

# Parameters monitored
params <- c("psi", "p", "Nsite", "Nocc.fs")

# MCMC settings
ni <- 2500   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3

# Call JAGS from R (ART 2.1 min)
out5 <- jagsUI::jags(win.data, inits, params, "Scripts/model5.txt", n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

#### Step 4: Look at the results ####
# Get mean occupancy for each species
mean.occu = data.frame(Species = names(dethist$Species),
                       Param = "Psi",
                       Mean = out5$mean$psi,
                       Lower_2.5 = out5$q2.5$psi,
                       Upper_97.5 = out5$q97.5$psi)
mean.det= data.frame(Species = names(dethist$Species),
                     Param = "P",
                     Mean = out5$mean$p,
                     Lower_2.5 = out5$q2.5$p,
                     Upper_97.5 = out5$q97.5$p)
plot.data = rbind(mean.occu, mean.det)

ggplot(plot.data) + 
  geom_pointrange(aes(x = reorder(Species, -Mean), y = Mean, ymin=Lower_2.5, ymax = Upper_97.5, colour=Param), position = position_dodge(width = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



#### Step 5: Specify an Occupancy JSDM ####

# Create the detection/nondetection (1/0) array
y <- multispp.data
class(y) <- 'numeric'

y[y > 1] <- 1  ## 'Y' replaced with 'y'
str(y)

# Calculate sum of max counts across sites and obs. number of occ. sites
tmp <- apply(cc, c(1,3), max, na.rm = TRUE)
summax <- apply(tmp, 2, sum) # p-ignorant estimate of Ntotal
tmp <- apply(yy, c(1,3), max, na.rm = TRUE)
nobs <- apply(tmp, 2, sum) # p-ignorant estimate of sum(z)
sort(nobs) # Look at species ordered by nobs


# Determine sample sizes
(nsites <- dim(y)[1])
(nspec <- dim(y)[3])
table(nreps <- apply(!is.na(y[,,1]), 1, sum))

  
# Provide for different numbers of LVs: could try 2, 5, 10, 15
# Reccomendation is 2-5 is useful, but maybe n/2? Not sure how to decide... 
NLV <- c(2, 5, 10, 15) # Here we will just take 2

# Bundle data (incl. choice of number of LVs)
str(bdata <- list(y = aperm(y, c(1,3,2)), 
                  nlv = NLV[2], 
                  nsites = nsites,
                  nspec = nspec, 
                  nreps = nreps))


# Specify model in BUGS language
cat(file = "Scripts/JSDMocc.txt", "
model{
  # Community priors for occupancy
  mu.beta0 <- logit(mean.psi0) # Intercept
  mean.psi0 ~ dunif(0, 1)
  tau.beta0 <- pow(sd.beta0, -2)
  sd.beta0 ~ dunif(0, 2)

  # Community priors for detection
  mu.alpha0 <- logit(mean.p) # Intercept
  mean.p ~ dunif(0, 1)
  tau.alpha0 <- pow(sd.alpha0, -2)
  sd.alpha0 ~ dunif(0, 1)

  # Define species-specific random effects for all coefficients
  for (k in 1:nspec) {
    # Random species effects in the occupancy model
    beta0[k] ~ dnorm(mu.beta0, tau.beta0) # Intercepts

    # Random effects for detection
    alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0) # Intercepts
  }
  
  # Priors for latent variables: standard Normal rv
  for(i in 1:nsites) {
    for(l in 1:nlv){
      LV[i,l] ~ dnorm(0, 1)
    }
  }
  
  # Latent variable coefficients with constraints
  # Diagonal elements positive, upper diagonal equal to 0
  for(l in 1:(nlv-1)){
    for(l2 in (l+1):nlv){
      lv.coef[l,l2] <- 0
    }
  }
  
  ## Sign constraints on diagonal elements
  for(l in 1:nlv) {
    lv.coef[l,l] ~ dunif(0, 1)
  }
  
  # Lower diagonal free
  for(l in 2:nlv){
    for(l2 in 1:(l-1)){
      lv.coef[l,l2] ~ dunif(-1, 1)
    }
  }
  
  # Other elements free
  for(l in (nlv+1):nspec) {
    for(l2 in 1:nlv){
      lv.coef[l,l2] ~ dunif(-1, 1)
    }
  }
  
  # Define the multi-species occupancy model
  for (i in 1:nsites) {                          # Loop over sites
    for (k in 1:nspec) {                         # Loop over species
      
      # Probit link GLM for occupancy via auxiliary variable approach
      eta[i,k] <- beta0[k] + inprod(lv.coef[k,], LV[i, ])
      
      # Draw Gaussian auxiliary variable, with variance constrained to 1
      mu.psi[i,k] ~ dnorm(eta[i,k], 1/(1-sum(lv.coef[k,1:nlv]^2)))
      z[i,k] <- step(mu.psi[i,k])
      
      # Bernoulli GLM for detection
      for (j in 1:nreps[i]) {                    # Loop over 2 or 3 surveys
        logit(p[i,k,j]) <- alpha0[k]
        y[i,k,j] ~ dbern(p[i,k,j]*z[i,k])
      }
    }
  }
  
  # Derived quantities
  for(k in 1:nspec){
    mean.occ[k] <- mean(z[,k])
  }
  
}
")


#### Step 6: Run the Occupancy JSDM ####

# 'Blind' initial values
inits <- function() { # nsites = nsites; nspec = nspec;
  lv.coef <- matrix(0, nspec, bdata$nlv)
  lv.coef[1:bdata$nlv, 1:bdata$nlv] <- 0
  for(l in 1:bdata$nlv-1){ lv.coef[l, (l+1):bdata$nlv] <- NA}
  LV <- matrix(rnorm(bdata$nlv * nsites), nsites, bdata$nlv)
  lv.coef <- matrix(runif(bdata$nlv * nspec, -sqrt(1/(bdata$nlv+1)),
                          sqrt(1/(bdata$nlv+1))), nspec, bdata$nlv) * lv.coef
  mu.psi <- array(1, dim = c(nsites, nspec))         # yields psi = 0.84
  list(LV = LV, lv.coef = lv.coef , mu.psi = mu.psi)
}

# Parameters to be monitored
params <- c('mean.psi0', 'mu.beta0', 'sd.beta0', 'z', 'mean.occ',
            'mean.p', 'mu.alpha0', 'sd.alpha0', 'beta0',
            'alpha0', 'LV', 'lv.coef') # could add 'z'

# This model will probably suck to fit
## - Trade off between number of LVs and run time
## - May need to re-fit model using previous runs
## - Will take days to weeks to fit

# MCMC settings
# na <- 10000 ; ni <- 250000 ; nt <- 100 ; nb <- 150000 ; nc <- 3  # 3 days
na <- 1000 ; ni <- 250 ; nt <- 1 ; nb <- 150 ; nc <- 3 # ~~~ for testing, 35 mins

# Call JAGS (ART is long !), check convergence and summarize posteriors
parallel:::setDefaultClusterOptions(setup_strategy = "sequential") # do this so parallel works on MacOSX

out1 <- jagsUI::jags(bdata, inits, params, "Scripts/JSDMocc.txt", n.adapt = na,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Check convergance
conv = data.frame(subset(out1$summary, out1$summary[,8] > 1.1))

#### Step 7: Look at occJSDM results ####

## MEAN OCCUPANCY & DETECTABILITY
# Get mean occupancy for each species (need to sort this)
plot.data.occJSDM = data.frame(Species = names(dethist$Species),
                       Param = "Psi",
                       Mean = out1$mean$mean.occ,
                       Lower_2.5 = out1$q2.5$mean.occ,
                       Upper_97.5 = out1$q97.5$mean.occ)
ggplot(plot.data.occJSDM) + 
  geom_pointrange(aes(x = reorder(Species, -Mean), y = Mean, ymin=Lower_2.5, ymax = Upper_97.5, colour=Param), position = position_dodge(width = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Comparison with previous model
mean.occu$Model = "Basic"
plot.data.occJSDM$Model = "occuJSDM"
plot.data = rbind(mean.occu, plot.data.occJSDM)

ggplot(plot.data) + 
  geom_pointrange(aes(x = reorder(Species, -Mean), y = Mean, ymin=Lower_2.5, ymax = Upper_97.5, colour=Model), position = position_dodge(width = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## CORRELATION MATRIX FOR SPECIES
# Compute the posterior mean of the correlation matrix
R <- getLVcorrMat(lv.coef = out1$sims.list$lv.coef) # From AHMbook
colnames(R) <- rownames(R) <- dimnames(multispp.data)$sps

# Plot the correlation matrix (Fig. 8.16)
corrplot(R, type = "lower", diag = FALSE, mar = c(1,0.5,5,1), tl.col = 'black',
         tl.pos = 'ld', tl.srt = 45, xpd = TRUE, main = 'Residual correlations')

#### Step 8: Covariate relationships ####


#### Step 9: Map Predictions #### 



###################################
# Read in data
dethist = readRDS("Data_Processing/camtrapR.dethist.UpperWarren.multispp.21112021.RData")

# Get info about the data
nsite = nrow(dethist$Sites) # number of sites
nrep = max(sapply(dethist$Species, ncol)) # maximum number of replicate surveys per season
nspec = length(dethist$Species)

# Convert det histories to a 3d array
multispp.data = dethist$Species
# Select just the mammals for now

multispp.data = multispp.data[c('Chuditch','Koomal', "Numbat", "Ngwayir", "Tammar", "Phascogale",
                                'Quenda', 'Vulpes', 'Woylie')]
spp.names = names(multispp.data)
multispp.data = array(unlist(multispp.data), 
                      dim = c(dim(multispp.data[[1]]), length(multispp.data)),
                      dimnames = list(site = rownames(multispp.data[[1]]), rep = colnames(multispp.data[[1]]), sps = names(multispp.data)))

dim(multispp.data) == c(nsite, nrep, nspec) # Check the dimensions match


#### Step 5: Specify an Occupancy JSDM ####

# Create the detection/nondetection (1/0) array
y <- multispp.data
class(y) <- 'numeric'

y[y > 1] <- 1  ## 'Y' replaced with 'y'
str(y)

# Calculate sum of max counts across sites and obs. number of occ. sites
tmp <- apply(cc, c(1,3), max, na.rm = TRUE)
summax <- apply(tmp, 2, sum) # p-ignorant estimate of Ntotal
tmp <- apply(yy, c(1,3), max, na.rm = TRUE)
nobs <- apply(tmp, 2, sum) # p-ignorant estimate of sum(z)
sort(nobs) # Look at species ordered by nobs


# Determine sample sizes
(nsites <- dim(y)[1])
(nspec <- dim(y)[3])
table(nreps <- apply(!is.na(y[,,1]), 1, sum))

# Prep covariates
xocc <- as.matrix(dethist$Sites[,13:15])
str(xocc <- cbind(xocc, xocc[,2]^2))
xocc <- scale(xocc)                   # Scale column-wise
xocc[is.na(xocc)] <- 0


# Provide for different numbers of LVs: could try 2, 5, 10, 15
# Reccomendation is 2-5 is useful, but maybe n/2? Not sure how to decide... 
NLV <- c(2, 5, 10, 15) # Here we will just take 2

# Bundle data (incl. choice of number of LVs)
str(bdata <- list(y = aperm(y, c(1,3,2)), 
                  Xocc = xocc,
                  ncov.occ = ncol(xocc),
                  nlv = NLV[1], 
                  nsites = nsites,
                  nspec = nspec, 
                  nreps = nreps))


# Specify model in BUGS language
cat(file = "Scripts/JSDMocc.txt", "
model{
  # Community priors for occupancy
  mu.beta0 <- logit(mean.psi0) # Intercept
  mean.psi0 ~ dunif(0, 1)
  tau.beta0 <- pow(sd.beta0, -2)
  sd.beta0 ~ dunif(0, 2)
  for(v in 1:ncov.occ) { # Coefficients
    mu.beta[v] ~ dnorm(0, 0.2)
    tau.beta[v] <- pow(sd.beta[v], -2)
    sd.beta[v] ~ dunif(0, 2)
  }

  # Community priors for detection
  mu.alpha0 <- logit(mean.p) # Intercept
  mean.p ~ dunif(0, 1)
  tau.alpha0 <- pow(sd.alpha0, -2)
  sd.alpha0 ~ dunif(0, 1)

  # Define species-specific random effects for all coefficients
  for (k in 1:nspec) {
    # Random species effects in the occupancy model
    beta0[k] ~ dnorm(mu.beta0, tau.beta0) # Intercepts
    for(v in 1:ncov.occ) { # Coefficients
      beta[k, v] ~ dnorm(mu.beta[v], tau.beta[v])
    }

    # Random effects for detection
    alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0) # Intercepts
  }
  
  # Priors for latent variables: standard Normal rv
  for(i in 1:nsites) {
    for(l in 1:nlv){
      LV[i,l] ~ dnorm(0, 1)
    }
  }
  
  # Latent variable coefficients with constraints
  # Diagonal elements positive, upper diagonal equal to 0
  for(l in 1:(nlv-1)){
    for(l2 in (l+1):nlv){
      lv.coef[l,l2] <- 0
    }
  }
  
  ## Sign constraints on diagonal elements
  for(l in 1:nlv) {
    lv.coef[l,l] ~ dunif(0, 1)
  }
  
  # Lower diagonal free
  for(l in 2:nlv){
    for(l2 in 1:(l-1)){
      lv.coef[l,l2] ~ dunif(-1, 1)
    }
  }
  
  # Other elements free
  for(l in (nlv+1):nspec) {
    for(l2 in 1:nlv){
      lv.coef[l,l2] ~ dunif(-1, 1)
    }
  }
  
  # Define the multi-species occupancy model
  for (i in 1:nsites) {                          # Loop over sites
    for (k in 1:nspec) {                         # Loop over species
      
      # Probit link GLM for occupancy via auxiliary variable approach
      eta[i,k] <- beta0[k] + inprod(beta[k, ], Xocc[i, ]) +
      inprod(lv.coef[k,], LV[i, ])
      
      # Draw Gaussian auxiliary variable, with variance constrained to 1
      mu.psi[i,k] ~ dnorm(eta[i,k], 1/(1-sum(lv.coef[k,1:nlv]^2)))
      z[i,k] <- step(mu.psi[i,k])
      
      # Bernoulli GLM for detection
      for (j in 1:nreps[i]) {                    # Loop over 2 or 3 surveys
        logit(p[i,k,j]) <- alpha0[k]
        y[i,k,j] ~ dbern(p[i,k,j]*z[i,k])
      }
    }
  }
  
  # Derived quantities
  for(k in 1:nspec){
    mean.occ[k] <- mean(z[,k])
  }
  
}
")


#### Step 6: Run the Occupancy JSDM ####

# 'Blind' initial values
inits <- function() { # nsites = nsites; nspec = nspec;
  lv.coef <- matrix(0, nspec, bdata$nlv)
  lv.coef[1:bdata$nlv, 1:bdata$nlv] <- 0
  for(l in 1:bdata$nlv-1){ lv.coef[l, (l+1):bdata$nlv] <- NA}
  LV <- matrix(rnorm(bdata$nlv * nsites), nsites, bdata$nlv)
  lv.coef <- matrix(runif(bdata$nlv * nspec, -sqrt(1/(bdata$nlv+1)),
                          sqrt(1/(bdata$nlv+1))), nspec, bdata$nlv) * lv.coef
  mu.psi <- array(1, dim = c(nsites, nspec))         # yields psi = 0.84
  list(LV = LV, lv.coef = lv.coef , mu.psi = mu.psi)
}

# Parameters to be monitored
params <- c('mean.psi0', 'mu.beta0', 'sd.beta0', 'z', 'mean.occ',
            'mean.p', 'mu.alpha0', 'sd.alpha0', 'beta0', 'beta',
            'alpha0', 'LV', 'lv.coef') # could add 'z'

# This model will probably suck to fit
## - Trade off between number of LVs and run time
## - May need to re-fit model using previous runs
## - Will take days to weeks to fit

# MCMC settings
na <- 100 ; ni <- 2500 ; nt <- 1 ; nb <- 1500 ; nc <- 3  # Took 2952.331 min for mammal only dataset
# na <- 1000 ; ni <- 25000 ; nt <- 10 ; nb <- 15000 ; nc <- 3 # ~~~ for testing, 171 minutes for mammals only dataset

# Call JAGS (ART is long !), check convergence and summarize posteriors
parallel:::setDefaultClusterOptions(setup_strategy = "sequential") # do this so parallel works on MacOSX

out1 <- jagsUI::jags(bdata, inits, params, "Scripts/JSDMocc.txt", n.adapt = na,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

#### PLOTS ####
R <- getLVcorrMat(lv.coef = out1$sims.list$lv.coef) # From AHMbook
colnames(R) <- rownames(R) <- dimnames(multispp.data)$sps

# Plot the correlation matrix (Fig. 8.16)
corrplot(R, type = "lower", diag = FALSE, mar = c(1,0.5,5,1), tl.col = 'black',
         tl.pos = 'ld', tl.srt = 45, xpd = TRUE, main = 'Residual correlations')


## Plot model coefficients and 95% credible intervals for each species
elev.cov = data.frame(Species = spp.names,
                      Covariate = "Elevation",
                        Mean = out1$mean$beta[,1],
                        Lower_2.5 = out1$q2.5$beta[,1],
                        Upper_97.5 = out1$q97.5$beta[,1])
nv.cov = data.frame(Species = spp.names,
                      Covariate = "Distance to Edge",
                      Mean = out1$mean$beta[,2],
                      Lower_2.5 = out1$q2.5$beta[,2],
                      Upper_97.5 = out1$q97.5$beta[,2])
tsf.cov = data.frame(Species = spp.names,
                      Covariate = "Time Since Fire",
                      Mean = out1$mean$beta[,3],
                      Lower_2.5 = out1$q2.5$beta[,3],
                      Upper_97.5 = out1$q97.5$beta[,3])
tsf2.cov = data.frame(Species = spp.names,
                     Covariate = "Time Since Fire^2",
                     Mean = out1$mean$beta[,4],
                     Lower_2.5 = out1$q2.5$beta[,4],
                     Upper_97.5 = out1$q97.5$beta[,4])

coef.data = rbind(elev.cov, nv.cov, tsf.cov, tsf2.cov)

ggplot(coef.data) + 
  geom_pointrange(aes(x = Species, y = Mean, ymin = Lower_2.5, ymax = Upper_97.5)) + 
  geom_hline(yintercept = 0, linetype='dashed') + ylab("Coeffcient and 95% CIs") +
  coord_flip() + facet_wrap(~Covariate, ncol=3, scales='free_x') + theme_bw()

### Look at creating rasters
xocc <- as.matrix(dethist$Sites[,13:15])
str(xocc <- cbind(xocc, xocc^2))

mean.ele = mean(xocc[,1], na.rm=TRUE)
sd.ele = sd(xocc[,1], na.rm=TRUE)
ele = raster("Data_Clean/covariate_StudyRegion_DEM.tif") 
ele.df = data.frame(rasterToPoints(ele))
names(ele.df) <- c("X","Y","elev")
ele.df$elev = (ele.df$elev - mean.ele)/sd.ele

mean.nv = mean(xocc[,3], na.rm=TRUE)
sd.nv = sd(xocc[,3], na.rm=TRUE)
nv = raster("Data_Clean/covariate_StudyRegion_NVdist.tif")
nv.df = data.frame(rasterToPoints(nv))
names(nv.df) <- c("X","Y","nv.dist")
nv.df$nv.dist = (nv.df$nv.dist - mean.nv)/sd.nv

mean.tsf = mean(xocc[,2], na.rm=TRUE)
sd.tsf = sd(xocc[,2], na.rm=TRUE)
tsf = raster("Data_Clean/covariate_StudyRegion_tsf.start.Oct16.tif")
tsf.df = data.frame(rasterToPoints(tsf))
names(tsf.df) <- c("X","Y","tsf.point")
tsf.df$tsf.point = (tsf.df$tsf.point - mean.tsf)/sd.tsf

pred.nd = left_join(ele.df, nv.df)
pred.nd = left_join(pred.nd, tsf.df)

sp.id = which(spp.names=="Ngwayir")
pred.nd$logit.psi.sp1 = with(out1$mean, 
                       beta0[sp.id] + beta[sp.id,1]*pred.nd$elev + 
                         beta[sp.id,2]*pred.nd$nv.dist +
                         beta[sp.id,3]*pred.nd$tsf.point + beta[sp.id,4]*pred.nd$tsf.point^2)
pred.nd$psi.sp1 = 1 / (1 + exp(-pred.nd$logit.psi.sp1))
pred.ras = rasterFromXYZ(data.frame(pred.nd$X, pred.nd$Y, pred.nd$psi.sp1))
plot(pred.ras)
plot(trap.points, add=TRUE, color='black')

# Mean relationships with covariates (spaghetti plots)
o.ele = seq(from = min(dethist$Sites$elev), to = max(dethist$Sites$elev), length.out = 500)
mean.ele = mean(dethist$Sites$elev)
sd.ele = sd(dethist$Sites$elev)
ele.pred = (o.ele - mean.ele)/sd.ele

o.nv = seq(min(dethist$Sites$nv.dist), max(dethist$Sites$nv.dist), length.out= 500)
mean.nv = mean(dethist$Sites$nv.dist)
sd.nv = sd(dethist$Sites$nv.dist)
nv.pred = (o.nv - mean.nv)/sd.nv

o.tsf = seq(min(dethist$Sites$tsf.point, na.rm=TRUE), max(dethist$Sites$tsf.point, na.rm=TRUE), length.out=500)
mean.tsf = mean(dethist$Sites$tsf.point, na.rm=TRUE)
sd.tsf = sd(dethist$Sites$tsf.point, na.rm=TRUE)
tsf.pred = (o.tsf - mean.tsf)/sd.tsf

tmp <- out1$sims.list
nsamp <- length(tmp[[1]])
predS <- array(NA, dim=c(500,nspec, 3))
psi.coef <- data.frame(out1$mean$beta)
names(psi.coef) <- c('betapsi1', 'betapsi2', 'betapsi3', 'betapsi4')
psi.coef$betapsi0 <- out1$mean$beta0

for(i in 1:nspec){
  predS[,i,1] <- plogis(psi.coef$betapsi0[i] + psi.coef$betapsi1[i]*ele.pred)
  predS[,i,2] <- plogis(psi.coef$betapsi0[i] + psi.coef$betapsi3[i]*nv.pred)
  predS[,i,3] <- plogis(psi.coef$betapsi0[i] + psi.coef$betapsi2[i]*tsf.pred + psi.coef$betapsi4[i]*tsf.pred^2)
}

plot(o.ele, predS[,1,3], lwd=3, type='l', lty=1, frame=F, ylim=c(0,1), xlab="Elevation", ylab="Occupancy Probability")
for(i in 1:31){
  lines(o.ele, predS[,i,1], col=i, lwd=3)
}

plot(o.nv, predS[,2,3], lwd=3, type='l', lty=1, frame=F, ylim=c(0,1), xlab="Distance to Edge", ylab="Occupancy Probability")
for(i in 1:31){
  lines(o.nv, predS[,i,2], col=i, lwd=3)
}

plot(o.tsf, predS[,3,3], lwd=3, type='l', lty=1, frame=F, ylim=c(0,1), xlab="Time Since Fire", ylab="Occupancy Probability")
for(i in 1:31){
  lines(o.tsf, predS[,i,3], col=i, lwd=3)
}

sp.id = which(spp.names=="Vulpes")
plot(o.tsf, predS[,sp.id,3], type='l', col='red', lwd=3)




