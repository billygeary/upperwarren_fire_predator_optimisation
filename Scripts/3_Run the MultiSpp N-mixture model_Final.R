##++++++++++++++++++++++++++++++++++++##
## 3: Run the multispp Nmixture model ##
##++++++++++++++++++++++++++++++++++++##
## Billy Geary
## March 2021

## Upper Warren Mammal Project

#### Step 1: Set Up ####
library(dplyr)
library(AHMbook)
library(corrplot)
library(ggplot2)
library(nimble)
library(parallel)
library(tidybayes)
library(ggmcmc)
library(coda)

# Read in data
dethist = readRDS("Data_Processing/camtrapR.counthist.UpperWarren.multispp.RData")

# Select species we want to include in the Multi Spp Occ model
names(dethist$Species)
dethist$Species =  dethist$Species[c("Woylie", "Chuditch", "Koomal", "Quenda", "Roo",
                                     "Vulpes", "Tammar", "Numbat", "Western Brush Wallaby", "Dunnart")]
spp = names(dethist$Species)

# Collapse nightly detection history to weekly occasions
# truncdethist = dethist
# for (s in spp){
#   d = truncdethist$Species[[s]]
#   trunc.d=data.frame(o1 = rowSums(d[,1:7]),
#                      o2 = rowSums(d[,8:14]),
#                      o3 = rowSums(d[,15:21]),
#                      o4 = rowSums(d[,22:28]))
#   truncdethist$Species[[s]] <- trunc.d
#   }
# 
# dethist = truncdethist

# Check for colinearity in the covariates
covs = dethist$Sites %>% 
  dplyr::select(X,Y,10:31) %>% dplyr::select(-c(Study.AreaName, Lag))
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

#### Step 2: Setup Model Inputs ####
y <- multispp.data
#y <- y[,1:28,] # Select only first 28 days so even sampling
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
                 east = covs_scaled$X,
                 north = covs_scaled$Y,
                 propsev = covs_scaled$PropSevere500,
                 bait = covs_scaled$Mean_Intensity_400,
                 tsf=covs_scaled$tsf.point,
                 agdist = covs_scaled$dist_to_ag,
                 propnv = covs_scaled$prop_nv_3km,
                 twi = covs_scaled$twi,
                 hydro = covs_scaled$dist_to_majorhydro,
                 road = covs_scaled$prop_filtered_roads_3km,
                 date=covs_scaled$date,
                 landscape = covs_scaled$forest_position,
                 R = Rmat, 
                 df = df)

#### Step 3: Fit Model ####
source("Scripts/3b_ModelFunctions_Final.R")

# Initial values
Nst <- maxC + 1
Nst[is.na(Nst)] <- 1
modelInits <- list(N = Nst, 
                   mean.p = rep(0.2, nspec), 
                   Omega = diag(nspec))

#ni <- 300000 ; nb <- 100000 ; nt <- 200 ; na = 10000; nc = 3
# ni <- 200 ; nb <- 100 ; nt <- 1 ; nc=3 ; na = 100
 ni <- 2000 ; nb <- 1000 ; nt <- 1 ; nc=3 ; na = 1000

# Fit Model
start = Sys.time()
model_output = run_nmix_fire_jagscode(data = bdata, nt = nt, ni = ni, nb = nb, nc=nc, na = na)
end = Sys.time()
(time.taken = end - start)


# Save model output
saveRDS(model_output, "Data_Processing/nmix_modeloutput.RData") 

##### Step 4: Assess convergence ####
model_output = readRDS("Data_Processing/nmix_modeloutput_19042022.RData")

# Rhat values
rhats = model_output$Rhat
lapply(rhats, FUN = function(x){which(x >1.1)}) # Which parameters are >1.1

#### Step 8: Predicted Community Measures ####

# Species Richness
richness = 
  model.samples %>% 
  tidy_draws() %>%
  gather_variables() %>%
  filter(.variable %in% paste0("Nsite[",1:nsites,"]")) %>%
  mutate(Site = readr::parse_number(.variable)) %>%
  group_by(Site, .variable) %>% 
  summarise(Rich_Mean = mean(.value), Rich_Lower_2.5 = quantile(.value, probs = 0.025), Rich_Upper_97.5 = quantile(.value, probs = 0.975))

# Geometric Mean Abundance
site.spp = expand.grid(Sites = 1:nsites, Species = 1:nspec)
gma = 
  model.samples %>%
  tidy_draws %>%
  gather_variables() %>%
  filter(.variable %in% paste0("N[",site.spp$Sites,", ", site.spp$Species,"]"))
gma$Site = readr::parse_number(gma$.variable)
gma$Species = readr::parse_number(gsub(".*\\,", "", gma$.variable))
gma = gma %>%
  group_by(Species, Site) %>% summarise(MeanAbundance = mean(.value))

gma2 = gma %>%
  pivot_wider(id_cols = Site, names_from = Species, values_from = MeanAbundance)

gma.metrics = gma %>%
  group_by(Site) %>% summarise(GMA = exp(mean(log(MeanAbundance))), 
                               SumSppAbundance = sum(MeanAbundance))
  
gma2 = left_join(gma2, gma.metrics, by = "Site")

community.measures = left_join(richness, gma2, by="Site")
