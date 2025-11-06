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
multispp.data[multispp.data>0] <- 1 # Convert to binary presence-absence
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
covs_scaled = covs %>% 
  dplyr::select(-c(landscape_position, forest_position)) %>%
  mutate(tsf_sqrt = sqrt(tsf.point),
         bait_log = log(Mean_Intensity_400)) %>%
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
              bait_log = covs_scaled$bait_log,
              tsf=covs_scaled$tsf.point,
              tsf_sqrt = covs_scaled$tsf_sqrt,
              #tsf_spline = as.matrix(tsf_spline),
              #ag = covs_scaled$prop_ag_3km,
              #propnv = covs_scaled$prop_nv_3km,
              twi = covs_scaled$twi,
              #hydro = covs_scaled$dist_to_majorhydro,
              #road = covs_scaled$prop_filtered_roads_3km,
              rainfall = covs_scaled$rainfall,
              date= yday(dethist$Sites$Start)/365, # Date to integer and convert to fraction
              #landscape = covs_scaled$forest_position,
              R = Rmat, 
              df = df)

### RUN MODELS
source("Scripts/model_functions.R")
linear_model_out = run_linear_model(bdata)
linear_model_noint_out = run_linear_noint_model(bdata)
sqrt_model_out = run_sqrt_model(bdata)
poly_model_out = run_poly_model(bdata)
spline_model_out = run_spline_model(bdata)

saveRDS(linear_model_noint_out, "Data_Clean/nmix_nimblemodel_final_linear_noint.RDS")
saveRDS(linear_model_out, "Data_Clean/nmix_nimblemodel_final_linear.RDS")
saveRDS(sqrt_model_out, "Data_Clean/nmix_nimblemodel_final_sqrt.RDS")
saveRDS(poly_model_out, "Data_Clean/nmix_nimblemodel_final_poly.RDS")
saveRDS(spline_model_out, "Data_Clean/nmix_nimblemodel_final_spline.RDS")


linear_model_out = readRDS("Data_Clean/nmix_nimblemodel_final_linear.RDS")
linear_model_noint_out = readRDS("Data_Clean/nmix_nimblemodel_final_linear_noint.RDS")
sqrt_model_out = readRDS("Data_Clean/nmix_nimblemodel_final_sqrt.RDS")
poly_model_out = readRDS("Data_Clean/nmix_nimblemodel_final_poly.RDS")
spline_model_out = readRDS("Data_Clean/nmix_nimblemodel_final_spline.RDS")

# Check the wAICs
waics <- data.frame(Model = c("Linear","Linear - No Int", "Sqrt", "Poly", "Spline"),
                    wAIC = c(linear_model_out$WAIC$WAIC, linear_model_noint_out$WAIC$WAIC, sqrt_model_out$WAIC$WAIC, 
                             poly_model_out$WAIC$WAIC, spline_model_out$WAIC$WAIC),
                    pwAIC = c(linear_model_out$WAIC$pWAIC,  linear_model_noint_out$WAIC$pWAIC, sqrt_model_out$WAIC$pWAIC, 
                             poly_model_out$WAIC$pWAIC, spline_model_out$WAIC$pWAIC))
waics

# Check the model convergence
library(MCMCvis)
library(tidybayes)
model_out <- readRDS("Data_Clean/nmix_nimblemodel_final_linear_noint.RDS")

sums = MCMCsummary(model_out$samples)
summary(is.na(sums$Rhat))
check = filter(sums, Rhat>1.1) # which actual values have high Rhats
check

params.to.check = rownames(check)
MCMCvis::MCMCtrace(model_out,params=params.to.check, ISB=FALSE, Rhat = TRUE)
# ESS and Rhats look good

# Posterior Predictive checks
draws = model_out$samples %>% tidy_draws()
# Bayesian p value
# Interpretation: A p-value close to 0.5 indicates a good fit, 
# while values close to 0 or 1 suggest poor fit.
p_value_N <- mean(draws$chi2_N_sim_total > draws$chi2_N_obs_total)
p_value_N 

# Plot the results
chi2_values = draws %>% 
  dplyr::select(chi2_N_obs_total, chi2_N_sim_total) %>%
  ggplot() + 
  geom_point(aes(x = chi2_N_obs_total, y = chi2_N_sim_total), colour="#414487FF") + 
  geom_abline(intercept =0, slope = 1) +
  labs(x = expression("Observed " ~ chi^2), y = expression("Simulated " ~ chi^2)) +
  theme_bw() + xlim(2000, 4250) + ylim(2000, 4250) 
chi2_values

ggsave(plot = chi2_values, filename = "Outputs/model_ppcheck_v2.pdf", dpi = 300, width = 6, height=5, scale = 1)

# MAE
N_samples <- draws %>%  gather_variables() %>% filter(grepl("N", .variable)) %>% filter(!grepl("chi", .variable)) %>% 
  group_by(.variable) %>% summarise(meanN = mean(.value))
N_samples <- N_samples %>%  mutate(Site = parse_number(.variable),
         Species = as.numeric(str_extract_all(N_samples$.variable, "\\d+", simplify=TRUE)[,2]))

lambda_samples <- draws %>% gather_variables() %>% filter(grepl("lambda", .variable)) %>% filter(!grepl("mean", .variable)) %>% 
  group_by(.variable) %>% summarise(meanLambda = mean(.value)) %>% 
  mutate(Site = parse_number(.variable),
         Species = as.numeric(str_extract_all(N_samples$.variable, "\\d+", simplify=TRUE)[,2]))

abundances = left_join(N_samples, lambda_samples, by = c("Site", "Species"))

# Calculate Mean Absolute Error (MAE) 
abundances$error = abundances$meanLambda - abundances$meanN
mae <- mean(abs(abundances$error), na.rm = TRUE)  
mae

# Calculate Mean Relative Error (MRE) 
abundances$relative_error = abundances$error / abundances$meanLambda*100
mre <- mean(abs(abundances$relative_error), na.rm = TRUE) 
mre
