####### MASTER DATA TIDYING - MONTE CARLO  ########
###### 04-25-2024 #####
### This script runs the data preparations in a Monte Carlo fashion ##
library(tidyverse)

set.seed(123456789) 

###### -------------- BASE DATA PREPARATION ------------------------- ###################
###### Define parameters used for particle property estimates ####
R.ave.water.marine <- 0.77 # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
R.ave.water.freshwater <- 0.67
R.ave.sediment.marine <- 0.75
R.ave.sediment.freshwater <- 0.70
beta_log10_body_length <- 0.9341
body_length_intercept <- 1.1200

### --- LOAD FUNCTIONS --- ###
# each script is now a function. Source files to load functions #
source("scripts/monte carlo/RDAmaker_functions.R") #get ToMEx1.0 fxn
source("scripts/monte carlo/ToMEx2.0_Data_Tidying_functions.R") 

####### ---- RUN FUNCTIONS ---- ###
### Static with base parameters to ensure everything is working ###
## generate ToMEx 1.0 dataset ###
aoc_setup <- ToMEx1.0fxn(R.ave.water.marine = R.ave.water.marine,
                         R.ave.water.freshwater = R.ave.water.freshwater,
                         R.ave.sediment.marine = R.ave.sediment.marine,
                         R.ave.sediment.freshwater = R.ave.sediment.freshwater,
                         beta_log10_body_length = beta_log10_body_length,
                         body_length_intercept = body_length_intercept) 

## generate ToMEx 2.0 dataset ###
tomex2.0_aoc_z_final <- ToMEx2.0fxn(aoc_setup = aoc_setup,
                                    R.ave.water.marine = R.ave.water.marine,
                                    R.ave.water.freshwater = R.ave.water.freshwater,
                                    R.ave.sediment.marine = R.ave.sediment.marine,
                                    R.ave.sediment.freshwater = R.ave.sediment.freshwater,
                                    beta_log10_body_length = beta_log10_body_length,
                                    body_length_intercept = body_length_intercept)

### Check if these files are identical to the legacy ones produced in RDAmaker.R and the data tidying script 
aoc_setup_OG <- read_rds("data/input/aoc_setup.RDS")
tomex2.0_aoc_z_final_OG <- read_rds("data/input/aoc_z_tomex2.RDS")

### these must return TRUE
# identical(aoc_setup, aoc_setup_OG)
# 
# tomex2.0_aoc_z_final_noRAve <- tomex2.0_aoc_z_final %>% ungroup() %>% select(-R.ave) 
#                                 
# identical(tomex2.0_aoc_z_final_noRAve, #rowid is assigned randomly 
#           tomex2.0_aoc_z_final_OG)
# 
# differences <- tomex2.0_aoc_z_final_noRAve  != tomex2.0_aoc_z_final_OG
# 
# # Assuming 'differences' is your logical matrix with TRUE/FALSE/NA values
# # Set NAs to FALSE (if NA means no difference for your case)
# differences[is.na(differences)] <- FALSE
# 
# # Get rows and columns indices where differences are TRUE
# diff_indices <- which(differences, arr.ind = TRUE)
# 
# # Get unique rows and columns involved in differences
# diff_rows <- unique(diff_indices[, "row"])
# diff_cols <- unique(diff_indices[, "col"])
# 
# # If differences are found, print them out
# if (length(diff_rows) > 0) {
#   cat("Differences found in the following rows and columns:\n")
#   for (i in seq_along(diff_rows)) {
#     for (j in seq_along(diff_cols)) {
#       if (differences[diff_rows[i], diff_cols[j]]) {  # Check if TRUE
#         # Using column names instead of indices
#         cat(sprintf("Difference at Row: %d, Column: %s\n",
#                     diff_rows[i], colnames(differences)[diff_cols[j]]))
#       }
#     }
#   }
# } else {
#   cat("No differences found.\n")
# }

#differences are due to different R.Ave use based on environment - need to fix this up in script

###### -------------- ALIGNMENTS ------------------------- ###################
##### ------ STATIC test ---- Monte Carlo is aftert this ------ ####
### Define parameters
x1M_set <- 1 #um lower size for all alignments
x1D_set <- 1 #um lower size for all alignments
x2D_set <- 5000 #um
upper.tissue.trans.size.um <- 83 #10 #um #set size for x2M
alpha = 2.07 #table s4 for marine surface water. length
a.sa = 1.5 #marine surface area power law
a.v = 1.48 #a_V for marine surface water volume
a.m = 1.32 # upper limit fora_m for mass for marine surface water in table S4 
a.ssa = 1.98 # A_SSA for marine surface water
#define additional parameters for calculations based on averages in the environment
R.ave.water.marine <- 0.77 # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
R.ave.water.freshwater <- 0.67
R.ave.sediment.marine <- 0.75
R.ave.sediment.freshwater <- 0.70
p.ave = 1.10 #average density in marine surface water

### perform alignments
aoc_aligned_test <- tomex2.0_aoc_z_final  %>% 
  ungroup() %>% 
  rename(environment = env_f) %>% 
  ### First filter the data ####
## First filter data with global filters
  filter(!environment %in% c("Terrestrial", "Not Reported"),
         Group != "Bacterium",
         Group != "Plant",
         effect.metric != "HONEC",
         tier_zero_tech_f == "Red Criteria Passed",
         tier_zero_risk_f == "Red Criteria Passed", #All thresholds must pass technical and risk red criteria
         risk.13 != 0 #Drop studies that received a score of 0 for endpoints criteria (this also drops studies that have not yet been scored) - KEEP THIS AFTER THE RED CRITERIA FILTERS  
  ) %>% 
  #Remove 26C temperature treatment data from Jaimukar et al. 2018
  filter(!(article == 42 & media.temp == 26)) %>% 
  mutate(max.size.ingest.um = 1000 * max.size.ingest.mm) %>%  #makes it less confusing below
  # calculate ERM for each species
  #### TISSUE TRANSLOCATION ####
# define upper size length for Translocation 
#set to 83um for upper limit or max size ingest, whichever is smaller
mutate(x2M_trans = case_when(is.na(max.size.ingest.um) ~ upper.tissue.trans.size.um, 
                             max.size.ingest.um  < upper.tissue.trans.size.um ~  max.size.ingest.um,
                             max.size.ingest.um  > upper.tissue.trans.size.um ~ upper.tissue.trans.size.um)) %>% 
  
  # calculate effect threshold for particles
  mutate(EC_mono_p.particles.mL_trans = dose.particles.mL.master) %>% 
  mutate(mu.p.mono = 1) %>% #mu_x_mono is always 1 for particles to particles
  mutate(mu.p.poly_trans = mux_polyfnx_generalizable(a.x = alpha, #alpha for particles
                                       x_UL= x2M_trans, #upper ingestible size limit (width of particle)
                                       x_LL = x1M_set)) %>% 
  # polydisperse effect threshold for particles
  mutate(EC_poly_p.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.p.mono)/mu.p.poly_trans) %>% 
  #calculate CF_bio for all conversions
  mutate(CF_bio_trans = CFfnx(x1M = x1M_set,#lower size bin
                              x2M = x2M_trans, #upper translocatable
                              x1D = x1D_set, #default
                              x2D = x2D_set,  #default
                              a = alpha)) %>%  
  ## Calculate environmentally relevant effect threshold for particles
  mutate(EC_env_p.particles.mL_trans = EC_poly_p.particles.mL_trans * CF_bio_trans) %>%  #aligned particle effect concentraiton (1-5000 um)
  
  #### Surface area ERM ####
##--- environmental calculations ---###
#calculate lower translocatable surface area
mutate(x_LL_sa_trans = SAfnx(a = 0.5 * x1D_set, #length
                             b = 0.5 * x1D_set, #0.5 * R.ave * x1D_set, #width
                             c = 0.5 * x1D_set  #0.5 * R.ave * 0.67 * x1D_set #height
)) %>%  
  #calculate upper translocatable surface area
  mutate(x_UL_sa_trans = SAfnx(a = 0.5 * x2M_trans, 
                               b = 0.5 * x2M_trans, #width #0.5 * R.ave * x2M, 
                               c = 0.5 * x2M_trans #heigth #0.5 * R.ave * 0.67 * x2M
  )) %>%  
  #calculate mu_x_poly (env) for surface area
  mutate(mu.sa.poly_trans = mux_polyfnx_generalizable(a.sa, x_UL_sa_trans, x_LL_sa_trans)) %>% 
  
  ##--- laboratory calculations ---###
  ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
  #(note that if mixed particles were used, a different equation must be used)
  mutate(mu.sa.mono = case_when(
    polydispersity == "monodisperse" ~ particle.surface.area.um2, # use reported surface area in monodisperse
    polydispersity == "polydisperse" ~  mux_polyfnx_generalizable(a.x = a.sa, 
                                                    x_LL = particle.surface.area.um2.min,
                                                    x_UL = particle.surface.area.um2.max))) %>% 
  
  #calculate polydisperse effect concentration for surface area (particles/mL)
  mutate(EC_poly_sa.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.sa.mono)/mu.sa.poly_trans) %>%  
  #calculate environmentally realistic effect threshold
  mutate(EC_env_sa.particles.mL_trans = EC_poly_sa.particles.mL_trans * CF_bio_trans) %>% 
  
  ##### FOOD DILUTION ####
# define upper size length for ingestion 
mutate(x2M_ingest = case_when(is.na(max.size.ingest.um) ~ x2D_set, 
                              max.size.ingest.um < x2D_set ~ max.size.ingest.um,
                              max.size.ingest.um > x2D_set ~ x2D_set
)) %>%  #set to 5,000 as upper limit or max size ingest, whichever is smaller
  # calculate effect threshold for particles
  mutate(EC_mono_p.particles.mL_ingest = dose.particles.mL.master) %>% 
  mutate(mu.p.mono = 1) %>% #mu_x_mono is always 1 for particles to particles
  mutate(mu.p.poly_ingest = mux_polyfnx_generalizable(a.x = alpha, #alpha for particles
                                        x_UL= x2M_ingest, #upper ingestible size limit
                                        x_LL = x1M_set)) %>% 
  # polydisperse effect threshold for particles
  mutate(EC_poly_p.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.p.mono)/mu.p.poly_ingest) %>% 
  #calculate CF_bio for all conversions
  mutate(CF_bio_ingest = CFfnx(x1M = x1M_set,#lower size bin
                               x2M = x2M_ingest, #upper ingestible length
                               x1D = x1D_set, #default
                               x2D = x2D_set,  #default upper size range
                               a = alpha)) %>%  
  ## Calculate environmentally relevant effect threshold for particles
  mutate(EC_env_p.particles.mL_ingest = EC_poly_p.particles.mL_ingest * CF_bio_ingest) %>%  #aligned particle effect concentraiton (1-5000 um)
  
  
  #### volume ERM ####
##--- environmental calculations ---###
#calculate lower ingestible volume 
mutate(x_LL_v_ingest = volumefnx_poly(length = x1D_set,
                                      width = x1D_set)) %>% 
  #calculate maximum ingestible volume 
  mutate(x_UL_v_ingest = volumefnx_poly(length = x2M_ingest, # length-limited
                                        #x2D_set, #upper definiton (accouunts for fibers) CONSERVATIVE
                                        width = x2M_ingest)) %>% #ingestion-limited
  # calculate mu.v.poly
  mutate(mu.v.poly_ingest = mux_polyfnx_generalizable(a.v, x_UL_v_ingest, x_LL_v_ingest)) %>% 
  ##--- laboratory calculations ---###
  ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
  #(note that if mixed particles were used, a different equation must be used)
  mutate(mu.v.mono = case_when(
    polydispersity == "monodisperse" ~ particle.volume.um3, # use reported volume in monodisperse
    polydispersity == "polydisperse" ~ mux_polyfnx_generalizable(a.x = a.v, 
                                                   x_LL = particle.volume.um3.min,
                                                   x_UL = particle.volume.um3.max))) %>% 
  
  #calculate polydisperse effect concentration for volume (particles/mL)
  mutate(EC_poly_v.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.v.mono)/mu.v.poly_ingest) %>%  
  #calculate environmentally realistic effect threshold
  mutate(EC_env_v.particles.mL_ingest = EC_poly_v.particles.mL_ingest * CF_bio_ingest) %>% 
  
  ###### CLEANUP #####
mutate(particles.mL.ox.stress = EC_env_sa.particles.mL_trans,
       particles.mL.food.dilution = EC_env_v.particles.mL_ingest) %>% 
  select(-rowid) %>% 
  mutate(unique_id = row_number()) %>% 
 # mutate(unique_id = digest::digest(paste(across(everything()), collapse = "-"), algo = "md5")) %>%
  ungroup()

#### CHECKS TO MAKE SURE IT WORKS AS PLANNED ####
# load aligned dataset from ToMEx2.0_EcoToxRisk project (data/output/aoc_final.RDS)
aoc_aligned_ref <-read_rds("scripts/monte carlo/ref data/aoc_final.RDS") %>% select(-c(rowid, alpha, a.sa, a.v, a.m, a.ssa, R.ave, p.ave))

identical(aoc_aligned_test, 
          aoc_aligned_ref)
### Cool! ##



####### MONTE CARLO TIME, BABY!!! ####
##### STEP 1: DERIVE VALUES TO RUN PROBABILISTICALLY #####
source("scripts/monte carlo/ssd_functions.R")
nboot = 10 #ssd bootstraps
n_sim <- 100 #monte carlo simulations

# particle properties
R.ave.water.marine <- 0.77 # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
R.ave.water.marine.sd = 0.29 #Tablse s3. Marine surface water
R.ave.water.marine_samples <- rnorm(n_sim, mean = R.ave.water.marine, sd = R.ave.water.marine.sd)

R.ave.water.freshwater <- 0.67 # average length to width ratio of microplastics in freshwater environment (Kooi et al. 2021)
R.ave.water.freshwater.sd = 0.28 #Tablse s3. freshwater surface water
R.ave.water.freshwater_samples <- rnorm(n_sim, mean = R.ave.water.freshwater, sd = R.ave.water.freshwater.sd)

R.ave.sediment.marine <- 0.75 # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
R.ave.sediment.marine.sd = 0.30 #Tablse s3. Marine surface sediment
R.ave.sediment.marine_samples <- rnorm(n_sim, mean = R.ave.sediment.marine, sd = R.ave.sediment.marine.sd)

R.ave.sediment.freshwater <- 0.70 # average length to width ratio of microplastics in freshwater environment (Kooi et al. 2021)
R.ave.sediment.freshwater.sd = 0.33 #Tablse s3. freshwater surface sediment
R.ave.sediment.freshwater_samples <- rnorm(n_sim, mean = R.ave.sediment.freshwater, sd = R.ave.sediment.freshwater.sd)

#### MARINE ####
p.ave.marine = 1.10 #average density in marine surface water
p.ave.sd.marine = 0.14 #Tablse s3. Marine surfac water
p.ave.marine_samples <- rnorm(n_sim, mean = p.ave.marine, sd = p.ave.sd.marine)

# Alignment properties
alpha.marine = 2.07 #table s4 for marine surface water. length
alpha.sd.marine = 0.03 #table s4 for marine surface water. lengthj
alpha.marine_samples <- rnorm(n_sim, mean = alpha.marine, sd = alpha.sd.marine)

# define parameters for power law coefficients
a.sa.marine = 1.50 #marine surface area power law
a.sa.sd.marine = 0.009 #marine surface water surface area power law - table s4
a.sa.marine_samples <- rnorm(n_sim, mean = a.sa.marine, sd = a.sa.sd.marine)

a.v.marine = 1.48 #a_V for marine surface water volume
a.v.sd.marine = 0.063
a.v.marine_samples <- rnorm(n_sim, mean = a.v.marine, sd = a.v.sd.marine)

a.m.marine = 1.32 # upper limit fora_m for mass for marine surface water in table S4 
a.m.sd.marine = 0.009
a.m.marine_samples <- rnorm(n_sim, mean = a.m.marine, sd = a.m.sd.marine)

a.ssa.marine = 1.98 # A_SSA for marine surface water
a.ssa.sd.marine = 0.297
a.ssa.marine_samples <- rnorm(n_sim, mean = a.ssa.marine, sd = a.ssa.sd.marine)

#### FRESHWATER ####
p.ave.freshwater = 1.04 #average density in freshwater surface water
p.ave.sd.freshwater = 0.12 #Tablse s3. freshwater surfac water
p.ave.freshwater_samples <- rnorm(n_sim, mean = p.ave.freshwater, sd = p.ave.sd.freshwater)

# Alignment properties
alpha.freshwater = 2.64 #table s4 for freshwater surface water. length
alpha.sd.freshwater = 0.01 #table s4 for freshwater surface water. lengthj
alpha.freshwater_samples <- rnorm(n_sim, mean = alpha.freshwater, sd = alpha.sd.freshwater)

# define parameters for power law coefficients
a.sa.freshwater = 2.00 #freshwater surface area power law
a.sa.sd.freshwater = 0.065 #freshwater surface water surface area power law - table s4
a.sa.freshwater_samples <- rnorm(n_sim, mean = a.sa.freshwater, sd = a.sa.sd.freshwater)

a.v.freshwater = 1.68 #a_V for freshwater surface water volume
a.v.sd.freshwater = 0.081
a.v.freshwater_samples <- rnorm(n_sim, mean = a.v.freshwater, sd = a.v.sd.freshwater)

a.m.freshwater = 1.65 # upper limit fora_m for mass for freshwater surface water in table S4 
a.m.sd.freshwater = 0.071
a.m.freshwater_samples <- rnorm(n_sim, mean = a.m.freshwater, sd = a.m.sd.freshwater)

a.ssa.freshwater = 2.71 # A_SSA for freshwater surface water
a.ssa.sd.freshwater = 0.009
a.ssa.freshwater_samples <- rnorm(n_sim, mean = a.ssa.freshwater, sd = a.ssa.sd.freshwater)

# Coefficients
beta_log10_body_length <- 0.9341
body_length_intercept <- 1.1200
# Standard errors for these coefficients
se_beta_log10_body_length <- 0.1376  # SE generated from Jams et al R code (not rpoerted in paper or SI!) (https://github.com/fmwindsor/plastic-allometry)
se_body_length_intercept <- 0.3222      # SE generated from Jams et al R code (not rpoerted in paper or SI!) (https://github.com/fmwindsor/plastic-allometry)
#data
sim_beta_log10_body_length_samples <- rnorm(n_sim, mean = beta_log10_body_length, sd = se_beta_log10_body_length)
sim_body_length_intercept_samples <- rnorm(n_sim, mean = body_length_intercept, sd = se_body_length_intercept)

#### Tissue translocation length ##
# Parameters from the logistic regression# See scripts/translocation/translocation.Rmd for model
beta_0 <- 1.308344# simple$coefficients[[1]]  # 1.24  # Intercept
beta_1 <- -0.01468148 # simple$coefficients[[2]]   #-0.014  # Slope for particle length
se_beta_0 <- 0.3963612 # simple_summary$coefficients[[1,2]]  # 0.40  # SE of intercept
se_beta_1 <- 0.006657993 #simple_summary$coefficients[[2,2]]  # SE of slope

# Simulate beta_0 and beta_1
sim_beta_0 <- rnorm(n_sim * 1.2, mean = beta_0, sd = se_beta_0)
sim_beta_1 <- rnorm(n_sim * 1.2, mean = beta_1, sd = se_beta_1)

# Calculate X50 for each simulation
sim_X50 <- -sim_beta_0 / sim_beta_1

#truncate distribution to not fall below 0
upper.tissue.trans.size.um_samples <- sim_X50 %>% data.frame() %>% filter(.>0) %>% slice(1:n_sim)
upper.tissue.trans.size.um_samples <- as.numeric(upper.tissue.trans.size.um_samples$.)

#define param values
param_values <- data.frame(
  alpha.marine = alpha.marine_samples,
  a.sa.marine = a.sa.marine_samples,
  a.v.marine = a.v.marine_samples,
  a.m.marine = a.m.marine_samples,
  a.ssa.marine = a.ssa.marine_samples,
  R.ave.water.marine = R.ave.water.marine_samples,
  alpha.freshwater = alpha.freshwater_samples,
  a.sa.freshwater = a.sa.freshwater_samples,
  a.v.freshwater = a.v.freshwater_samples,
  a.m.freshwater = a.m.freshwater_samples,
  a.ssa.freshwater = a.ssa.freshwater_samples,
  R.ave.water.freshwater = R.ave.water.freshwater_samples,
  R.ave.water.freshwater = R.ave.water.freshwater_samples,
  R.ave.sediment.marine = R.ave.sediment.marine_samples,
  R.ave.sediment.freshwater = R.ave.sediment.freshwater_samples,
  sim_beta_log10_body_length = sim_beta_log10_body_length_samples,
  sim_body_length_intercept = sim_body_length_intercept_samples,
  upper.tissue.trans.size.um =  upper.tissue.trans.size.um_samples
)

#### Define Model Wrapper
model_wrapper <- function(params, iteration){
  
  # Print the current iteration
  print(paste("Iteration:", iteration))
  
  #unpack parameters
  # Ensure all extracted parameters are correctly coerced to numeric
  # Coerce to numeric directly while accessing the first element of potential list
  alpha.marine <- as.numeric(params$alpha.marine[1])
  a_sa.marine <- as.numeric(params$a.sa.marine[1])
  a_v.marine <- as.numeric(params$a.v.marine[1])
  a_m.marine <- as.numeric(params$a.m.marine[1])
  a_ssa.marine <- as.numeric(params$a.ssa.marine[1])
  alpha.freshwater <- as.numeric(params$alpha.freshwater[1])
  a_sa.freshwater <- as.numeric(params$a.sa.freshwater[1])
  a_v.freshwater <- as.numeric(params$a.v.freshwater[1])
  a_m.freshwater <- as.numeric(params$a.m.freshwater[1])
  a_ssa.freshwater <- as.numeric(params$a.ssa.freshwater[1])
  R_ave_water_marine <- as.numeric(params$R.ave.water.marine[1])
  R_ave_water_freshwater <- as.numeric(params$R.ave.water.freshwater[1])
  R_ave_sediment_marine <- as.numeric(params$R.ave.sediment.marine[1])
  R_ave_sediment_freshwater <- as.numeric(params$R.ave.sediment.freshwater[1])
  sim_beta_log10_body_length <- as.numeric(params$sim.beta.log10.body.length[1])
  sim_body_length_intercept <- as.numeric(params$sim.body.length.intercept[1])
  upper.tissue.trans.size.um <- as.numeric(params$upper.tissue.trans.size.um[1])
  
  
  ####### ---- RUN FUNCTIONS ---- ###
  ### Static with base parameters to ensure everything is working ###
  ## generate ToMEx 1.0 dataset ###
  aoc_setup <- ToMEx1.0fxn(R.ave.water.marine = R.ave.water.marine,
                           R.ave.water.freshwater = R.ave.water.freshwater,
                           R.ave.sediment.marine = R.ave.sediment.marine,
                           R.ave.sediment.freshwater = R.ave.sediment.freshwater,
                           beta_log10_body_length = beta_log10_body_length,
                           body_length_intercept = body_length_intercept) 
  
  ## generate ToMEx 2.0 dataset ###
  tomex2.0_aoc_z_final <- ToMEx2.0fxn(aoc_setup = aoc_setup,
                                      R.ave.water.marine = R.ave.water.marine,
                                      R.ave.water.freshwater = R.ave.water.freshwater,
                                      R.ave.sediment.marine = R.ave.sediment.marine,
                                      R.ave.sediment.freshwater = R.ave.sediment.freshwater,
                                      beta_log10_body_length = beta_log10_body_length,
                                      body_length_intercept = body_length_intercept)
  
  ### perform alignments
  aoc_MC_iter <- tomex2.0_aoc_z_final  %>% 
    ungroup() %>% 
    rename(environment = env_f) %>% 
    ### First filter the data ####
  ## First filter data with global filters
  filter(!environment %in% c("Terrestrial", "Not Reported"),
         Group != "Bacterium",
         Group != "Plant",
         effect.metric != "HONEC",
         tier_zero_tech_f == "Red Criteria Passed",
         tier_zero_risk_f == "Red Criteria Passed", #All thresholds must pass technical and risk red criteria
         risk.13 != 0 #Drop studies that received a score of 0 for endpoints criteria (this also drops studies that have not yet been scored) - KEEP THIS AFTER THE RED CRITERIA FILTERS  
  ) %>% 
    #Remove 26C temperature treatment data from Jaimukar et al. 2018
    filter(!(article == 42 & media.temp == 26)) %>% 
    mutate(max.size.ingest.um = 1000 * max.size.ingest.mm) %>%  #makes it less confusing below
    # calculate ERM for each species
    #### TISSUE TRANSLOCATION ####
  # define upper size length for Translocation 
  #set to 83um for upper limit or max size ingest, whichever is smaller
  mutate(x2M_trans = case_when(
    is.na(max.size.ingest.um) ~ upper.tissue.trans.size.um,
    TRUE ~ pmin(max.size.ingest.um, upper.tissue.trans.size.um)
  )) %>% 
    # define environment-specific alpha parameters #
    mutate(alpha = case_when(environment == "Marine" ~ alpha.marine,
                             environment == "Freshwater" ~ alpha.freshwater),
           a.sa = case_when(environment == "Marine" ~ a_sa.marine,
                            environment == "Freshwater" ~ a_sa.freshwater),
           a.v = case_when(environment == "Marine" ~ a_v.marine,
                            environment == "Freshwater" ~ a_v.freshwater),
           a.m = case_when(environment == "Marine" ~ a_m.marine,
                           environment == "Freshwater" ~ a_m.freshwater),
           a.ssa = case_when(environment == "Marine" ~ a_ssa.marine,
                           environment == "Freshwater" ~ a_ssa.freshwater),
           ) %>% 
    # calculate effect threshold for particles
    mutate(EC_mono_p.particles.mL_trans = dose.particles.mL.master) %>% 
    mutate(mu.p.mono = 1) %>% #mu_x_mono is always 1 for particles to particles
    mutate(mu.p.poly_trans = mux_polyfnx_generalizable(a.x = alpha, #alpha for particles
                                         x_UL= x2M_trans, #upper ingestible size limit (width of particle)
                                         x_LL = x1M_set)) %>% 
    # polydisperse effect threshold for particles
    mutate(EC_poly_p.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.p.mono)/mu.p.poly_trans) %>% 
    #calculate CF_bio for all conversions
    mutate(CF_bio_trans = CFfnx(x1M = x1M_set,#lower size bin
                                x2M = x2M_trans, #upper translocatable
                                x1D = x1D_set, #default
                                x2D = x2D_set,  #default
                                a = alpha)) %>%  
    ## Calculate environmentally relevant effect threshold for particles
    mutate(EC_env_p.particles.mL_trans = EC_poly_p.particles.mL_trans * CF_bio_trans) %>%  #aligned particle effect concentraiton (1-5000 um)
    
    #### Surface area ERM ####
  ##--- environmental calculations ---###
  #calculate lower translocatable surface area
  mutate(x_LL_sa_trans = SAfnx(a = 0.5 * x1D_set, #length
                               b = 0.5 * x1D_set, #0.5 * R.ave * x1D_set, #width
                               c = 0.5 * x1D_set  #0.5 * R.ave * 0.67 * x1D_set #height
  )) %>%  
    #calculate upper translocatable surface area
    mutate(x_UL_sa_trans = SAfnx(a = 0.5 * x2M_trans, 
                                 b = 0.5 * x2M_trans, #width #0.5 * R.ave * x2M, 
                                 c = 0.5 * x2M_trans #heigth #0.5 * R.ave * 0.67 * x2M
    )) %>%  
    #calculate mu_x_poly (env) for surface area
    mutate(mu.sa.poly_trans = mux_polyfnx_generalizable(a.sa, x_UL_sa_trans, x_LL_sa_trans)) %>% 
    
    ##--- laboratory calculations ---###
    ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
    #(note that if mixed particles were used, a different equation must be used)
    mutate(mu.sa.mono = case_when(
      polydispersity == "monodisperse" ~ particle.surface.area.um2, # use reported surface area in monodisperse
      polydispersity == "polydisperse" ~  mux_polyfnx_generalizable(a.x = a.sa, 
                                                      x_LL = particle.surface.area.um2.min,
                                                      x_UL = particle.surface.area.um2.max))) %>% 
    
    #calculate polydisperse effect concentration for surface area (particles/mL)
    mutate(EC_poly_sa.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.sa.mono)/mu.sa.poly_trans) %>%  
    #calculate environmentally realistic effect threshold
    mutate(EC_env_sa.particles.mL_trans = EC_poly_sa.particles.mL_trans * CF_bio_trans) %>% 
    
    ##### FOOD DILUTION ####
  # define upper size length for ingestion 
  mutate(x2M_ingest = case_when(is.na(max.size.ingest.um) ~ x2D_set, 
                                max.size.ingest.um < x2D_set ~ max.size.ingest.um,
                                max.size.ingest.um > x2D_set ~ x2D_set
  )) %>%  #set to 5,000 as upper limit or max size ingest, whichever is smaller
    # calculate effect threshold for particles
    mutate(EC_mono_p.particles.mL_ingest = dose.particles.mL.master) %>% 
    mutate(mu.p.mono = 1) %>% #mu_x_mono is always 1 for particles to particles
    mutate(mu.p.poly_ingest = mux_polyfnx_generalizable(a.x = alpha, #alpha for particles
                                          x_UL= x2M_ingest, #upper ingestible size limit
                                          x_LL = x1M_set)) %>% 
    # polydisperse effect threshold for particles
    mutate(EC_poly_p.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.p.mono)/mu.p.poly_ingest) %>% 
    #calculate CF_bio for all conversions
    mutate(CF_bio_ingest = CFfnx(x1M = x1M_set,#lower size bin
                                 x2M = x2M_ingest, #upper ingestible length
                                 x1D = x1D_set, #default
                                 x2D = x2D_set,  #default upper size range
                                 a = alpha)) %>%  
    ## Calculate environmentally relevant effect threshold for particles
    mutate(EC_env_p.particles.mL_ingest = EC_poly_p.particles.mL_ingest * CF_bio_ingest) %>%  #aligned particle effect concentraiton (1-5000 um)
    
    
    #### volume ERM ####
  ##--- environmental calculations ---###
  #calculate lower ingestible volume 
  mutate(x_LL_v_ingest = volumefnx_poly(length = x1D_set,
                                        width = x1D_set)) %>% 
    #calculate maximum ingestible volume 
    mutate(x_UL_v_ingest = volumefnx_poly(length = x2M_ingest, # length-limited
                                          #x2D_set, #upper definiton (accouunts for fibers) CONSERVATIVE
                                          width = x2M_ingest)) %>% #ingestion-limited
    # calculate mu.v.poly
    mutate(mu.v.poly_ingest = mux_polyfnx_generalizable(a.v, x_UL_v_ingest, x_LL_v_ingest)) %>% 
    ##--- laboratory calculations ---###
    ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
    #(note that if mixed particles were used, a different equation must be used)
    mutate(mu.v.mono = case_when(
      polydispersity == "monodisperse" ~ particle.volume.um3, # use reported volume in monodisperse
      polydispersity == "polydisperse" ~ mux_polyfnx_generalizable(a.x = a.v, 
                                                     x_LL = particle.volume.um3.min,
                                                     x_UL = particle.volume.um3.max))) %>% 
    
    #calculate polydisperse effect concentration for volume (particles/mL)
    mutate(EC_poly_v.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.v.mono)/mu.v.poly_ingest) %>%  
    #calculate environmentally realistic effect threshold
    mutate(EC_env_v.particles.mL_ingest = EC_poly_v.particles.mL_ingest * CF_bio_ingest) %>% 
    
    ###### CLEANUP #####
  mutate(particles.mL.ox.stress = EC_env_sa.particles.mL_trans,
         particles.mL.food.dilution = EC_env_v.particles.mL_ingest) %>% 
    #### annotate studies in which particles are too big to be ingested or translocated for those ERMs (to be filtered in ssd functions)
    mutate(translocatable = ifelse(size.length.um.used.for.conversions > x2M_trans, 
                                   "not translocatable", 
                                   "translocatable")) %>% 
    mutate(ingestible = ifelse(size.length.um.used.for.conversions > x2M_ingest, 
                               "not ingestible", 
                               "ingestible")) %>% 
    
    #rowwise() %>%
    mutate(unique_id = row_number()) %>% 
    # mutate(unique_id = digest::digest(paste(across(everything()), collapse = "-"), algo = "md5")) %>%
    ungroup()
  
  ##### Base Thresholds
  #filter out risk criteria (not done above)#
  aoc_risk_paper <- aoc_MC_iter %>%
    drop_na(effect.metric) %>%
    filter(tier_zero_tech_f == ("Red Criteria Passed"))
  
  # calculate thresholds for different environments
  marine_thresholds <- process_environment_data(aoc_risk_paper,
                                                "Marine", 
                                                upper.tissue.trans.size.um = upper.tissue.trans.size.um,
                                                x1D_set = x1D_set, 
                                                x2D_set = x2D_set)
  
  freshwater_thresholds <- process_environment_data(aoc_risk_paper, 
                                                    "Freshwater",
                                                    upper.tissue.trans.size.um = upper.tissue.trans.size.um,
                                                    x1D_set = x1D_set, 
                                                    x2D_set = x2D_set)
  
  freshwater_marine_thresholds <- process_environment_data(aoc_risk_paper, 
                                                           c("Freshwater", "Marine"),
                                                           upper.tissue.trans.size.um = upper.tissue.trans.size.um,
                                                           x1D_set = x1D_set, 
                                                           x2D_set = x2D_set)
  
  
  ##### SAVE OUTPUT OF MONTE CARLO ######
  
  MC_results[[i]] <- list(particles_mL_ox_stress = aoc_MC_iter$particles.mL.ox.stress,
                          particles_mL_food_dilution = aoc_MC_iter$particles.mL.food.dilution,
                          unique_id = aoc_MC_iter$unique_id,
                          base_thresholds = list(
                            marine = marine_thresholds,
                            freshwater = freshwater_thresholds,
                            freshwater_marine = freshwater_marine_thresholds
                          ))
} #close model_wrapper function

####### RUN MONTE CARLO #####

# Assume param_values is a dataframe where each row corresponds to a set of parameter samples
MC_results <- vector("list", nrow(param_values))

for (i in 1:nrow(param_values)) {
  # Extract a single row of parameter samples as a list
  param_set <- param_values[i, ]
  
  # Run the model function
  MC_results[[i]] <- model_wrapper(param_set, i)
}

# save results
saveRDS(MC_results, "scripts/monte carlo/output/MC_results.rds")


#### START HERE IF WORKING FROM SAVED FILE ###
MC_results <- readRDS("scripts/monte carlo/output/MC_results.rds")

##################################### ANALYSIS #########################################

# Extract and combine base_thresholds for each environment
all_thresholds_marine <- map(MC_results, ~ .x$base_thresholds$marine) %>% 
  bind_rows(.id = "simulation_id")

all_thresholds_freshwater <- map(MC_results, ~ .x$base_thresholds$freshwater) %>% 
  bind_rows(.id = "simulation_id")

# all_thresholds_freshwater_marine <- map(MC_results, ~ .x$base_thresholds$freshwater_marine) %>% 
#   bind_rows(.id = "simulation_id")

# Combine marine and freshwater data into one data frame with an additional "Environment" column
all_thresholds_combined <- bind_rows(
  mutate(all_thresholds_marine, Environment = "Marine"),
  mutate(all_thresholds_freshwater, Environment = "Freshwater"),
 # mutate(all_thresholds_freshwater_marine, Environment = "freshwater_marine"),
  .id = "simulation_id"
)

# Pivot data longer if the data structure requires it
# Assuming that your data might already be in a wide format and needs to be made long
all_thresholds_long <- all_thresholds_combined %>%
  pivot_longer(
    cols = -c(Tier, simulation_id, Environment),
    names_to = "Metric",
    values_to = "Value"
  )

# Calculate summary statistics
summary_stats_base_thresholds <- all_thresholds_long %>%
  group_by(Tier, Environment, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Std_Dev = sd(Value, na.rm = TRUE),
    N_sim = n(),
    .groups = 'drop'
  )

### Ggplot of SSD Tier Thresholds
# Prepare the data by calculating ERM and filtering
prepared_data <- all_thresholds_long %>%
  mutate(ERM = case_when(
    grepl("Tissue", Metric) ~ "Tissue Translocation",
    grepl("Food", Metric) ~ "Food Dilution"
  )) %>%
  filter(
    grepl("Default", Metric)
  )

# Calculate median values for each group
median_data <- prepared_data %>%
  group_by(Tier, ERM) %>%
  summarise(Median = median(Value), .groups = 'drop')

# Assign specific colors based on Tier levels
tier_colors <- setNames(c("#e2efd9", "#f9e39c", "#f0a95f", "#f0514b"), levels(prepared_data$Tier))

# define functino for making figures for marine and freshwater
MC_histograms_fxn <- function(environment){
  
  data <- prepared_data %>% 
    filter(Environment == environment)
  
  MC_histograms_base <- ggplot(data, aes(x = Value)) +
  geom_density(aes(color = Tier), size = 1) +  # Map 'Tier' for color
  geom_histogram(aes(y = ..density.., fill = Tier), bins = 150, color = "black", alpha = 0.6, linewidth = 0.05) +  # Map 'Tier' for fill
  geom_vline(data = median_data, aes(xintercept = Median),
             linetype = "dotted", color = "black", size = 1.2) +
#  geom_text(data = median_data, aes(x = Median, y = 0, label = paste0("Median:", scientific(Median, 2))),
        #    hjust = +1.5, vjust = -1.5, color = "black", size = 6, fontface = "italic") +
  xlab("Particles/L") +
  scale_x_log10(labels = scales::comma) +
  facet_wrap(~ Tier + ERM, scales = "free", ncol = 2) +
  scale_fill_manual(values = tier_colors) +  # Apply manual colors for fill
  scale_color_manual(values = tier_colors)  # Apply manual colors for lines

library(ggdark)
MC_histograms_dark <- MC_histograms_base +
  dark_theme_bw(base_size = 20) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "none",
    strip.text = element_blank(),  # This hides the facet wrap titles
    #  strip.text = element_text(face = "bold", size = 16, margin = margin(t = 1, b = 1)),
    panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )

MC_histograms_light <- MC_histograms_base +
  theme_bw(base_size = 20) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "none",
    strip.text = element_blank(),  # This hides the facet wrap titles
 #  strip.text = element_text(face = "bold", size = 16, margin = margin(t = 1, b = 1)),
    panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )

return(list(MC_histograms_light, MC_histograms_dark))
}

# run function 
MC_histograms_marine <- MC_histograms_fxn("Marine")
MC_histograms_freshwater <- MC_histograms_fxn("Freshwater")

MC_histograms_marine[1]
MC_histograms_freshwater[1]

ggsave("scripts/monte carlo/output/MC_histograms_freshwater.png",
       MC_histograms_freshwater[1], 
      dpi = 300,
      width = 12, height = 9, units = "in")

ggsave("scripts/monte carlo/output/MC_histograms_marine.png",
       MC_histograms_marine[1], 
       dpi = 300,
       width = 12, height = 9, units = "in")


ggsave("output/presentation_Figs/MC_histograms.png",
       MC_histograms_marine[2], 
       dpi = 300,
       width = 10, height = 5,
       units = "in")

### Ggplot of monte carlo SSD thresholds
all_thresholds_long %>% 
  filter(grepl("Default", Metric)) %>%
  mutate(ERM = case_when(
    grepl("Tissue", Metric) ~ "Tissue Translocation",
    grepl("Food", Metric) ~ "Food Dilution")) %>% 
  ggplot(aes(x = Value, y = Tier, color = Environment, fill = Environment)) +
  #geom_boxplot(alpha = 0.5) +
  geom_violin(alpha = 0.9,
              draw_quantiles = T) +
  cols4all::scale_fill_discrete_c4a_cat("carto.bold") +
  cols4all::scale_color_discrete_c4a_cat("carto.bold") +
  #geom_jitter() +
  scale_x_log10(name = "Particles/L (1 - 5000 um)",
                labels = scales::comma) +
  facet_grid(rows = vars(Environment), cols = vars(ERM),
             scales = "free") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

# Output the summary statistics
print(summary_stats_base_thresholds)

#save output
saveRDS(summary_stats_base_thresholds, "scripts/monte carlo/output/summary_stats_base_thresholds.rds")


##### Determine Coefficient of Variation for each output ###
### Step 1: Create a dataframe from the list
results_df <- do.call(rbind, lapply(MC_results, function(x) {
  data.frame(
    unique_id = x$unique_id,
    particles_L_ox_stress = x$particles_mL_ox_stress *1000,
    particles_L_food_dilution = x$particles_mL_food_dilution * 1000
  )
}))
skimr::skim(results_df)

## Step 2:
# Calculate CoV
cov_results <- results_df %>%
  group_by(unique_id) %>%
  summarise(sd_particles_L_ox_stress = sd(particles_L_ox_stress, na.rm = TRUE),
            mean_particles_L_ox_stress = mean(particles_L_ox_stress, na.rm = TRUE),
            median_particles_L_ox_stress = median(particles_L_ox_stress, na.rm = TRUE),
            CoV_ox_stress = sd_particles_L_ox_stress / mean_particles_L_ox_stress,
            #food dilution
            sd_particles_L_food_dilution = sd(particles_L_food_dilution, na.rm = TRUE),
            mean_particles_L_food_dilution = mean(particles_L_food_dilution, na.rm = TRUE),
            median_particles_L_food_dilution = median(particles_L_food_dilution, na.rm = TRUE),
            CoV_food_dilution = sd_particles_L_food_dilution / mean_particles_L_food_dilution
  )

# View the results
skimr::skim(cov_results)

# Step 3: merge with full dataset and deterministic dataset
aoc_MC <- merge(aoc_aligned_test, cov_results, by = "unique_id", all.x = TRUE)
aoc_MC <- left_join(aoc_MC, 
                    aoc_aligned_test %>% ungroup() %>%  
                      select(c(unique_id, particles.mL.food.dilution, particles.mL.ox.stress)) %>%
                      mutate(particles_L_food_dilution_deterministic = particles.mL.food.dilution * 1000,
                             particles_L_ox_stress_deterministic = particles.mL.ox.stress * 1000),
                    by = "unique_id"
) %>% 
  mutate(n_sim = n_sim)
#filter(size.length.um.used.for.conversions >= x1M_set) %>% #alignments not valid below 1 um, so filter them before they become problems later

#Save file for EcoTox Project
saveRDS(aoc_MC, "scripts/monte carlo/output/aoc_MC.rds")


##### CHECKS ####

# summary stats to ensure MC was succesful 
aoc_MC %>% 
  select(c(particles_L_ox_stress_deterministic, median_particles_L_ox_stress,
           particles_L_food_dilution_deterministic, median_particles_L_food_dilution)) %>% 
  mutate(ox_stress_ratio = particles_L_ox_stress_deterministic / median_particles_L_ox_stress,
         food_dilution_ratio = particles_L_food_dilution_deterministic / median_particles_L_food_dilution)

#define limits
min_val <- min(aoc_MC$CoV_food_dilution, na.rm = TRUE)
max_val <- max(aoc_MC$CoV_food_dilution, na.rm = TRUE)

# You may want to add some padding around the min and max
pad <- (max_val - min_val) * 0.05
xlims <- c(min_val - pad, max_val + pad)

hist_cv_food <- aoc_MC %>% 
  ggplot(aes(x = CoV_food_dilution, fill = polydispersity)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  xlim(xlims) +  # Set x limits
  scale_x_log10() +
  xlab("Coefficient of Variation (unitless)") +
  ylab("Toxicity Data Points")# +
#fill.type +
#  theme.light #+
#theme(axis.title.x = element_blank())

hist_cv_ox <- aoc_MC %>% 
  ggplot(aes(x = CoV_ox_stress, fill = polydispersity)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  #fill.type +
  scale_x_log10() +
  # xlim(xlims) +  # Set x limits
  xlab("Coefficient of Variation (unitless)") +
  ylab("Toxicity Data Points") +
  # theme.light +
  theme(legend.position = "none",
        axis.title.y = element_blank())

CoV_plots <-  ggpubr::ggarrange(hist_cv_food, hist_cv_ox,
                                labels = c("A", "B"),
                                ncol = 2,
                                common.legend = TRUE, legend = "top")

getwd()

ggsave(filename = "CoV_plots.jpg",
       dpi = 300,
       path = "output/Manuscript_Figs/", 
       plot = CoV_plots, width = 10, height = 8, units = "in")

CoV_plots



#############################################################################################################
############################################################################################################
############################################ Sensitivity Analysis ###########################################
#############################################################################################################
############################################################################################################

# For high number of iterations, need to simplify data saved to not overload the RAM
#### Define Model Wrapper
model_wrapper_sobol <- function(params, iteration, N, sobol_results, save_interval = 500){
  
  #report time
  start_time <- Sys.time()
  for (i in 1:N) {
    iter_start_time <- Sys.time()
  
  
  #unpack parameters
  # Ensure all extracted parameters are correctly coerced to numeric
  # Coerce to numeric directly while accessing the first element of potential list
    alpha.marine <- as.numeric(params$alpha.marine[1])
    a_sa.marine <- as.numeric(params$a.sa.marine[1])
    a_v.marine <- as.numeric(params$a.v.marine[1])
    a_m.marine <- as.numeric(params$a.m.marine[1])
    a_ssa.marine <- as.numeric(params$a.ssa.marine[1])
    alpha.freshwater <- as.numeric(params$alpha.freshwater[1])
    a_sa.freshwater <- as.numeric(params$a.sa.freshwater[1])
    a_v.freshwater <- as.numeric(params$a.v.freshwater[1])
    a_m.freshwater <- as.numeric(params$a.m.freshwater[1])
    a_ssa.freshwater <- as.numeric(params$a.ssa.freshwater[1])
    R_ave_water_marine <- as.numeric(params$R.ave.water.marine[1])
    R_ave_water_freshwater <- as.numeric(params$R.ave.water.freshwater[1])
    #R_ave_sediment_marine <- as.numeric(params$R.ave.sediment.marine[1])
 #  R_ave_sediment_freshwater <- as.numeric(params$R.ave.sediment.freshwater[1])
    sim_beta_log10_body_length <- as.numeric(params$sim.beta.log10.body.length[1])
    sim_body_length_intercept <- as.numeric(params$sim.body.length.intercept[1])
    upper.tissue.trans.size.um <- as.numeric(params$upper.tissue.trans.size.um[1])
    
  
  
  ####### ---- RUN FUNCTIONS ---- ###
  ### Static with base parameters to ensure everything is working ###
  ## generate ToMEx 1.0 dataset ###
  aoc_setup <- ToMEx1.0fxn(R.ave.water.marine = R.ave.water.marine,
                           R.ave.water.freshwater = R.ave.water.freshwater,
                           R.ave.sediment.marine = R.ave.sediment.marine,
                           R.ave.sediment.freshwater = R.ave.sediment.freshwater,
                           beta_log10_body_length = beta_log10_body_length,
                           body_length_intercept = body_length_intercept) 
  
  ## generate ToMEx 2.0 dataset ###
  tomex2.0_aoc_z_final <- ToMEx2.0fxn(aoc_setup = aoc_setup,
                                      R.ave.water.marine = R.ave.water.marine,
                                      R.ave.water.freshwater = R.ave.water.freshwater,
                                      R.ave.sediment.marine = R.ave.sediment.marine,
                                      R.ave.sediment.freshwater = R.ave.sediment.freshwater,
                                      beta_log10_body_length = beta_log10_body_length,
                                      body_length_intercept = body_length_intercept)
  
  ### perform alignments
  aoc_MC_iter <- tomex2.0_aoc_z_final  %>% 
    ungroup() %>% 
    rename(environment = env_f) %>% 
    ### First filter the data ####
  ## First filter data with global filters
  filter(!environment %in% c("Terrestrial", "Not Reported"),
         Group != "Bacterium",
         Group != "Plant",
         effect.metric != "HONEC",
         tier_zero_tech_f == "Red Criteria Passed",
         tier_zero_risk_f == "Red Criteria Passed", #All thresholds must pass technical and risk red criteria
         risk.13 != 0 #Drop studies that received a score of 0 for endpoints criteria (this also drops studies that have not yet been scored) - KEEP THIS AFTER THE RED CRITERIA FILTERS  
  ) %>% 
    #Remove 26C temperature treatment data from Jaimukar et al. 2018
    filter(!(article == 42 & media.temp == 26)) %>% 
    mutate(max.size.ingest.um = 1000 * max.size.ingest.mm) %>%  #makes it less confusing below
    # calculate ERM for each species
    #### TISSUE TRANSLOCATION ####
  # define upper size length for Translocation 
  #set to 83um for upper limit or max size ingest, whichever is smaller
  mutate(x2M_trans = case_when(
    is.na(max.size.ingest.um) ~ upper.tissue.trans.size.um,
    TRUE ~ pmin(max.size.ingest.um, upper.tissue.trans.size.um)
  )) %>% 
    # define environment-specific alpha parameters #
    mutate(alpha = case_when(environment == "Marine" ~ alpha.marine,
                             environment == "Freshwater" ~ alpha.freshwater),
           a.sa = case_when(environment == "Marine" ~ a_sa.marine,
                            environment == "Freshwater" ~ a_sa.freshwater),
           a.v = case_when(environment == "Marine" ~ a_v.marine,
                           environment == "Freshwater" ~ a_v.freshwater),
           a.m = case_when(environment == "Marine" ~ a_m.marine,
                           environment == "Freshwater" ~ a_m.freshwater),
           a.ssa = case_when(environment == "Marine" ~ a_ssa.marine,
                             environment == "Freshwater" ~ a_ssa.freshwater),
    ) %>% 
    
    # calculate effect threshold for particles
    mutate(EC_mono_p.particles.mL_trans = dose.particles.mL.master) %>% 
    mutate(mu.p.mono = 1) %>% #mu_x_mono is always 1 for particles to particles
    mutate(mu.p.poly_trans = mux_polyfnx_generalizable(a.x = alpha, #alpha for particles
                                         x_UL= x2M_trans, #upper ingestible size limit (width of particle)
                                         x_LL = x1M_set)) %>% 
    # polydisperse effect threshold for particles
    mutate(EC_poly_p.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.p.mono)/mu.p.poly_trans) %>% 
    #calculate CF_bio for all conversions
    mutate(CF_bio_trans = CFfnx(x1M = x1M_set,#lower size bin
                                x2M = x2M_trans, #upper translocatable
                                x1D = x1D_set, #default
                                x2D = x2D_set,  #default
                                a = alpha)) %>%  
    ## Calculate environmentally relevant effect threshold for particles
    mutate(EC_env_p.particles.mL_trans = EC_poly_p.particles.mL_trans * CF_bio_trans) %>%  #aligned particle effect concentraiton (1-5000 um)
    
    #### Surface area ERM ####
  ##--- environmental calculations ---###
  #calculate lower translocatable surface area
  mutate(x_LL_sa_trans = SAfnx(a = 0.5 * x1D_set, #length
                               b = 0.5 * x1D_set, #0.5 * R.ave * x1D_set, #width
                               c = 0.5 * x1D_set  #0.5 * R.ave * 0.67 * x1D_set #height
  )) %>%  
    #calculate upper translocatable surface area
    mutate(x_UL_sa_trans = SAfnx(a = 0.5 * x2M_trans, 
                                 b = 0.5 * x2M_trans, #width #0.5 * R.ave * x2M, 
                                 c = 0.5 * x2M_trans #heigth #0.5 * R.ave * 0.67 * x2M
    )) %>%  
    #calculate mu_x_poly (env) for surface area
    mutate(mu.sa.poly_trans = mux_polyfnx_generalizable(a.sa, x_UL_sa_trans, x_LL_sa_trans)) %>% 
    
    ##--- laboratory calculations ---###
    ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
    #(note that if mixed particles were used, a different equation must be used)
    mutate(mu.sa.mono = case_when(
      polydispersity == "monodisperse" ~ particle.surface.area.um2, # use reported surface area in monodisperse
      polydispersity == "polydisperse" ~  mux_polyfnx_generalizable(a.x = a.sa, 
                                                      x_LL = particle.surface.area.um2.min,
                                                      x_UL = particle.surface.area.um2.max))) %>% 
    
    #calculate polydisperse effect concentration for surface area (particles/mL)
    mutate(EC_poly_sa.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.sa.mono)/mu.sa.poly_trans) %>%  
    #calculate environmentally realistic effect threshold
    mutate(EC_env_sa.particles.mL_trans = EC_poly_sa.particles.mL_trans * CF_bio_trans) %>% 
    
    ##### FOOD DILUTION ####
  # define upper size length for ingestion 
  mutate(x2M_ingest = case_when(is.na(max.size.ingest.um) ~ x2D_set, 
                                max.size.ingest.um < x2D_set ~ max.size.ingest.um,
                                max.size.ingest.um > x2D_set ~ x2D_set
  )) %>%  #set to 5,000 as upper limit or max size ingest, whichever is smaller
    # calculate effect threshold for particles
    mutate(EC_mono_p.particles.mL_ingest = dose.particles.mL.master) %>% 
    mutate(mu.p.mono = 1) %>% #mu_x_mono is always 1 for particles to particles
    mutate(mu.p.poly_ingest = mux_polyfnx_generalizable(a.x = alpha, #alpha for particles
                                          x_UL= x2M_ingest, #upper ingestible size limit
                                          x_LL = x1M_set)) %>% 
    # polydisperse effect threshold for particles
    mutate(EC_poly_p.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.p.mono)/mu.p.poly_ingest) %>% 
    #calculate CF_bio for all conversions
    mutate(CF_bio_ingest = CFfnx(x1M = x1M_set,#lower size bin
                                 x2M = x2M_ingest, #upper ingestible length
                                 x1D = x1D_set, #default
                                 x2D = x2D_set,  #default upper size range
                                 a = alpha)) %>%  
    ## Calculate environmentally relevant effect threshold for particles
    mutate(EC_env_p.particles.mL_ingest = EC_poly_p.particles.mL_ingest * CF_bio_ingest) %>%  #aligned particle effect concentraiton (1-5000 um)
    
    
    #### volume ERM ####
  ##--- environmental calculations ---###
  #calculate lower ingestible volume 
  mutate(x_LL_v_ingest = volumefnx_poly(length = x1D_set,
                                        width = x1D_set)) %>% 
    #calculate maximum ingestible volume 
    mutate(x_UL_v_ingest = volumefnx_poly(length = x2M_ingest, # length-limited
                                          #x2D_set, #upper definiton (accouunts for fibers) CONSERVATIVE
                                          width = x2M_ingest)) %>% #ingestion-limited
    # calculate mu.v.poly
    mutate(mu.v.poly_ingest = mux_polyfnx_generalizable(a.v, x_UL_v_ingest, x_LL_v_ingest)) %>% 
    ##--- laboratory calculations ---###
    ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
    #(note that if mixed particles were used, a different equation must be used)
    mutate(mu.v.mono = case_when(
      polydispersity == "monodisperse" ~ particle.volume.um3, # use reported volume in monodisperse
      polydispersity == "polydisperse" ~ mux_polyfnx_generalizable(a.x = a.v, 
                                                     x_LL = particle.volume.um3.min,
                                                     x_UL = particle.volume.um3.max))) %>% 
    
    #calculate polydisperse effect concentration for volume (particles/mL)
    mutate(EC_poly_v.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.v.mono)/mu.v.poly_ingest) %>%  
    #calculate environmentally realistic effect threshold
    mutate(EC_env_v.particles.mL_ingest = EC_poly_v.particles.mL_ingest * CF_bio_ingest) %>% 
    
    ###### CLEANUP #####
  mutate(particles.mL.ox.stress = EC_env_sa.particles.mL_trans,
         particles.mL.food.dilution = EC_env_v.particles.mL_ingest) %>% 
    #### annotate studies in which particles are too big to be ingested or translocated for those ERMs (to be filtered in ssd functions)
    mutate(translocatable = ifelse(size.length.um.used.for.conversions > x2M_trans, 
                                   "not translocatable", 
                                   "translocatable")) %>% 
    mutate(ingestible = ifelse(size.length.um.used.for.conversions > x2M_ingest, 
                               "not ingestible", 
                               "ingestible")) %>% 
    
    #rowwise() %>%
    mutate(unique_id = row_number()) %>% 
    # mutate(unique_id = digest::digest(paste(across(everything()), collapse = "-"), algo = "md5")) %>%
    ungroup()
  
  ##### Base Thresholds
  #filter out risk criteria (not done above)#
  aoc_risk_paper <- aoc_MC_iter %>%
    drop_na(effect.metric) %>%
    filter(tier_zero_tech_f == ("Red Criteria Passed"))
  
  # calculate thresholds for different environments
  marine_thresholds <- process_environment_data(aoc_risk_paper,
                                                "Marine", 
                                                upper.tissue.trans.size.um = upper.tissue.trans.size.um,
                                                x1D_set = x1D_set, 
                                                x2D_set = x2D_set)
  
  freshwater_thresholds <- process_environment_data(aoc_risk_paper, 
                                                    "Freshwater",
                                                    upper.tissue.trans.size.um = upper.tissue.trans.size.um,
                                                    x1D_set = x1D_set, 
                                                    x2D_set = x2D_set)
  
  # freshwater_marine_thresholds <- process_environment_data(aoc_risk_paper, 
  #                                                          c("Freshwater", "Marine"),
  #                                                          upper.tissue.trans.size.um = upper.tissue.trans.size.um,
  #                                                          x1D_set = x1D_set, 
  #                                                          x2D_set = x2D_set)
  
  
  ##### SAVE OUTPUT OF MONTE CARLO ######
  
  sobol_results[[i]] <- list(
    # particles_mL_ox_stress = aoc_MC_iter$particles.mL.ox.stress,
    #                       particles_mL_food_dilution = aoc_MC_iter$particles.mL.food.dilution,
    #                       unique_id = aoc_MC_iter$unique_id,
                          base_thresholds = list(
                            marine = marine_thresholds,
                            freshwater = freshwater_thresholds#,
                           # freshwater_marine = freshwater_marine_thresholds
                          ))
  
  # Save results every save_interval iterations
  if (iteration %% save_interval == 0) {
    filename <- paste0("scripts/monte carlo/output/sobol_results_iteration_", iteration, ".rds")
    saveRDS(sobol_results, file = filename)
    cat(sprintf("Saved results for iteration %d\n", iteration))
    flush.console()
  }
  
  iter_end_time <- Sys.time()
  iter_time <- iter_end_time - iter_start_time
  elapsed_time <- iter_end_time - start_time
  remaining_iterations <- N - i
  remaining_time <- (elapsed_time / i) * remaining_iterations
  
  cat(sprintf("Iteration %d/%d: Time per iteration = %.2f secs, Estimated remaining time = %.2f secs\n",
              i, N, iter_time, remaining_time))
  }
  total_time <- Sys.time() - start_time
  cat(sprintf("Total time for %d iterations: %.2f secs\n", N, total_time))
  
  return(sobol_results)
} #close model_wrapper function





### static params ##
x1M_set <- 1 #um lower size for all alignments
x1D_set <- 1 #um lower size for all alignments
x2D_set <- 5000 #um

# Load the necessary packages
library(sensobol)
library(truncnorm)
### Application ###
# 1. Generate Samples for the Parameters: Generate samples for your parameters using a suitable sampling technique such as Sobol' sequences.
# Define the parameter names
params <- c("alpha.marine", "a.sa.marine", "a.v.marine", "a.m.marine", "a.ssa.marine", 
            "alpha.freshwater", "a.sa.freshwater", "a.v.freshwater", "a.m.freshwater", "a.ssa.freshwater", 
            "R.ave.water.marine", "R.ave.water.freshwater",
          #  "R.ave.sediment.marine", "R.ave.sediment.freshwater",
            "sim_beta_log10_body_length", "sim_body_length_intercept",
            "upper.tissue.trans.size.um")

# Number of samples
N <- 4
matrices <- c("A", "B", "AB", "BA")
first <- total <- "azzini"

# Generate the Sobol' sequence
mat <- sobol_matrices(N = N, 
                      params = params,
                      type = "LHS"
                      )

# Convert to data.table
# Convert to data.table
mat <- data.table::as.data.table(mat)



# Transform each column to the specified probability distribution
mat[, "alpha.marine" := rnorm(.N, mean = 2.07, sd = 0.07)]
mat[, "a.sa.marine" := rnorm(.N, mean = 1.5, sd = 0.009)]
mat[, "a.v.marine" := rnorm(.N, mean = 1.48, sd = 0.063)]
mat[, "a.m.marine" := rnorm(.N, mean = 1.32, sd = 0.009)]
mat[, "a.ssa.marine" := rnorm(.N, mean = 1.98, sd = 0.297)]
mat[, "alpha.freshwater" := rnorm(.N, mean = alpha.freshwater, sd = alpha.sd.freshwater)]
mat[, "a.sa.freshwater" := rnorm(.N, mean = a.sa.freshwater, sd = a.sa.sd.freshwater)]
mat[, "a.v.freshwater" := rnorm(.N, mean = a.v.freshwater, sd = a.v.sd.freshwater)]
mat[, "a.m.freshwater" := rnorm(.N, mean = a.m.freshwater, sd = a.m.sd.freshwater)]
mat[, "a.ssa.freshwater" := rnorm(.N, mean = a.ssa.freshwater, sd = a.ssa.sd.freshwater)]
mat[, "R.ave.water.marine" := rtruncnorm(.N, a = 0.0001, b= 0.9999, mean = 0.77, sd = 0.29)]
mat[, "R.ave.water.freshwater" := rtruncnorm(.N, a = 0.0001, b= 0.9999, mean = 0.67, sd = 0.28)]
#mat[, "R.ave.sediment.marine" := rtruncnorm(.N, a = 0.0001, b= 0.9999, mean = 0.75, sd = 0.30)]
#mat[, "R.ave.sediment.freshwater" := rtruncnorm(.N, a = 0.0001, b= 0.9999, mean = 0.70, sd = 0.33)]
mat[, "sim_beta_log10_body_length" := rnorm(.N, mean = 0.9341, sd = 0.1376)]
mat[, "sim_body_length_intercept" := rnorm(.N, mean = 1.1200, sd = 0.3222)]

# Parameters from logistic regression model
beta_0 <- 1.308344
beta_1 <- -0.01468148
se_beta_0 <- 0.3963612
se_beta_1 <- 0.006657993

# Simulate beta_0 and beta_1
sim_beta_0 <- rnorm(nrow(mat) * 1.2, mean = beta_0, sd = se_beta_0)
sim_beta_1 <- rnorm(nrow(mat) * 1.2, mean = beta_1, sd = se_beta_1)

# Calculate X50 for each simulation
sim_X50 <- -sim_beta_0 / sim_beta_1

# Truncate distribution to not fall below 0
upper.tissue.trans.size.um_samples <- sim_X50 %>% 
  data.frame() %>% 
  filter(. > x1M_set) %>% 
  slice(1:nrow(mat))
upper.tissue.trans.size.um_samples <- as.numeric(upper.tissue.trans.size.um_samples$.)

# Replace Sobol' matrix column for upper.tissue.trans.size.um
mat[, "upper.tissue.trans.size.um" := upper.tissue.trans.size.um_samples]

# Convert the data.table to a data.frame
mat <- as.data.frame(mat)

# save mat (necessary for stats later)
saveRDS(mat, file = "scripts/monte carlo/output/mat.rds")


## alternatrively, just use param_values derived above
#mat <- param_values[1:N,]

# 2. Run the Model for Each Sample Set: Run your model for each set of parameter samples.
# Initialize a list to store the results
# Initialize the results list and counter
sobol_results <- vector("list", nrow(mat))
N <- as.integer(nrow(mat))
N

# Loop through each row of the parameter matrix
for (i in 1:N) {
  param_set <- mat[i, ]
  
  # Perform the model wrapping with error handling
  sobol_results <- tryCatch({
    model_wrapper_sobol(param_set, i, N = N, sobol_results, save_interval = 2)
  }, error = function(e) {
    cat(sprintf("Error at iteration %d: %s\n", i, e$message))
    flush.console()
    sobol_results
  })
  
  # Print progress and flush the console output
  cat(sprintf("Processing iteration %d\n", i))
  flush.console()
}

# Save the final results
saveRDS(sobol_results, file = "scripts/monte carlo/output/sobol_results.rds")

# Print a message to indicate completion
cat("All iterations complete and final results saved.\n")
flush.console()

#### start here from saved file ###
sobol_results <- readRDS("scripts/monte carlo/output/sobol_results.rds")

# Extract and combine marine and freshwater thresholds
all_thresholds_marine_sobol <- map_dfr(sobol_results, ~ .x$base_thresholds$marine, .id = "simulation_id")
all_thresholds_freshwater_sobol <- map_dfr(sobol_results, ~ .x$base_thresholds$freshwater, .id = "simulation_id")

# 3. Extract the Output of Interest: Extract the relevant output from each model run (marine and freshwater thresholds).
# Convert results to a numeric vector
Y_marine_food_T3 <- as.numeric(all_thresholds_marine_sobol %>% filter(Tier == "Tier3") %>%  pull(`Food Dilution (Default)`))
Y_freshwater_food_T3 <- as.numeric(all_thresholds_freshwater_sobol %>% filter(Tier == "Tier3") %>%  pull(`Food Dilution (Default)`))

Y_marine_tissue_T3 <- as.numeric(all_thresholds_marine_sobol %>% filter(Tier == "Tier3") %>%  pull(`Tissue Translocation (Default)`))
Y_freshwater_tissue_T3 <- as.numeric(all_thresholds_freshwater_sobol %>% filter(Tier == "Tier3") %>%  pull(`Tissue Translocation (Default)`))


#histogram of results
hist(Y_marine_food_T3)
hist(Y_freshwater_food_T3)

#stopped MC at iteration = 2037, so need to subset mat to that value
#mat2 <- mat[1:2037,]



#scatter plot of parameter interactions
scatter_marine_food_t3 <- sensobol::plot_scatter(data = mat2, N = N, Y = Y_marine_food_T3, params = params, method = "bin") + ggtitle("Scatter Plot: Marine Food")
scatter_freshwater_food_t3 <- sensobol::plot_scatter(data = mat2, N = N, Y = Y_freshwater_food_T3, params = params, method = "bin") + ggtitle("Scatter Plot: Freshwater Food")
scatter_marine_tissue_t3 <- sensobol::plot_scatter(data = mat2, N = N, Y = Y_marine_tissue_T3, params = params, method = "bin") + ggtitle("Scatter Plot: Marine Tissue")
scatter_freshwater_tissue_t3 <- sensobol::plot_scatter(data = mat2, N = N, Y = Y_freshwater_tissue_T3, params = params, method = "bin") + ggtitle("Scatter Plot: Freshwater Tissue")
 
# Arrange the plots into a 2x2 matrix with a common legend
scatterplots <- combined_plot <- ggpubr::ggarrange(
  scatter_marine_food_t3, scatter_freshwater_food_t3,
  scatter_marine_tissue_t3, scatter_freshwater_tissue_t3,
  ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"
)

scatterplots


ggsave(filename = "sobol_scatterplots.jpg",
       dpi = 300,
       path = "output/Manuscript_Figs", 
       plot = scatterplots, width = 14, height = 10, units = "in")



# multiscatter
multi_food <- plot_multiscatter(data = mat2, 
                  N = N,
                  Y = Y_marine_food_T3, 
                  params = c("alpha.marine", "a.v.marine", "R.ave.water.marine", "sim_body_length_intercept"))

multi_food

# multiscatter
multi_tissue <- plot_multiscatter(data = mat2, 
                                N = N,
                                Y = Y_marine_food_T3, 
                                params = c("alpha.marine", "a.sa.marine", "upper.tissue.trans.size.um"))

multi_tissue

#### checks ##
# Check for NAs in parameter values
any(is.na(mat))

# Check for NAs in model output
any(is.na(Y_marine_food_T3))
any(is.na(Y_freshwater_food_T3))
any(is.na(Y_marine_tissue_T3))
any(is.na(Y_freshwater_tissue_T3))

# Check variance of parameters
apply(mat, 2, var)

# Check variance of outputs
var(Y_marine_food_T3)
var(Y_freshwater_tissue_T3)

# Check for NAs in model output
valid_indices <- complete.cases(Y_marine_food_T3, Y_marine_tissue_T3, Y_freshwater_tissue_T3, Y_freshwater_food_T3)
Y_marine_food_T3_valid <- Y_marine_food_T3[valid_indices]
Y_marine_tissue_T3_valid <- Y_marine_tissue_T3[valid_indices]
Y_freshwater_food_T3_valid <- Y_freshwater_food_T3[valid_indices]
Y_freshwater_tissue_T3_valid <- Y_freshwater_tissue_T3[valid_indices]
mat_valid <- mat[valid_indices, ]


# Update N to the number of valid samples
N_valid <- as.numeric(nrow(mat_valid))
#N must be integer, with Y = 14 * N
N <- N_valid / 14
N <- as.integer(N)
N
length_y_rows <- N * 14
#139

R <- 2000 #bootstrap iterations

ind.marine_food_T3 <- sobol_indices(Y = Y_marine_food_T3_valid[1:length_y_rows], N = N, params = params, boot = TRUE, R = R, type = type, conf = conf)
ind.marine_tissue_T3 <- sobol_indices(Y = Y_marine_tissue_T3_valid[1:length_y_rows], N = N, params = params, boot = TRUE, R = R, type = type, conf = conf)
ind.freshwater_food_T3 <- sobol_indices(Y = Y_freshwater_food_T3_valid[1:length_y_rows], N = N, params = params, boot = TRUE, R = R, type = type, conf = conf)
ind.freshwater_tissue_T3 <- sobol_indices(Y = Y_freshwater_tissue_T3_valid[1:length_y_rows], N = N, params = params, boot = TRUE, R = R, type = type, conf = conf)


#The output of sobol_indices() is an S3 object of class sensobol, with the results stored in the component results. To improve the visualization of the object, we set the number of digits in each numerical column to 3:
cols <- colnames(ind.marine_food_T3$results)[1:5]
ind.marine_food_T3$results[,(cols):=round(.SD,3),.SDcols=(cols)]

ind.marine_food_T3
#The output informs of the first and total-order estimators used in the calculation, thetotal numberofmodel runsandthesumof thefirst-order indices: if (Pk i=1Si)<1, themodel is non-additive.

#WecanalsocomputetheSobol’ indicesofadummyparameter, i.e., aparameter thathas no influence on themodel output, to estimate thenumerical approximationerror. This will beused later onto identifyparameterswhose contributionto theoutputvariance is less than the approximationerror, and therefore cannot be considered influential. Like sobol_indices(),thefunctionsobol_dummy()allowstoobtainpointestimates(thedefault) orbootstrapestimates. Inthisexampleweusethelatteroption:
ind.marine_food_T3_dummy <- sobol_dummy(Y = Y_marine_food_T3_valid[1:length_y_rows], N = N, params = params, boot = TRUE, R = R)
ind.marine_tissue_T3_dummy <- sobol_dummy(Y = Y_marine_tissue_T3_valid[1:length_y_rows], N = N, params = params, boot = TRUE, R = R)
ind.freshwater_food_T3_dummy <- sobol_dummy(Y = Y_freshwater_food_T3_valid[1:length_y_rows], N = N, params = params, boot = TRUE, R = R)
ind.freshwater_tissue_T3_dummy <- sobol_dummy(Y = Y_freshwater_tissue_T3_valid[1:length_y_rows], N = N, params = params, boot = TRUE, R = R)


plot(ind.marine_food_T3, dummy = ind.marine_food_T3_dummy)
plot(ind.marine_tissue_T3, dummy = ind.marine_tissue_T3_dummy)

# 4. Calculate Sobol' Indices: Calculate the Sobol' indices using the results.
# Function to plot Sobol' indices
plot_indices <- function(indices, title) {
  indices$results %>% 
    mutate(sensitivity = case_when(sensitivity == "Si" ~ "First-Order Sobol Index",
                                   sensitivity == "Ti" ~ "Total-Order Sobol Index")) %>% 
  ggplot(aes(y = parameters, x = original, fill = sensitivity)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbarh(aes(xmin = low.ci, xmax = high.ci), 
                   position = position_dodge(0.9), height = 0.25,
                   color = "black") +
    labs(#title = title,
         #y = "Parameters",
         x = "Sobol' Index") +
    scale_y_discrete(position = "right") +  # Duplicate y-axis on the right side
    theme_minimal(base_size = 14) +
    theme(legend.title = element_blank())
}

# Plot Sobol' indices for marine thresholds
sobol_plot_marine_food <-  plot_indices(ind.marine_food_T3, "Sobol' Indices for Marine Food Thresholds")

# Plot Sobol' indices for freshwater thresholds
sobol_plot_freshwater_tissue <-  plot_indices(ind.freshwater_tissue_T3, "Sobol' Indices for Freshwater Tissue Thresholds")

# Plot Sobol' indices for freshwater thresholds
sobol_plot_marine_tissue <-  plot_indices(ind.marine_tissue_T3, "Sobol' Indices for Marine Tissue Thresholds")

# Plot Sobol' indices for freshwater thresholds
sobol_plot_freshwater_food <-  plot_indices(ind.freshwater_food_T3, "Sobol' Indices for Freshwater Food Thresholds")


library(ggpubr)
indices <- ggarrange(sobol_plot_marine_food + labs(title = "Food Dilution (Marine)") +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)),
          sobol_plot_freshwater_food + labs(title = "Food Dilution (Freshwater)") +
            theme(axis.text.y = element_blank(), axis.title.y = element_blank(),axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)),
          sobol_plot_marine_tissue + labs(title = "Tissue Translocation (Marine)") + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)),# + labs(y = "Tissue Translocation"),
          sobol_plot_freshwater_tissue + labs(title = "Tissue Translocation (Freshwater)") + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)),
          common.legend = TRUE, legend = "bottom",
          labels = c("A", "B", "C", "D"))

ggsave(filename = "sobol_indices.jpg",
       dpi = 300,
       path = "output/Manuscript_Figs", 
       plot = indices, width = 12, height = 9, units = "in")




#### TUTORIAL ###

N <- 100
params <- c("$r$", "$K$", "$N_0$")
matrices <- c("A", "B", "AB", "BA")
first <- total <- "azzini"
order <- "second"
R <- 10 ^ 3
type <- "percent" 
conf <- 0.95

population_growth <- function (r, K, X0) {
  X<- X0 
  for (i in 0:20) {
    X <- X + r * X * (1- X / K)}
  return(X)
}

population_growth_run <- function (dt) {return(mapply(population_growth, dt[, 1], dt[, 2], dt[, 3]))}
mat <- sobol_matrices(matrices = matrices, N = N, params = params, order = order, type = "LHS")

mat[, "$r$"] <- qnorm(mat[, "$r$"], 1.7, 0.3)
mat[, "$K$"] <- qnorm(mat[, "$K$"], 40, 1) 
mat[, "$N_0$"] <- qunif(mat[, "$N_0$"], 10, 50)

y <- population_growth_run(mat)
plot_uncertainty(Y = y, N = N) + labs(y = "Counts", x = "$y$")
plot_scatter(data = mat, N = N, Y = y, params = params, method = "bin")
plot_multiscatter(data = mat, N = N, Y = y, params = params, smpl = 2^11)
ind <- sobol_indices(matrices = matrices, Y = y, N = N, params = params, first = first, 
                     total = total, order = order, boot = TRUE, R = R, 
                     parallel = "no", type = type, conf = conf)
cols <- colnames(ind$results)[1:5]
ind$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
ind
ind.dummy<-sobol_dummy(Y=y,N=N,params=params,boot=TRUE,R=R)
plot(ind,dummy=ind.dummy)
plot(ind, order = "second")

# quasi-random numbers are assumed to be the safest bet when selecting a sampling algorithm for a function of unknown behavior.
# Generate samples
set.seed(123)
 
# example from chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/sensobol/vignettes/sensobol.pdf
#N <- 2 ^ 10 #sample size N of the base sample matrix
N = 100
k <- 12 #number of uncertaint parameters
params <- paste("$x_", 1:k, "$", sep = "") #vector with the parameters' name
R <- 10^3 #bootsrap number
type <- "norm"  #bootstrap confidence interval type
conf <- 0.95 #95% CI

#create the sample matrix
mat <- sobol_matrices(N = N, params = params, order = "first", type="QRN") #uses A, B, A^(i)B design with Sobnol's QRN numbers to compute first and total-order indices

## Once the sample matrix is defined we can run our model. Note that in mat each column is a model input and each row a sample point, so the model has to be coded as to run rowwise. This is already the case of the Sobol’ G function included in sensobol:
y <- sobol_Fun(mat)
# swiftly visualize model output by plotting a histogram of the model output obtained from the A matrix
plot_uncertainty(Y = y, N = N) + labs(y = "Counts", x = "$y$")

##Before computing Sobol’ indices, it is recommended to explore how the model output maps onto the model input space. 
# sensobol includes two functions to that aim, plot_scatter() and plot_multiscatter(). 
# The first displays the model output y against xi while showing the mean y value (i.e., as in Figures 1–2), and allows the user to identify patterns denoting sensitivity (Pianosi, Beven, Freer, Hall, Rougier, Stephenson, and Wagener 2016).
plot_scatter(data = mat, N = N, Y = y, params = params)

### The scatter plots in Figure 6 evidence that x1, x2 and x3 have more “shape” than the rest and thus have a higher influence on y than (x4,...,x8). However, scatter plots do not always permit to detect which parameters have a joint effect on the model output.
# To gain a first insight on these interactions, the function plot_multiscatter() plots xi against xj and maps the resulting coordinate to its respective model output value. Interactions are then visible by the emergence of colored patterns.

plot_multiscatter(data = mat, 
                  N = N,
                  Y = y, 
                  params = paste("$x_", 1:4, "$", sep = ""))
#The results suggest that x1 might interact with x2 given the colored pattern of the (x1,x2) facet: the highest values of the model output are concentrated in the corners of the (x1,x2) 
# input space and thus result from combinations of high/low x1 values with high/low x2 values.
# In case the analyst is interested in assessing the exact weight of this high-order interaction, the computation of second-order indices would be required.

##The last step is the computation of Sobol’ indices. In this specific case, we set boot = TRUE to bootstrap the Sobol’ indices and get confidence intervals:
ind <- sobol_indices(Y = y, N = N, params = params, boot = TRUE, R = R, type = type, conf = conf)
#The output of sobol_indices() is an S3 object of class sensobol, with the results stored in the component results. To improve the visualization of the object, we set the number of digits in each numerical column to 3:
cols <- colnames(ind$results)[1:5]
ind$results[,(cols):=round(.SD,3),.SDcols=(cols)]
ind
#The output informs of the first and total-order estimators used in the calculation, thetotal numberofmodel runsandthesumof thefirst-order indices: if (Pk i=1Si)<1, themodel is non-additive.

#WecanalsocomputetheSobol’ indicesofadummyparameter, i.e., aparameter thathas no influence on themodel output, to estimate thenumerical approximationerror. This will beused later onto identifyparameterswhose contributionto theoutputvariance is less than the approximationerror, and therefore cannot be considered influential. Like sobol_indices(),thefunctionsobol_dummy()allowstoobtainpointestimates(thedefault) orbootstrapestimates. Inthisexampleweusethelatteroption:
ind.dummy <- sobol_dummy(Y = y, N = N, params = params, boot = TRUE, R = R)

plot(ind, dummy = ind.dummy)


#############################################################################################################
############################################################################################################
######################################################## Parallel processing ##############################
##########################################################################################################
#########################################################################################################
# library(parallel)
# 
# # Identify loaded packages before creating the cluster
# loaded_packages <- loadedNamespaces()
# 
# numCores <- detectCores() - 2
# cl <- makeCluster(numCores)
# 
# # Export necessary variables and functions to the cluster
# other_vars <- c("param_values")
# all_functions <- Filter(is.function, mget(ls(globalenv())))
# function_names <- names(all_functions)
# vars_to_export <- c(other_vars, function_names, "loaded_packages")  # Include 'loaded_packages'
# clusterExport(cl, varlist = vars_to_export)
# 
# # Load all packages in each worker node
# clusterEvalQ(cl, lapply(loaded_packages, library, character.only = TRUE))
# 
# # Convert param_values rows to list of lists for input to parLapply
# param_list <- split(param_values, seq(nrow(param_values)))
# 
# # Execute model_wrapper function in parallel
# sobol_results <- parLapply(cl, param_list, function(params) {
#   model_wrapper(params)
# })
# 
# # Stop and close the cluster after use
# stopCluster(cl)


