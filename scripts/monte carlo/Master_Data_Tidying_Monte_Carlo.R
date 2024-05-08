####### MASTER DATA TIDYING - MONTE CARLO  ########
###### 04-25-2024 #####
### This script runs the data preparations in a Monte Carlo fashion ##

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
identical(aoc_setup, aoc_setup_OG)

tomex2.0_aoc_z_final_noRAve <- tomex2.0_aoc_z_final %>% ungroup() %>% select(-R.ave) 
                                
identical(tomex2.0_aoc_z_final_noRAve, #rowid is assigned randomly 
          tomex2.0_aoc_z_final_OG)

differences <- tomex2.0_aoc_z_final_noRAve  != tomex2.0_aoc_z_final_OG

# Assuming 'differences' is your logical matrix with TRUE/FALSE/NA values
# Set NAs to FALSE (if NA means no difference for your case)
differences[is.na(differences)] <- FALSE

# Get rows and columns indices where differences are TRUE
diff_indices <- which(differences, arr.ind = TRUE)

# Get unique rows and columns involved in differences
diff_rows <- unique(diff_indices[, "row"])
diff_cols <- unique(diff_indices[, "col"])

# If differences are found, print them out
if (length(diff_rows) > 0) {
  cat("Differences found in the following rows and columns:\n")
  for (i in seq_along(diff_rows)) {
    for (j in seq_along(diff_cols)) {
      if (differences[diff_rows[i], diff_cols[j]]) {  # Check if TRUE
        # Using column names instead of indices
        cat(sprintf("Difference at Row: %d, Column: %s\n",
                    diff_rows[i], colnames(differences)[diff_cols[j]]))
      }
    }
  }
} else {
  cat("No differences found.\n")
}

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
  mutate(mu.p.poly_trans = mux.polyfnx(a.x = alpha, #alpha for particles
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
  mutate(mu.sa.poly_trans = mux.polyfnx(a.sa, x_UL_sa_trans, x_LL_sa_trans)) %>% 
  
  ##--- laboratory calculations ---###
  ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
  #(note that if mixed particles were used, a different equation must be used)
  mutate(mu.sa.mono = case_when(
    polydispersity == "monodisperse" ~ particle.surface.area.um2, # use reported surface area in monodisperse
    polydispersity == "polydisperse" ~  mux.polyfnx(a.x = a.sa, 
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
  mutate(mu.p.poly_ingest = mux.polyfnx(a.x = alpha, #alpha for particles
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
  mutate(mu.v.poly_ingest = mux.polyfnx(a.v, x_UL_v_ingest, x_LL_v_ingest)) %>% 
  ##--- laboratory calculations ---###
  ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
  #(note that if mixed particles were used, a different equation must be used)
  mutate(mu.v.mono = case_when(
    polydispersity == "monodisperse" ~ particle.volume.um3, # use reported volume in monodisperse
    polydispersity == "polydisperse" ~ mux.polyfnx(a.x = a.v, 
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
n_sim <- 500 #monte carlo simulations

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


p.ave = 1.10 #average density in marine surface water
p.ave.sd = 0.14 #Tablse s3. Marine surfac water
p.ave_samples <- rnorm(n_sim, mean = p.ave, sd = p.ave.sd)

# Alignment properties
alpha = 2.07 #table s4 for marine surface water. length
alpha.sd = 0.07 #table s4 for marine surface water. lengthj
alpha_samples <- rnorm(n_sim, mean = alpha, sd = alpha.sd)

# define parameters for power law coefficients
a.sa = 1.5 #marine surface area power law
a.sa.sd = 0.009 #marine surface water surface area power law - table s4
a.sa_samples <- rnorm(n_sim, mean = a.sa, sd = a.sa.sd)

a.v = 1.48 #a_V for marine surface water volume
a.v.sd = 0.063
a.v_samples <- rnorm(n_sim, mean = a.v, sd = a.v.sd)

a.m = 1.32 # upper limit fora_m for mass for marine surface water in table S4 
a.m.sd = 0.009
a.m_samples <- rnorm(n_sim, mean = a.m, sd = a.m.sd)

a.ssa = 1.98 # A_SSA for marine surface water
a.ssa.sd = 0.297
a.ssa_samples <- rnorm(n_sim, mean = a.ssa, sd = a.ssa.sd)

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
  alpha = alpha_samples,
  a.sa = a.sa_samples,
  a.v = a.v_samples,
  a.m = a.m_samples,
  a.ssa = a.ssa_samples,
  R.ave.water.marine = R.ave.water.marine_samples,
  R.ave.water.freshwater = R.ave.water.freshwater_samples,
  R.ave.sediment.marine = R.ave.sediment.marine_samples,
  R.ave.sediment.freshwater = R.ave.sediment.freshwater_samples,
  sim_beta_log10_body_length = sim_beta_log10_body_length_samples,
  sim_body_length_intercept = sim_body_length_intercept_samples,
  upper.tissue.trans.size.um =  upper.tissue.trans.size.um_samples
)

#### Define Model Wrapper
model_wrapper <- function(params){
  
  #unpack parameters
  # Ensure all extracted parameters are correctly coerced to numeric
  # Coerce to numeric directly while accessing the first element of potential list
  alpha <- as.numeric(params$alpha[1])
  a_sa <- as.numeric(params$a.sa[1])
  a_v <- as.numeric(params$a.v[1])
  a_m <- as.numeric(params$a.m[1])
  a_ssa <- as.numeric(params$a.ssa[1])
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
    
    # calculate effect threshold for particles
    mutate(EC_mono_p.particles.mL_trans = dose.particles.mL.master) %>% 
    mutate(mu.p.mono = 1) %>% #mu_x_mono is always 1 for particles to particles
    mutate(mu.p.poly_trans = mux.polyfnx(a.x = alpha, #alpha for particles
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
    mutate(mu.sa.poly_trans = mux.polyfnx(a.sa, x_UL_sa_trans, x_LL_sa_trans)) %>% 
    
    ##--- laboratory calculations ---###
    ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
    #(note that if mixed particles were used, a different equation must be used)
    mutate(mu.sa.mono = case_when(
      polydispersity == "monodisperse" ~ particle.surface.area.um2, # use reported surface area in monodisperse
      polydispersity == "polydisperse" ~  mux.polyfnx(a.x = a.sa, 
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
    mutate(mu.p.poly_ingest = mux.polyfnx(a.x = alpha, #alpha for particles
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
    mutate(mu.v.poly_ingest = mux.polyfnx(a.v, x_UL_v_ingest, x_LL_v_ingest)) %>% 
    ##--- laboratory calculations ---###
    ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
    #(note that if mixed particles were used, a different equation must be used)
    mutate(mu.v.mono = case_when(
      polydispersity == "monodisperse" ~ particle.volume.um3, # use reported volume in monodisperse
      polydispersity == "polydisperse" ~ mux.polyfnx(a.x = a.v, 
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
  MC_results[[i]] <- model_wrapper(param_set)
}



##################################### ANALYSIS #########################################

# Extract and combine base_thresholds for each environment
all_thresholds_marine <- map(MC_results, ~ .x$base_thresholds$marine) %>% 
  bind_rows(.id = "simulation_id")

all_thresholds_freshwater <- map(MC_results, ~ .x$base_thresholds$freshwater) %>% 
  bind_rows(.id = "simulation_id")

all_thresholds_freshwater_marine <- map(MC_results, ~ .x$base_thresholds$freshwater_marine) %>% 
  bind_rows(.id = "simulation_id")

# Combine marine and freshwater data into one data frame with an additional "Environment" column
all_thresholds_combined <- bind_rows(
  mutate(all_thresholds_marine, Environment = "Marine"),
  mutate(all_thresholds_freshwater, Environment = "Freshwater"),
  mutate(all_thresholds_freshwater_marine, Environment = "freshwater_marine"),
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
    grepl("Default", Metric),
    Environment == "Marine"
  )

# Calculate median values for each group
median_data <- prepared_data %>%
  group_by(Tier, ERM) %>%
  summarise(Median = median(Value), .groups = 'drop')

# Assign specific colors based on Tier levels
tier_colors <- setNames(c("#e2efd9", "#f9e39c", "#f0a95f", "#f0514b"), levels(prepared_data$Tier))

MC_histograms <- ggplot(prepared_data, aes(x = Value)) +
  geom_density(aes(color = Tier), size = 1) +  # Map 'Tier' for color
  geom_histogram(aes(y = ..density.., fill = Tier), bins = 150, color = "black", alpha = 0.6, linewidth = 0.05) +  # Map 'Tier' for fill
  geom_vline(data = median_data, aes(xintercept = Median),
             linetype = "dotted", color = "black", size = 1.2) +
#  geom_text(data = median_data, aes(x = Median, y = 0, label = paste0("Median:", scientific(Median, 2))),
        #    hjust = +1.5, vjust = -1.5, color = "black", size = 6, fontface = "italic") +
  xlab("Particles/L") +
  scale_x_log10(labels = comma) +
  facet_wrap(~ Tier + ERM, scales = "free", ncol = 2) +
  scale_fill_manual(values = tier_colors) +  # Apply manual colors for fill
  scale_color_manual(values = tier_colors) +  # Apply manual colors for lines
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

MC_histograms

ggsave("scripts/monte carlo/output/MC_histograms.png",
       MC_histograms, 
      dpi = 300,
      width = 12, height = 9, units = "in")

### Ggplot of monte carlo SSD thresholds
all_thresholds_long %>% 
  filter(grepl("Default", Metric)) %>%
  mutate(ERM = case_when(
    grepl("Tissue", Metric) ~ "Tissue Translocation",
    grepl("Food", Metric) ~ "Food Dilution")) %>% 
  ggplot(aes(x = Value, y = Tier, color = Environment)) +
  geom_boxplot() +
 # geom_jitter() +
  scale_x_log10() +
  facet_grid(rows = vars(Environment), cols = vars(ERM)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

# Output the summary statistics
print(summary_stats_base_thresholds)

#save output
saveRDS(MC_results, "scripts/monte carlo/output/MC_results.rds")
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
# library(sensitivity)
# 
# #### Step 2: Generate Samples for Sensitivity Analysis ###
# # Define parameter ranges based on your existing data (assuming a normal distribution around your samples)
# param.dists <- list(
#   alpha = runif(100, min(alpha_samples), max(alpha_samples)),
#   a.sa = runif(100, min(a.sa_samples), max(a.sa_samples)),
#   a.v = runif(100, min(a.v_samples), max(a.v_samples)),
#   a.m = runif(100, min(a.m_samples), max(a.m_samples)),
#   a.ssa = runif(100, min(a.ssa_samples), max(a.ssa_samples)),
#   R.ave.water.marine = runif(100, min(R.ave.water.marine_samples), max(R.ave.water.marine_samples)),
#   R.ave.water.freshwater = runif(100, min(R.ave.water.freshwater_samples), max(R.ave.water.freshwater_samples)),
#   R.ave.sediment.marine = runif(100, min(R.ave.sediment.marine_samples), max(R.ave.sediment.marine_samples)),
#   R.ave.sediment.freshwater = runif(100, min(R.ave.sediment.freshwater_samples), max(R.ave.sediment.freshwater_samples)),
#   sim_beta_log10_body_length = runif(100, min(sim_beta_log10_body_length_samples), max(sim_beta_log10_body_length_samples)),
#   sim_body_length_intercept = runif(100, min(sim_body_length_intercept_samples), max(sim_body_length_intercept_samples))
# )
# 
# param_ranges <- list(
#   alpha = seq(min(alpha_samples), max(alpha_samples), length.out = 100),
#   a_sa = seq(min(a.sa_samples), max(a.sa_samples), length.out = 100),
#   a_v = seq(min(a.v_samples), max(a.v_samples), length.out = 100),
#   a_m = seq(min(a.m_samples), max(a.m_samples), length.out = 100),
#   a_ssa = seq(min(a.ssa_samples), max(a.ssa_samples), length.out = 100),
#   R_ave_water_marine = seq(min(R.ave.water.marine_samples), max(R.ave.water.marine_samples), length.out = 100),
#   R_ave_water_freshwater = seq(min(R.ave.water.freshwater_samples), max(R.ave.water.freshwater_samples), length.out = 100),
#   R_ave_sediment_marine = seq(min(R.ave.sediment.marine_samples), max(R.ave.sediment.marine_samples), length.out = 100),
#   R_ave_sediment_freshwater = seq(min(R.ave.sediment.freshwater_samples), max(R.ave.sediment.freshwater_samples), length.out = 100),
#   sim_beta_log10_body_length = seq(min(sim_beta_log10_body_length_samples), max(sim_beta_log10_body_length_samples), length.out = 100),
#   sim_body_length_intercept = seq(min(sim_body_length_intercept_samples), max(sim_body_length_intercept_samples), length.out = 100)
# )
# 
# # Generate sampling matrices
# library(randtoolbox)
# n <- 100  # number of samples
# X1 <- sobol(n, length(param.dists), scrambling = 0)
# X2 <- sobol(n, length(param.dists), scrambling = 0)
# 
# # Scale X1 and X2 to the parameter ranges
# scale01 <- function(x, min, max) min + (max - min) * x
# X1 <- apply(X1, 2, scale01, min = sapply(param.dists, min), max = sapply(param.dists, max))
# X2 <- apply(X2, 2, scale01, min = sapply(param.dists, min), max = sapply(param.dists, max))
# 
# #Step 3: Run the Model with Generated Samples
# # Define the names corresponding to each index in your param.dists
# param_names <- names(unlist(param.dists))
# 
# # Wrap the call in another function to set names on the vector
# Y1 <- apply(X1, 1, function(x) model_wrapper(setNames(x, param_names)))
# Y2 <- apply(X2, 1, function(x) model_wrapper(setNames(x, param_names)))



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
# MC_results <- parLapply(cl, param_list, function(params) {
#   model_wrapper(params)
# })
# 
# # Stop and close the cluster after use
# stopCluster(cl)


