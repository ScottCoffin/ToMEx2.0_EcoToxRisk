###############################################
#### MINIMAL REPREX — PSSD++ FRAMEWORK ####
###############################################

# PURPOSE
# -------
# This script provides a minimal, end-to-end reproducible example (REPREX)
# of the PSSD++ probabilistic ecotoxicological threshold framework used to
# derive microplastic (MP) hazard thresholds across environments and
# Effect–Response Metrics (ERMs).
#
# The workflow integrates:
#   (i) Monte Carlo (MC) uncertainty propagation,
#   (ii) probabilistic Species Sensitivity Distributions (pSSDs),
#   (iii) alignment of laboratory toxicity data to environmentally relevant
#        MP particle distributions, and
#   (iv) parallel computation with caching and progress reporting.
#
# The implementation mirrors the analytical pipeline described in
# Coffin et al. (2025) and related manuscripts, while remaining compact
# enough for debugging, method development, and demonstration.

# OVERVIEW OF WORKFLOW
# --------------------
# 1. Load required R packages (data wrangling, plotting, MC simulation,
#    parallelization, progress tracking).
#
# 2. Enforce a specific version of ssdtools (v0.3.7) to ensure consistency
#    with published SSD fitting behavior.
#
# 3. Source all custom utility, SSD, plotting, and PSSD functions required
#    for the analysis.
#
# 4. Import a preprocessed microplastics ecotoxicity dataset (ToMEx 2.0–derived),
#    apply quality-control filters, compute composite assessment factors (AFs),
#    and retain only studies meeting Tier 0 technical and risk criteria.
#
# 5. Generate synthetic uncertainty parameter sets using Latin Hypercube
#    Sampling (LHS) to support Monte Carlo simulation of ERM alignment
#    parameters (particle size bounds, truncation limits, uncertainty factors).
#
# 6. Perform MC-based alignment of toxicity data to environmentally relevant
#    MP exposure metrics (MC_sim_align_parallel), producing a simulation-ready
#    dataset spanning ERMs, environments, and uncertainty realizations.
#
# 7. Construct ERM-specific datasets for:
#       - Food Dilution
#       - Tissue Translocation
#    and stratify these by environment (Freshwater vs Marine) and threshold
#    tier (Tier 3/4 only in this REPREX).
#
# 8. Run the PSSD++ framework using make_all_pSSDs(), which:
#       - Iterates over tiers, environments, and ERMs
#       - Fits probabilistic SSDs using MC-simulated NOEC distributions
#       - Derives HC5 and HC10-based PNECs
#       - Generates and saves publication-ready figures
#       - Caches results to disk to avoid redundant recomputation
#       - Executes in parallel with live progress reporting
#
# 9. Aggregate resulting PNECs across all combinations into a summary table
#    suitable for visualization, reporting, or downstream risk assessment.
#
# KEY FEATURES
# ------------
# - Fully probabilistic (uncertainty propagated end-to-end)
# - Parallelized across Tier × Environment × ERM combinations
# - Deterministic caching for reproducibility and efficiency
# - Minimal simulation size (sim = 10) for demonstration purposes
#   (increase to ≥300 for production analyses)
#
# OUTPUTS
# -------
# - Cached PSSD results (.rds) for each Tier × Environment × ERM combination
# - Saved PSSD and PNEC figures in the specified output directory
# - A summarized table of PNECs across all model runs
#
# INTENDED USE
# ------------
# - Method development and debugging
# - Demonstration of the PSSD++ framework
# - Reproducible example accompanying manuscripts or technical reports
#
# NOTES
# -----
# - This script assumes all sourced helper functions are available locally.
# - Parallel progress is reported using the progressr + future framework.
# - Console verbosity is intentionally minimal for batch execution.
#
###############################################

# libraries
pkgs <- c(
  "tidyverse",
  "sensobol",
  "truncnorm",
  "ggpubr",
  "gtsummary",
  "doParallel",
  "doSNOW",
  "tictoc",
  "crayon",
  "Matrix",
  "reshape2",
  "future",
  "future.apply",
  "progressr"
)

invisible(lapply(pkgs, function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

#check SSDtools package version is 0.3.7#
source("scripts/utils/check_pckg.R")
check_and_install_version("ssdtools", "0.3.7")
library(ssdtools)

# import ssd functions
source("scripts/utils/ssd_functions.R")
#import custom functions
source("scripts/utils/functions.R")
# import plotting functions
source("scripts/utils/plotting_functions.R")
# PSSD functions
source("scripts/PSSD/rmore.r")
source("scripts/PSSD/do.pssd.r")

# import data
aoc_risk_paper <- readRDS("data/input/aoc_z_tomex2.RDS") %>% #dataset is prepped in shiny repo
  dplyr::rename(environment = env_f) %>%
  dplyr::mutate(AF.total = af.time * af.noec) %>% #calculate composite AF
  tidyr::drop_na(effect.metric) %>%
  filter(
    !environment %in% c("Terrestrial", "Not Reported"),
    Group != "Bacterium",
    Group != "Plant",
    effect.metric != "HONEC",
    tier_zero_tech_f == "Red Criteria Passed",
    tier_zero_risk_f == "Red Criteria Passed", #All thresholds must pass technical and risk red criteria
    risk.13 != 0 #Drop studies that received a score of 0 for endpoints criteria (this also drops studies that have not yet been)
  )

# generate synthetic parameter values for MC-SIM using Latin-Hypercube Sampling (uses parameter default values used in paper - modify as needed)
mat <- matrix_function(
  n = 10,
  params = param_default_values,
  upper.tissue.truncation.limit = 500,
  x1M_set = 1,
  x2D_set = 5000
)

# generate histograms to visualize parameter distributions
parameter_histograms <- parameter_histograms_function(mat)
parameter_histograms$alpha_combined_plot
parameter_histograms$tissue_body

### MC-SIM Align Data ###
MC_sim_df <- MC_sim_align_parallel(
  tox_data = aoc_risk_paper,
  param_matrix = mat,
  n_sim = 3,
  x1D_set = 1,
  x2D_set = 5000,
  num_cores = "auto"
)

### PSSD++ Simulation ##
########## food dilution ################
results_df_food <- MC_sim_df %>%
  dplyr::filter(
    ingestible != "not ingestible",
    particles_L_food_dilution > 0,
    Group != "Algae"
  ) %>%
  dplyr::mutate(dose_new_particles_L = particles_L_food_dilution) %>%
  tidyr::drop_na(particles_L_food_dilution)
# Tier 3/4
results_df_food_t3_t4 <- results_df_food %>%
  dplyr::filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
##### marine food ####
# T1/2
results_df_food_marine <- results_df_food %>%
  dplyr::filter(environment == "Marine")
# T3/4
results_df_food_t3_t4_marine <- results_df_food_t3_t4 %>%
  dplyr::filter(environment == "Marine")
##### freshwater food ####
# T1/2
results_df_food_freshwater <- results_df_food %>%
  dplyr::filter(environment == "Freshwater")
# T3/4
results_df_food_t3_t4_freshwater <- results_df_food_t3_t4 %>%
  dplyr::filter(environment == "Freshwater")
########## Tissue Translocation ################
results_df_tissue <- MC_sim_df %>%
  dplyr::filter(
    translocatable != "not translocatable",
    particles_L_ox_stress > 0
  ) %>%
  dplyr::mutate(dose_new_particles_L = particles_L_ox_stress) %>%
  tidyr::drop_na(particles_L_ox_stress)
# Tier 3/4
results_df_tissue_t3_t4 <- results_df_tissue %>%
  dplyr::filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
##### marine tissue ####
# T1/2
results_df_tissue_marine <- results_df_tissue %>%
  dplyr::filter(environment == "Marine")
# T3/4
results_df_tissue_t3_t4_marine <- results_df_tissue_t3_t4 %>%
  dplyr::filter(environment == "Marine")
##### freshwater tissue ####
# T1/2
results_df_tissue_freshwater <- results_df_tissue %>%
  filter(environment == "Freshwater")
# T3/4
results_df_tissue_t3_t4_freshwater <- results_df_tissue_t3_t4 %>%
  filter(environment == "Freshwater")

# ------------------------------------------------------------------
# ERM registry
# ------------------------------------------------------------------
erm_registry <- list(
  "Food Dilution" = list(
    base = results_df_food,
    t3_t4 = results_df_food_t3_t4
  ),
  "Tissue Translocation" = list(
    base = results_df_tissue,
    t3_t4 = results_df_tissue_t3_t4
  )
)

##### perform pSSD #####
pSSDs <- make_all_pSSDs(
  MC_sim_df, # MC-sim DF
  tiers = c(3), #mehinto et al (2022) threshold tiers to loop throughj
  environments = c("Freshwater", "Marine"), #environments to loop through
  erms = c("Food Dilution", "Tissue Translocation"), #ERMs to loop through
  sim = 3, # number of PSSD simulations
  cv_uf = 0.5, # coefficient of variation for uncertainty factors
  rmore_method = "step", # method to handle pSSD distribution building (options = 'step' i.e., original trapezoidal method in Wigger et al., 2020, or 'lognormal' - shortcut developed in Coffin et al. (2025))
  quantile_type = 8, # quantile type for distribution fitting (8 is default, and use of others may cause unexpected or unreliable results)
  debug_option = FALSE, # debugging option
  parallel = TRUE, # enable parallel processing
  workers = parallel::detectCores() - 1, # number of worker cores to use
  base_cache_dir = "package/test_output/pssd_cache/", # directory for caching PSSD results
  base_output_path = "package/test_output/figures/", # directory for saving PSSD figures
  overwrite_cache = FALSE # whether to overwrite existing cache files
)

# summarize PNECs into wide plot
PNEC_summary <- summarize_PNECs(pSSDs)

# Display the table
print(PNEC_summary)
