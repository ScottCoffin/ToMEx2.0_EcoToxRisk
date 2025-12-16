# Minimal beta test for ToMEx Monte Carlo alignment and PSSD++ utilities
# Runs a tiny subset with reduced nboot to validate the stitched functions.

library(tidyverse)
library(truncnorm)
library(ssdtools)

source("package/tomex_functions.R")

# Settings for quick validation
nboot <- 10
x1D_set <- 1
x2D_set <- 5000

# Inputs
aoc_z <- readRDS("data/input/aoc_z_tomex2.RDS")
aoc_sample <- aoc_z #%>% slice_head(n = 400) # small subset for speed

# Parameter draws (use first row for a deterministic-esque run)
params <- generate_data(n_sim = 5)
param_row <- params[1, ]

# Align data with sampled parameters
aligned <- align_data(
  aoc_sample,
  upper.tissue.trans.size.um_input = param_row$upper.tissue.trans.size.um,
  body_length_intercept_input = param_row$sim_body_length_intercept,
  beta_log10_body_length_input = param_row$sim_beta_log10_body_length,
  R.ave.marine_input = param_row$R.ave.water.marine,
  H_W_ratio.marine_input = param_row$H_W_ratio.marine,
  alpha.marine_input = param_row$alpha.marine,
  a.sa.marine_input = param_row$a.sa.marine,
  a.v.marine_input = param_row$a.v.marine,
  a.m.marine_input = param_row$a.m.marine,
  R.ave.freshwater_input = param_row$R.ave.water.freshwater,
  H_W_ratio.freshwater_input = param_row$H_W_ratio.freshwater,
  alpha.freshwater_input = param_row$alpha.freshwater,
  a.sa.freshwater_input = param_row$a.sa.freshwater,
  a.v.freshwater_input = param_row$a.v.freshwater,
  a.m.freshwater_input = param_row$a.m.freshwater
)

# Derive thresholds for both environments
thresholds_base <- process_environment_data(
  aligned,
  env_filter = c("Marine", "Freshwater"),
  upper.tissue.trans.size.um = param_row$upper.tissue.trans.size.um,
  x1D_set = x1D_set,
  x2D_set = x2D_set,
  nboot = nboot
) %>%
  convert_units_fxn(environment = "Marine", params = param_row, x1D_set = x1D_set, x2D_set = x2D_set)

print("Threshold summary (particles/L and converted units):")
print(thresholds_base)

# Full wrapper smoke test (Monte Carlo alignment + SSD per Mehinto tiers)
sobol_like <- model_wrapper_sobol_parallel(
  data = aoc_sample,
  params_row = param_row,
  simulation_id = "sim1",
  x1D_set = x1D_set,
  x2D_set = x2D_set,
  nboot = nboot
)

print("Wrapper output keys:")
print(names(sobol_like))
