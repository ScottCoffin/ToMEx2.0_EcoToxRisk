#################################
#### MINIMAL REPREX PSSD++ #####
#################################

# libraries
library(tidyverse)
library(sensobol)
library(truncnorm)
library(ggpubr)
library(ggdark)
library(gtsummary)
library(doParallel)
library(doSNOW)
library(tictoc)
library(crayon)
library(Matrix)
library(reshape2)
library(future)
library(future.apply)
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
  mutate(AF.total = af.time * af.noec) %>% #calculate composite AF
  drop_na(effect.metric) %>%
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
  n_sobol = 10,
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
  n_sim = 10,
  x1D_set = 1,
  x2D_set = 5000,
  num_cores = "auto"
)

### PSSD++ Simulation ##

########## food dilution ################
results_df_food <- MC_sim_df %>%
  filter(
    ingestible != "not ingestible",
    particles_L_food_dilution > 0,
    Group != "Algae"
  ) %>%
  mutate(dose_new_particles_L = particles_L_food_dilution) %>%
  drop_na(particles_L_food_dilution)

# Tier 3/4
results_df_food_t3_t4 <- results_df_food %>%
  filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))

##### marine food ####
# T1/2
results_df_food_marine <- results_df_food %>%
  filter(environment == "Marine")
# T3/4
results_df_food_t3_t4_marine <- results_df_food_t3_t4 %>%
  filter(environment == "Marine")

##### freshwater food ####
# T1/2
results_df_food_freshwater <- results_df_food %>%
  filter(environment == "Freshwater")

# T3/4
results_df_food_t3_t4_freshwater <- results_df_food_t3_t4 %>%
  filter(environment == "Freshwater")


########## Tissue Translocation ################
results_df_tissue <- MC_sim_df %>%
  filter(translocatable != "not translocatable", particles_L_ox_stress > 0) %>%
  mutate(dose_new_particles_L = particles_L_ox_stress) %>%
  drop_na(particles_L_ox_stress)

# Tier 3/4
results_df_tissue_t3_t4 <- results_df_tissue %>%
  filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))

##### marine tissue ####
# T1/2
results_df_tissue_marine <- results_df_tissue %>%
  filter(environment == "Marine")
# T3/4
results_df_tissue_t3_t4_marine <- results_df_tissue_t3_t4 %>%
  filter(environment == "Marine")

##### freshwater tissue ####
# T1/2
results_df_tissue_freshwater <- results_df_tissue %>%
  filter(environment == "Freshwater")

# T3/4
results_df_tissue_t3_t4_freshwater <- results_df_tissue_t3_t4 %>%
  filter(environment == "Freshwater")

# ------------------------------------------------------------------
# ERM registry: canonical source of truth
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

# function to generate PSSDs and plots with caching and optional debugging
make_all_pSSDs <- function(
  MC_sim_df = NULL, # MC-sim dataframe
  tiers = c(3), #mehinto et al tiers
  environments = c("Freshwater", "Marine"), #environments to loop through
  erms = c("Food Dilution", "Tissue Translocation"), #ERMs to loop through
  sim = 10,
  cv_uf = 0.5,
  rmore_method = "step",
  quantile_type = 8,
  debug_option = FALSE,
  parallel = TRUE,
  workers = parallel::detectCores() - 1,
  base_cache_dir = "package/test_output/pssd_cache/",
  base_output_path = "package/test_output/figures/",
  overwrite_cache = FALSE
) {
  # ---- Argument validation ----
  if (is.null(MC_sim_df)) {
    stop("MC_sim_df must be provided explicitly.", call. = FALSE)
  }
  # Create output directories if they don't exist
  output_path <- paste0(base_output_path, rmore_method, "_", sim, "sims")
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  cache_dir <- paste0(base_cache_dir, rmore_method, "_", sim, "sims")
  # Create cache directory if it doesn't exist
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = T)
  }
  # generate color palette for unique levels of species/groups
  all_species <- MC_sim_df %>%
    filter(
      environment %in% environments,
      !Group %in% c("Insect", "Annelida") #not used in this analysis
    ) %>%
    distinct(Species, Group, environment) %>%
    droplevels() %>%
    arrange(environment, Group, Species)

  # Generate a consistent color palette
  global_color_palette <- generate_color_palette(all_species)

  # Initialize tracking
  total_iterations <- length(tiers) * length(environments) * length(erms)
  iteration_times <- numeric(total_iterations)
  iteration_count <- 0
  tictoc::tic("PSSD++ For Loop Begins...")

  # Start loop
  results <- list()
  for (tier in tiers) {
    for (environment in environments) {
      for (erm in erms) {
        iteration_count <- iteration_count + 1
        iteration_start_time <- Sys.time()

        combo_id <- paste0("Tier", tier, "_", environment, "_", erm)
        cache_file <- file.path(cache_dir, paste0(combo_id, ".rds"))

        if (file.exists(cache_file) && !overwrite_cache) {
          cat(blue(sprintf(
            "\n🔁 [%s] Using cached result: %s\n",
            format(iteration_start_time, "%H:%M:%S"),
            combo_id
          )))
          results[[combo_id]] <- readRDS(cache_file)
          next
        }

        if (file.exists(cache_file) && overwrite_cache) {
          cat(yellow(sprintf(
            "\n♻️  [%s] Overwriting cached result: %s\n",
            format(iteration_start_time, "%H:%M:%S"),
            combo_id
          )))
        }

        cat(green(sprintf(
          "\n⏳ [%s] Starting: Tier %d | Environment: %s | ERM: %s\n",
          format(iteration_start_time, "%H:%M:%S"),
          tier,
          environment,
          erm
        )))

        # Wrap in tryCatch for error handling
        tryCatch(
          {
            results_df <- get_results_df(tier, environment, erm)

            result <- generate_plots_and_summary(
              tier = tier,
              environment = environment,
              erm = erm,
              color_palette = global_color_palette,
              results_df = results_df,
              sim = sim,
              num_iterations = num_iterations,
              quantile_type = quantile_type,
              cv_dp = NULL,
              cv_uf = cv_uf,
              rmore_method = rmore_method,
              species_data_source = MC_sim_df,
              output_path = output_path,
              presentation_path = presentation_path,
              debug = debug_option
            )

            saveRDS(result, cache_file)
            results[[combo_id]] <- result

            iteration_end_time <- Sys.time()
            iteration_time_sec <- as.numeric(difftime(
              iteration_end_time,
              iteration_start_time,
              units = "secs"
            ))
            iteration_times[iteration_count] <- iteration_time_sec

            running_average_sec <- mean(iteration_times[1:iteration_count])
            est_time_remaining_sec <- running_average_sec *
              (total_iterations - iteration_count)

            cat(yellow(sprintf(
              "✅ Completed: %s in %.2f hr | Est. remaining: %.2f hr\n",
              combo_id,
              iteration_time_sec / 3600,
              est_time_remaining_sec / 3600
            )))
          },
          error = function(e) {
            cat(red(sprintf("❌ ERROR in %s: %s\n", combo_id, e$message)))
          }
        )
      }
    }
  }
  return(results)

  # Wrap up
  pSSD_time <- tictoc::toc()
  total_runtime_sec <- pSSD_time$toc - pSSD_time$tic
  cat(blue(sprintf(
    "\n🎉 PSSD++ complete for all combinations! Total time: %.2f hours\n",
    total_runtime_sec / 3600
  )))
}

pSSDs <- make_all_pSSDs(MC_sim_df, overwrite_cache = TRUE)

# Initialize an empty list to store the standardized summaries
PNEC_summaries <- list()

# Loop through the results list to extract and standardize PNEC summaries
for (combination in names(results)) {
  # Extract the PNEC summary statistics for HC5 (hcx = 0.05)
  PNEC_stats_05 <- results[[combination]]$summary_05$stats
  PNEC_stats_05_long <- data.frame(
    Statistic = rownames(PNEC_stats_05),
    Value = as.numeric(PNEC_stats_05[, 1]), # Extract the first column (values)
    HCX = "HC5", # Add a column to indicate HC5
    Tier = "Tier3" # Assign Tier3 for HC5
  )

  # Extract the PNEC summary statistics for HC10 (hcx = 0.10)
  PNEC_stats_10 <- results[[combination]]$summary_10$stats
  PNEC_stats_10_long <- data.frame(
    Statistic = rownames(PNEC_stats_10),
    Value = as.numeric(PNEC_stats_10[, 1]), # Extract the first column (values)
    HCX = "HC10", # Add a column to indicate HC10
    Tier = "Tier4" # Assign Tier4 for HC10
  )

  # Combine HC5 and HC10 into a single data frame
  PNEC_stats_long <- rbind(PNEC_stats_05_long, PNEC_stats_10_long)

  # Add the combination name and split it into environment and ERM
  PNEC_stats_long$Combination <- combination
  PNEC_stats_long$Environment <- sub(".*_(.*)_(.*)", "\\1", combination) # Extract environment
  PNEC_stats_long$ERM <- sub(".*_(.*)", "\\1", combination) # Extract ERM

  # Append to the list
  PNEC_summaries[[combination]] <- PNEC_stats_long
}

# Combine all standardized summaries into a single data frame
PNEC_summary_table <- do.call(rbind, PNEC_summaries)

# Reorder columns for clarity
PNEC_summary_table <- PNEC_summary_table[, c(
  "Tier",
  "Environment",
  "ERM",
  "HCX",
  "Statistic",
  "Value"
)]

# Sort the table by Environment, ERM, Tier, and HCX
PNEC_summary_table <- PNEC_summary_table[
  order(
    PNEC_summary_table$Environment,
    PNEC_summary_table$ERM,
    PNEC_summary_table$Tier,
    PNEC_summary_table$HCX
  ),
]

# Remove row names
rownames(PNEC_summary_table) <- NULL

# Pivot the table wider
PNEC_summary_table_wide <- PNEC_summary_table %>%
  pivot_wider(
    names_from = c(Statistic), # Use both Statistic and HCX for new column names
    values_from = Value # Use the Value column for the values
  ) %>%
  select(-HCX)

# Save the table as an RDS file
saveRDS(PNEC_summary_table_wide, "../data/output/PNEC_summary_table_wide.rds")

# Display the table in R
print(PNEC_summary_table_wide)
