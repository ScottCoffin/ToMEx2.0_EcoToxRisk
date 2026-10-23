# Effect of the number of bootstrap samples (nboot) on the SSD thresholds.
#
# The analysis (and the ToMEx 2.0 app) uses nboot = 10. Bootstrapping affects
# only the confidence limits: Tier 1 (= lower 95% confidence limit of the
# HC5) and the 95% intervals of the traditional SSD thresholds. HC point
# estimates (Tiers 2-4) do not depend on it. Running the full MC + SSD
# analysis with nboot = 1000 is too costly, so this script tests the effect
#   A. on the traditional (deterministic) SSD: all thresholds and intervals
#      with nboot = 10 vs 1000, and the spread of Tier 1 over 20 seeds at 10;
#   B. on MC + SSD Tier 1: paired values with nboot = 10 vs 1000 for a
#      random subsample of Monte Carlo simulations.
#
# Usage (repo root): Rscript docs/correction/R/nboot_impact.R [n_mc_sims] [workers] [out.rds]
# Requires PSSDplusplus >= 0.4.0 and ssdtools 0.3.7.

args <- commandArgs(trailingOnly = TRUE)
k_sims <- if (length(args) >= 1) as.integer(args[1]) else 20L
workers <- if (length(args) >= 2) as.integer(args[2]) else 4L
out <- if (length(args) >= 3) args[3] else "docs/correction/results/nboot_impact.rds"

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(PSSDplusplus) })
stopifnot(as.character(utils::packageVersion("ssdtools")) == "0.3.7")

tox <- readRDS("data/input/aoc_z_tomex2.RDS") %>%
  ssd_eligible() %>% mutate(environment = env_f) %>% drop_na(effect.metric) %>%
  filter(!environment %in% c("Terrestrial", "Not Reported"), Group != "Bacterium", Group != "Plant",
         effect.metric != "HONEC", tier_zero_tech_f == "Red Criteria Passed",
         tier_zero_risk_f == "Red Criteria Passed", risk.13 != 0)

# Tier 1/2 inputs exactly as process_environment_data() builds them
tier12_inputs <- function(al, env, upper_tt = 88) {
  tt <- al %>% mutate(dose_new = particles.mL.ox.stress / (af.time * af.noec) * 1000) %>% drop_na(dose_new) %>%
    filter(dose_new > 0, between(size.length.um.used.for.conversions, 1, upper_tt),
           translocatable != "not translocatable", shape_f != "Not Reported", poly_f != "Not Reported",
           environment %in% env, Group != "Bacterium", Group != "Plant", effect.metric != "HONEC")
  fd <- al %>% filter(Group != "Algae") %>%
    mutate(dose_new = particles.mL.food.dilution / (af.time * af.noec) * 1000) %>% drop_na(dose_new) %>%
    filter(dose_new > 0, between(size.length.um.used.for.conversions, 1, 5000), ingestible != "not ingestible",
           poly_f != "Not Reported", shape_f != "Not Reported", environment %in% env,
           Group != "Bacterium", Group != "Plant", effect.metric != "HONEC")
  list(TT = tt, FD = fd)
}
t1 <- function(d, nboot, seed) PSSDplusplus:::SSD_function_t1(d, hcxlcl = 5, nboot = nboot, seed = seed)$hc5lcl

## A. deterministic SSD (reused from an earlier run of this script if present)
al <- suppressWarnings(align_data(tox, verbose = FALSE))
prev <- if (file.exists(out)) readRDS(out) else NULL
if (!is.null(prev$deterministic)) {
  det <- prev$deterministic
  seed_spread <- prev$seed_spread
} else {
det <- bind_rows(lapply(c("Freshwater", "Marine"), function(env) {
  bind_rows(lapply(c(10, 1000), function(nb) {
    PSSDplusplus:::process_environment_data(al, env, upper.tissue.trans.size.um = 88, nboot = nb, seed = 99) %>%
      select(Tier, contains("Default"), contains("5th"), contains("95th")) %>%
      pivot_longer(-Tier, names_to = "metric", values_to = "value") %>%
      mutate(environment = env, nboot = nb)
  }))
}))
seed_spread <- bind_rows(lapply(c("Freshwater", "Marine"), function(env) {
  inp <- tier12_inputs(al, env)
  bind_rows(lapply(names(inp), function(erm) tibble(
    environment = env, ERM = erm, seed = c(99, 1:19),
    t1_nboot10 = sapply(c(99, 1:19), function(s) t1(inp[[erm]], 10, s)))))
}))
}

## B. MC + SSD Tier 1 on a subsample of simulations ----------------------------
set.seed(20261020)
pm <- matrix_function(n_sobol = ceiling(k_sims / 16) + 1, params = param_default_values,
                      upper.tissue.truncation.limit = 500, x1M_set = 1, x2D_set = 5000)[seq_len(k_sims), ]
one_sim <- function(i) {
  p <- pm[i, ]
  a <- suppressWarnings(PSSDplusplus::align_data(tox, verbose = FALSE,
    upper.tissue.trans.size.um_input = p$upper.tissue.trans.size.um,
    beta_log10_body_length_input = p$sim.beta.log10.body.length,
    body_length_intercept_input = p$sim.body.length.intercept,
    R.ave.marine_input = p$R.ave.water.marine, H_W_ratio.marine_input = p$H_W_ratio.marine,
    alpha.marine_input = p$alpha.marine, a.sa.marine_input = p$a.sa.marine,
    a.v.marine_input = p$a.v.marine, a.m.marine_input = p$a.m.marine,
    R.ave.freshwater_input = p$R.ave.water.freshwater, H_W_ratio.freshwater_input = p$H_W_ratio.freshwater,
    alpha.freshwater_input = p$alpha.freshwater, a.sa.freshwater_input = p$a.sa.freshwater,
    a.v.freshwater_input = p$a.v.freshwater, a.m.freshwater_input = p$a.m.freshwater))
  bind_rows(lapply(c("Freshwater", "Marine"), function(env) {
    inp <- tier12_inputs(a, env, upper_tt = p$upper.tissue.trans.size.um)
    bind_rows(lapply(names(inp), function(erm) tibble(
      sim = i, environment = env, ERM = erm,
      t1_nboot10 = t1(inp[[erm]], 10, 99), t1_nboot1000 = t1(inp[[erm]], 1000, 99))))
  }))
}
cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("tox", "pm", "tier12_inputs", "t1"))
parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(PSSDplusplus) }))
mc <- bind_rows(parallel::parLapply(cl, seq_len(k_sims), one_sim))

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(deterministic = det, seed_spread = seed_spread, mc = mc, k_sims = k_sims, date = Sys.time()), out)
message("saved ", out)
