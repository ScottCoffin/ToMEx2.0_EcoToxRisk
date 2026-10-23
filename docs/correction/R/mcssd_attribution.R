# Attribution of the change in the MC + SSD thresholds to the two Monte Carlo
# fixes (E4: body-length coefficient uncertainty; E5: particle-shape
# uncertainty for tested particles), on a random subsample of simulations.
# 2 x 2 design on the current data: coefficients fixed / sampled x shape
# frozen (re-aligned already-aligned data, as published) / propagated.
#
# Usage (repo root): Rscript docs/correction/R/mcssd_attribution.R [n_sims] [workers] [out.rds]
# Requires PSSDplusplus >= 0.4.0 and ssdtools 0.3.7.

args <- commandArgs(trailingOnly = TRUE)
k_sims <- if (length(args) >= 1) as.integer(args[1]) else 100L
workers <- if (length(args) >= 2) as.integer(args[2]) else 2L
out <- if (length(args) >= 3) args[3] else "docs/correction/results/mcssd_attribution.rds"

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(PSSDplusplus) })
stopifnot(as.character(utils::packageVersion("ssdtools")) == "0.3.7")

tox <- readRDS("data/input/aoc_z_tomex2.RDS") %>%
  ssd_eligible() %>% mutate(environment = env_f) %>% drop_na(effect.metric) %>%
  filter(!environment %in% c("Terrestrial", "Not Reported"), Group != "Bacterium", Group != "Plant",
         effect.metric != "HONEC", tier_zero_tech_f == "Red Criteria Passed",
         tier_zero_risk_f == "Red Criteria Passed", risk.13 != 0)
tox_prealigned <- suppressWarnings(align_data(tox, verbose = FALSE)) # as the published MC did

set.seed(20261021)
pm <- matrix_function(n_sobol = ceiling(k_sims / 16) + 1, params = param_default_values,
                      upper.tissue.truncation.limit = 500, x1M_set = 1, x2D_set = 5000)[seq_len(k_sims), ]

configs <- expand.grid(coef = c("fixed", "sampled"), shape = c("frozen", "propagated"), stringsAsFactors = FALSE)

one <- function(i) {
  p <- pm[i, ]
  bind_rows(lapply(seq_len(nrow(configs)), function(j) {
    cf <- configs[j, ]
    d <- if (cf$shape == "propagated") tox else tox_prealigned
    a <- suppressWarnings(align_data(d, verbose = FALSE,
      upper.tissue.trans.size.um_input = p$upper.tissue.trans.size.um,
      beta_log10_body_length_input = if (cf$coef == "sampled") p$sim.beta.log10.body.length else 0.9341,
      body_length_intercept_input = if (cf$coef == "sampled") p$sim.body.length.intercept else 1.12,
      R.ave.marine_input = p$R.ave.water.marine, H_W_ratio.marine_input = p$H_W_ratio.marine,
      alpha.marine_input = p$alpha.marine, a.sa.marine_input = p$a.sa.marine,
      a.v.marine_input = p$a.v.marine, a.m.marine_input = p$a.m.marine,
      R.ave.freshwater_input = p$R.ave.water.freshwater, H_W_ratio.freshwater_input = p$H_W_ratio.freshwater,
      alpha.freshwater_input = p$alpha.freshwater, a.sa.freshwater_input = p$a.sa.freshwater,
      a.v.freshwater_input = p$a.v.freshwater, a.m.freshwater_input = p$a.m.freshwater))
    bind_rows(lapply(c("Freshwater", "Marine"), function(env) {
      t <- suppressWarnings(suppressMessages(PSSDplusplus:::process_environment_data(
        a, env, upper.tissue.trans.size.um = p$upper.tissue.trans.size.um, nboot = 10, seed = 99)))
      tibble(sim = i, coef = cf$coef, shape = cf$shape, environment = env, Tier = t$Tier,
             TT = t[[2]], FD = t[[3]], R_fw = p$R.ave.water.freshwater, R_m = p$R.ave.water.marine)
    }))
  }))
}

cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("tox", "tox_prealigned", "pm", "configs"))
parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(PSSDplusplus) }))
res <- bind_rows(parallel::parLapply(cl, seq_len(k_sims), one))

summ <- res %>% pivot_longer(c(TT, FD), names_to = "ERM", values_to = "value") %>%
  group_by(coef, shape, environment, ERM, Tier) %>%
  summarise(median = median(value, na.rm = TRUE), q05 = quantile(value, 0.05, na.rm = TRUE),
            q95 = quantile(value, 0.95, na.rm = TRUE), .groups = "drop")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(per_sim = res, summary = summ, k_sims = k_sims, date = Sys.time()), out)
message("saved ", out)
