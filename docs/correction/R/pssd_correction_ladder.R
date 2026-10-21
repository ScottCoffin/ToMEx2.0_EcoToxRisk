# Stepwise attribution of the corrections to the PSSD++ thresholds
# (Coffin et al. 2026, J. Hazard. Mater. 503:141021, Table 2).
#
# Starting from the published configuration, each step changes one thing and
# re-runs the Monte Carlo alignment and PSSD++ with the same parameter matrix
# and seeds, so the change in the thresholds is attributable to that step.
#
#   S0 As published: published input data (commit 0d44be7), missing mouth size
#      = 5 mm, in vitro records included, body-length coefficients fixed,
#      particle shape frozen at default ratios
#   S1 + no 5 mm default for missing mouth sizes            (code error)
#   S2 + in vitro records excluded                          (code error)
#   S3 + reported mouth sizes restored (commit 3bcc930)     (data-preparation error)
#   S4 + ToMEx 2.0 database v2.0.1 and record-level body-length / in vivo
#      corrections (current data/input)                     (data refinement)
#   S5 + body-length model uncertainty propagated to predicted mouth sizes
#                                                           (code error)
#   S6 + particle-shape uncertainty propagated to tested particles
#                                                           (code error)
#
# All steps use the log-normal species distributions of the published
# analysis (rmore_method = "lognormal"; described in the SI).
#
# Usage (repo root):
#   Rscript docs/correction/R/pssd_correction_ladder.R [n_mc] [sim] [steps]
#   small check : n_mc = 30,   sim = 20   (minutes)
#   production  : n_mc = 9424, sim = 300  (~3+ h and ~60 GB per step: run
#                 one step at a time, e.g. "S0", with nothing else running)
# Each run writes docs/correction/results/ladder_nmc<n>_sim<s>[_<steps>].rds;
# the correction document combines all files with the same n_mc and sim.
#
# Checkpointing: every finished step is saved to
# docs/correction/results/checkpoints/nmc<n>_sim<s>/<step>.rds and skipped on a
# re-run; within a step the MC alignment (every 500 simulations) and each
# PSSD++ combination are checkpointed, so a crashed run resumes where it
# stopped. Re-run the same command to resume.
# Requires PSSDplusplus >= 0.4.0 and git (to read the published input).

args <- commandArgs(trailingOnly = TRUE)
n_mc <- if (length(args) >= 1) as.integer(args[1]) else 30L
sim <- if (length(args) >= 2) as.integer(args[2]) else 20L
only <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else NULL
out <- sprintf("docs/correction/results/ladder_nmc%d_sim%d%s.rds", n_mc, sim,
               if (is.null(only)) "" else paste0("_", paste(only, collapse = "-")))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(PSSDplusplus)
})
stopifnot(utils::packageVersion("PSSDplusplus") >= "0.4.0",
          "checkpoint_dir" %in% names(formals(MC_sim_align_parallel)))
root <- normalizePath(".")
stopifnot(file.exists(file.path(root, "ToMEx2.0_EcoToxRisk.Rproj")))

read_input <- function(commit = NULL) {
  if (is.null(commit)) return(readRDS(file.path(root, "data/input/aoc_z_tomex2.RDS")))
  f <- tempfile(fileext = ".RDS")
  system2("git", c("-C", shQuote(root), "show", paste0(commit, ":data/input/aoc_z_tomex2.RDS")), stdout = f)
  readRDS(f)
}
data_versions <- list(
  published = read_input("0d44be7"),          # input of the published analysis
  reported_mouth = read_input("3bcc930"),     # + reported (Koelmans 2020) mouth sizes
  current = read_input()                      # + v2.0.1 + record-level corrections
)

# manuscript filters (ToMEx2_EcoTox.Rmd / EcoTox_MonteCarlo.Rmd)
prep_tox <- function(d, exclude_in_vitro) {
  if (exclude_in_vitro) d <- ssd_eligible(d)
  d %>%
    mutate(environment = env_f) %>%
    drop_na(effect.metric) %>%
    filter(!environment %in% c("Terrestrial", "Not Reported"), Group != "Bacterium",
           Group != "Plant", effect.metric != "HONEC",
           tier_zero_tech_f == "Red Criteria Passed",
           tier_zero_risk_f == "Red Criteria Passed", risk.13 != 0)
}

steps <- tribble(
  ~step, ~label,                                           ~category,               ~data,            ~missing_mouth, ~in_vitro_out, ~coef_uncert, ~shape_uncert, ~method,
  "S0", "As published",                                    "baseline",              "published",      "x2D",          FALSE,         FALSE,        FALSE,         "lognormal",
  "S1", "No 5 mm default for missing mouth sizes",         "code error",            "published",      "exclude",      FALSE,         FALSE,        FALSE,         "lognormal",
  "S2", "In vitro records excluded",                       "code error",            "published",      "exclude",      TRUE,          FALSE,        FALSE,         "lognormal",
  "S3", "Reported mouth sizes restored",                   "data-preparation error","reported_mouth", "exclude",      TRUE,          FALSE,        FALSE,         "lognormal",
  "S4", "Database v2.0.1 + record-level corrections",      "data refinement",       "current",        "exclude",      TRUE,          FALSE,        FALSE,         "lognormal",
  "S5", "Body-length model uncertainty propagated",        "code error",            "current",        "exclude",      TRUE,          TRUE,         FALSE,         "lognormal",
  "S6", "Particle-shape uncertainty propagated",           "code error",            "current",        "exclude",      TRUE,          TRUE,         TRUE,          "lognormal"
)
if (!is.null(only)) steps <- steps[steps$step %in% only, ]

set.seed(20261020)
pm_sampled <- matrix_function(n_sobol = ceiling(n_mc / 16) + 1, params = param_default_values,
                              upper.tissue.truncation.limit = 500, x1M_set = 1, x2D_set = 5000)[seq_len(n_mc), ]
pm_fixed <- pm_sampled %>% mutate(sim.beta.log10.body.length = 0.9341, sim.body.length.intercept = 1.12)

ckpt_root <- file.path("docs/correction/results/checkpoints", sprintf("nmc%d_sim%d", n_mc, sim))
dir.create(ckpt_root, recursive = TRUE, showWarnings = FALSE)
save_atomic <- function(x, f) { tmp <- paste0(f, ".tmp"); saveRDS(x, tmp); if (file.exists(f)) unlink(f); stopifnot(file.rename(tmp, f)) }

run_step <- function(s) {
  work <- file.path(ckpt_root, s$step)
  tox <- prep_tox(data_versions[[s$data]], s$in_vitro_out)
  # published MC re-aligned already-aligned data: widths/heights stay at the
  # default R.ave / H:W ratio (shape uncertainty not propagated)
  tox_in <- if (s$shape_uncert) tox else
    suppressWarnings(align_data(tox, missing_mouth = s$missing_mouth, verbose = FALSE))
  pm <- if (s$coef_uncert) pm_sampled else pm_fixed
  invisible(utils::capture.output(mc <- suppressWarnings(suppressMessages(MC_sim_align_parallel(
    tox_data = tox_in, param_matrix = pm, n_sim = n_mc, x1D_set = 1, x2D_set = 5000,
    num_cores = 1, missing_mouth = s$missing_mouth,
    checkpoint_dir = file.path(work, "mc"), checkpoint_every = 500)))))
  # results_df_sobol post-processing (EcoTox_MonteCarlo.Rmd) and PSSD++ filters
  mc <- mc %>%
    mutate(ingestible = ifelse(!is.na(ingestible_poly), as.character(ingestible_poly), as.character(ingestible)),
           translocatable = ifelse(!is.na(translocatable_poly), as.character(translocatable_poly), as.character(translocatable))) %>%
    filter(shape_f != "Not Reported", poly_f != "Not Reported", Group != "Bacterium",
           Group != "Plant", effect.metric != "HONEC")
  food <- mc %>% filter(ingestible != "not ingestible", particles_L_food_dilution > 0, Group != "Algae") %>%
    mutate(dose_new_particles_L = particles_L_food_dilution) %>% drop_na(particles_L_food_dilution)
  tissue <- mc %>% filter(translocatable != "not translocatable", particles_L_ox_stress > 0) %>%
    mutate(dose_new_particles_L = particles_L_ox_stress) %>% drop_na(particles_L_ox_stress)
  t34 <- function(d) d %>% filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
  reg <- list("Food Dilution" = list(base = food, t3_t4 = t34(food)),
              "Tissue Translocation" = list(base = tissue, t3_t4 = t34(tissue)))
  # each combination is cached as it finishes (and reused only for identical inputs)
  set.seed(7)
  invisible(utils::capture.output(p <- suppressWarnings(suppressMessages(make_all_pSSDs(
    MC_sim_df = mc, tiers = 3, environments = c("Freshwater", "Marine"),
    erms = c("Food Dilution", "Tissue Translocation"), erm_registry = reg, sim = sim,
    cv_uf = 0.5, rmore_method = s$method, quantile_type = 8, parallel = FALSE,
    base_cache_dir = file.path(work, "pssd"), base_output_path = file.path(tempdir(), "ladder_figs"),
    overwrite_cache = FALSE)))))
  summarize_PNECs(p) %>%
    transmute(step = s$step, Tier, Environment, ERM, Median, Q5, Q95,
              n_records_fd = n_distinct(food$unique_id), n_records_tt = n_distinct(tissue$unique_id))
}

t0 <- Sys.time()
res <- bind_rows(lapply(seq_len(nrow(steps)), function(i) {
  s <- steps[i, ]
  done <- file.path(ckpt_root, paste0(s$step, ".rds"))
  if (file.exists(done)) {
    message(s$step, ": reusing finished step")
    return(readRDS(done))
  }
  t1 <- Sys.time()
  r <- run_step(s)
  save_atomic(r, done)
  unlink(file.path(ckpt_root, s$step), recursive = TRUE) # intermediate checkpoints no longer needed
  message(sprintf("%s done in %.1f min", s$step, as.numeric(Sys.time() - t1, units = "mins")))
  r
})) %>% left_join(steps %>% select(step, label, category), by = "step")

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
save_atomic(list(results = res, n_mc = n_mc, sim = sim, steps = steps,
                 package_version = as.character(utils::packageVersion("PSSDplusplus")),
                 run_time_min = as.numeric(Sys.time() - t0, units = "mins"), date = Sys.time()), out)
message("saved ", out)
