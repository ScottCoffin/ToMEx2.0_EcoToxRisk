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
#   S7 + polydisperse particle geometry fixed (ToMEx2.0_EcoToxRisk#8: smallest
#      particle's width lost for spheres/fibers/fragments; fiber widths, also the
#      monodisperse 15 um default capped at the length) and SSD size filters by
#      bioaccessibility labels instead of the midpoint length
#      [particle counts of mass-dosed tests as converted by the database]
#                                                           (code errors, E7, E10)
#   S8 + particle counts of mass-dosed tests from the dosed mass with the
#      alignment's own geometry and size distribution      (code error, E8)
#
# S0-S6 use the geometry of PSSDplusplus 0.4.0 and must be run with 0.4.0
# (GitHub release v0.4.0); S7-S8 need 0.4.1 or later
# (align_data(particle_count_source = "database" / "mass")). The script checks this.
#
# All steps use the log-normal species distributions of the published
# analysis (rmore_method = "lognormal"; described in the SI).
#
# Usage (repo root):
#   Rscript docs/correction/R/pssd_correction_ladder.R [n_mc] [sim] [steps]
#   small check : n_mc = 30,   sim = 20   (minutes)
#   correction  : n_mc = 1000, sim = 300 (one step per process: ~30-110 min
#                 and up to ~70 GB for the early steps; run one at a time).
#                 At n_mc = 9424 the results agree within 4% (S0-S2 checked)
#                 but S0 takes ~12 h.
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
  "S6", "Particle-shape uncertainty propagated",           "code error",            "current",        "exclude",      TRUE,          TRUE,         TRUE,          "lognormal",
  "S7", "Particle geometry and size filters fixed",        "code error",            "current",        "exclude",      TRUE,          TRUE,         TRUE,          "lognormal",
  "S8", "Mass-dosed particle counts from the dosed mass",  "code error",            "current",        "exclude",      TRUE,          TRUE,         TRUE,          "lognormal"
)
# the geometry fix (S7) is in the code, not a switch: S0-S6 need 0.4.0, S7-S8 >= 0.4.1;
# S7 keeps the database's particle counts of mass-dosed tests, S8 recomputes them
steps$package <- ifelse(steps$step %in% c("S7", "S8"), ">= 0.4.1", "0.4.0")
steps$counts <- c(S7 = "database", S8 = "mass")[steps$step] # NA for S0-S6 (argument absent in 0.4.0)
if (!is.null(only)) steps <- steps[steps$step %in% only, ]
pkg_v <- utils::packageVersion("PSSDplusplus")
for (i in seq_len(nrow(steps))) {
  done <- file.path("docs/correction/results/checkpoints", sprintf("nmc%d_sim%d", n_mc, sim), paste0(steps$step[i], ".rds"))
  ok <- if (steps$package[i] == "0.4.0") pkg_v == "0.4.0" else pkg_v >= "0.4.1"
  if (!ok && !file.exists(done)) {
    stop(sprintf("step %s needs PSSDplusplus %s (installed: %s)", steps$step[i], steps$package[i], pkg_v), call. = FALSE)
  }
}

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
    suppressWarnings(do.call(align_data, c(list(tox, missing_mouth = s$missing_mouth, verbose = FALSE),
                                            if (!is.na(s$counts)) list(particle_count_source = s$counts))))
  pm <- if (s$coef_uncert) pm_sampled else pm_fixed
  mc_args <- list(tox_data = tox_in, param_matrix = pm, n_sim = n_mc, x1D_set = 1, x2D_set = 5000,
                  num_cores = 1, missing_mouth = s$missing_mouth,
                  checkpoint_dir = file.path(work, "mc"), checkpoint_every = 500)
  if (!is.na(s$counts)) mc_args$particle_count_source <- s$counts
  invisible(utils::capture.output(mc <- suppressWarnings(suppressMessages(
    do.call(MC_sim_align_parallel, mc_args)))))
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
  # one combination at a time (memory is freed in between; a single call for
  # all four grew to >100 GB at production size), each with its own seed as in
  # ToMEx2_EcoTox.Rmd, and cached as it finishes (reused only for identical inputs)
  combos <- expand.grid(environment = c("Freshwater", "Marine"),
                        erm = c("Food Dilution", "Tissue Translocation"), stringsAsFactors = FALSE)
  pnecs <- bind_rows(lapply(seq_len(nrow(combos)), function(k) {
    combo_id <- paste0("Tier3_", combos$environment[k], "_", combos$erm[k])
    set.seed(2026 + sum(utf8ToInt(combo_id)))
    invisible(utils::capture.output(p <- suppressWarnings(suppressMessages(make_all_pSSDs(
      MC_sim_df = mc, tiers = 3, environments = combos$environment[k], erms = combos$erm[k],
      erm_registry = reg, sim = sim, cv_uf = 0.5, rmore_method = s$method, quantile_type = 8,
      parallel = FALSE, base_cache_dir = file.path(work, "pssd"),
      base_output_path = file.path(tempdir(), "ladder_figs"), overwrite_cache = FALSE)))))
    out <- summarize_PNECs(p)
    rm(p)
    invisible(gc())
    out
  }))
  pnecs %>%
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
