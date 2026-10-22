# Long runs must survive a crash: finished work is saved and reused on re-run,
# and never reused when the inputs differ.

quiet <- function(expr) {
  utils::capture.output(res <- suppressWarnings(suppressMessages(expr)))
  res
}

mini_inputs <- function(n = 5) {
  set.seed(11)
  pm <- matrix_function(n_sobol = 1, params = param_default_values)
  pm <- pm[seq_len(n), ]
  list(tox = PSSDplusplus::tomex2_mini, pm = pm)
}

test_that("MC_sim_align_parallel() checkpoints chunks and resumes after a crash", {
  x <- mini_inputs(5)
  dir <- withr::local_tempdir()
  ref <- quiet(MC_sim_align_parallel(x$tox, x$pm, n_sim = 5, num_cores = 1))
  run <- function() {
    msgs <- NULL
    utils::capture.output(msgs <- utils::capture.output(type = "message", suppressWarnings(
      MC_sim_align_parallel(x$tox, x$pm, n_sim = 5, num_cores = 1,
                            checkpoint_dir = dir, checkpoint_every = 2))))
    msgs
  }

  # the run crashes in simulation 3 (second chunk)
  real_wrapper <- MC_sim_align_wrapper
  local({
    local_mocked_bindings(MC_sim_align_wrapper = function(tox_data, params, simulation_id, ...) {
      if (identical(simulation_id, "sim3")) stop("simulated crash")
      real_wrapper(tox_data, params, simulation_id, ...)
    })
    expect_error(suppressWarnings(suppressMessages(utils::capture.output(
      MC_sim_align_parallel(x$tox, x$pm, n_sim = 5, num_cores = 1,
                            checkpoint_dir = dir, checkpoint_every = 2)))), "simulated crash")
  })
  expect_length(list.files(dir, "^mc_align_.*[.]rds$"), 1) # simulations 1-2 kept
  expect_length(list.files(dir, "[.]tmp$"), 0) # no partial files

  # re-run: the finished chunk is reused, the rest computed, result unchanged
  msgs <- run()
  expect_equal(sum(grepl("reusing", msgs)), 1)
  expect_equal(sum(grepl("saved", msgs)), 2)
  out <- quiet(MC_sim_align_parallel(x$tox, x$pm, n_sim = 5, num_cores = 1,
                                     checkpoint_dir = dir, checkpoint_every = 2))
  expect_equal(out, ref, ignore_attr = TRUE)

  # different inputs never reuse these checkpoints
  other <- x$pm
  other$alpha.freshwater <- other$alpha.freshwater + 0.1
  utils::capture.output(msgs <- utils::capture.output(type = "message", suppressWarnings(
    MC_sim_align_parallel(x$tox, other, n_sim = 5, num_cores = 1,
                          checkpoint_dir = dir, checkpoint_every = 2))))
  expect_false(any(grepl("reusing", msgs)))
})

test_that("make_all_pSSDs() resumes from its cache and never reuses stale results", {
  x <- mini_inputs(3)
  mc <- quiet(MC_sim_align_parallel(x$tox, x$pm, n_sim = 3, num_cores = 1))
  food <- mc[mc$ingestible != "not ingestible" & !is.na(mc$particles_L_food_dilution) &
               mc$particles_L_food_dilution > 0 & mc$Group != "Algae", ]
  food$dose_new_particles_L <- food$particles_L_food_dilution
  tissue <- mc[mc$translocatable != "not translocatable" & !is.na(mc$particles_L_ox_stress) &
                 mc$particles_L_ox_stress > 0, ]
  tissue$dose_new_particles_L <- tissue$particles_L_ox_stress
  reg <- list("Food Dilution" = list(base = food, t3_t4 = food),
              "Tissue Translocation" = list(base = tissue, t3_t4 = tissue))
  cache <- withr::local_tempdir()
  figs <- withr::local_tempdir()
  run <- function(sim, seed = 1) {
    set.seed(seed)
    out <- utils::capture.output(p <- suppressWarnings(suppressMessages(make_all_pSSDs(
      MC_sim_df = mc, tiers = 3, environments = "Freshwater",
      erms = c("Food Dilution", "Tissue Translocation"), erm_registry = reg,
      sim = sim, rmore_method = "lognormal", parallel = FALSE,
      base_cache_dir = cache, base_output_path = figs))))
    list(p = p, log = out)
  }
  first <- run(2)
  files <- list.files(cache, recursive = TRUE, pattern = "\\.rds$")
  skip_if(length(files) == 0, "no PSSD++ combination ran on the mini data")

  # a second call reuses every finished combination (different seed, same result)
  second <- run(2, seed = 99)
  expect_true(any(grepl("Using cached", second$log)))
  expect_false(any(grepl("recomputing", second$log)))
  expect_equal(summarize_PNECs(second$p), summarize_PNECs(first$p))

  # different inputs (sim) must not reuse the cache
  third <- run(3)
  expect_true(any(grepl("different inputs; recomputing", third$log)))
  expect_false(any(grepl("Using cached", third$log)))
  expect_length(list.files(cache, recursive = TRUE, pattern = "\\.tmp$"), 0)
})

test_that("make_all_pSSDs() seeds each combination: order- and caller-seed independent", {
  x <- mini_inputs(3)
  mc <- quiet(MC_sim_align_parallel(x$tox, x$pm, n_sim = 3, num_cores = 1))
  mk <- function(col) {
    d <- mc[!is.na(mc[[col]]) & mc[[col]] > 0, ]
    d$dose_new_particles_L <- d[[col]]
    d
  }
  food <- mk("particles_L_food_dilution"); food <- food[food$Group != "Algae", ]
  tissue <- mk("particles_L_ox_stress")
  reg <- list("Food Dilution" = list(base = food, t3_t4 = food),
              "Tissue Translocation" = list(base = tissue, t3_t4 = tissue))
  run <- function(erms, seed, combo_seed = TRUE) {
    set.seed(seed)
    before <- .Random.seed
    utils::capture.output(p <- suppressWarnings(suppressMessages(make_all_pSSDs(
      MC_sim_df = mc, tiers = 3, environments = "Freshwater", erms = erms, erm_registry = reg,
      sim = 2, rmore_method = "lognormal", parallel = FALSE, combo_seed = combo_seed,
      base_cache_dir = withr::local_tempdir(.local_envir = parent.frame(2)),
      base_output_path = withr::local_tempdir(.local_envir = parent.frame(2))))))
    list(s = summarize_PNECs(p), rng_unchanged = identical(.Random.seed, before))
  }
  a <- run(c("Food Dilution", "Tissue Translocation"), seed = 1)
  skip_if(nrow(a$s) == 0, "no PSSD++ combination ran on the mini data")
  b <- run(c("Tissue Translocation", "Food Dilution"), seed = 999)
  key <- function(s) s[order(s$ERM, s$Tier), , drop = FALSE]
  expect_equal(key(a$s), key(b$s), ignore_attr = TRUE)
  expect_true(a$rng_unchanged)
  # combo_seed = FALSE draws from the caller's stream (results depend on the seed)
  c1 <- run(c("Food Dilution", "Tissue Translocation"), seed = 1, combo_seed = FALSE)
  c2 <- run(c("Food Dilution", "Tissue Translocation"), seed = 2, combo_seed = FALSE)
  expect_false(isTRUE(all.equal(key(c1$s), key(c2$s))))
})
