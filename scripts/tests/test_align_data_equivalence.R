# Equivalence of the script copy of align_data() (scripts/utils/functions.R)
# with the package version (PSSDplusplus::align_data()) for surface water.
# Run from the repo root:
#   Rscript -e "testthat::test_file('scripts/tests/test_align_data_equivalence.R')"

library(testthat)
# testthat runs from this file's directory; the repo root is two levels up
root <- if (file.exists("ToMEx2.0_EcoToxRisk.Rproj")) "." else normalizePath("../..")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  pkgload::load_all(file.path(root, "package"), quiet = TRUE, export_all = FALSE)
})

script_env <- new.env()
suppressMessages(suppressWarnings(
  sys.source(file.path(root, "scripts/utils/functions.R"), envir = script_env, chdir = FALSE)
))

aoc_z <- readRDS(file.path(root, "data/input/aoc_z_tomex2.RDS")) %>%
  filter(env_f %in% c("Freshwater", "Marine"),
         !grepl("sediment", exposure.route, ignore.case = TRUE))

compare_cols <- c(
  "x2M_ingest", "x2M_trans", "ingestible", "translocatable",
  "particles.mL.food.dilution", "particles.mL.ox.stress"
)

align_both <- function(df) {
  list(
    script = suppressWarnings(script_env$align_data(df)),
    package = suppressWarnings(PSSDplusplus::align_data(df, verbose = FALSE))
  )
}

test_that("script and package align_data() agree on the publication data", {
  out <- align_both(aoc_z)
  key <- intersect(out$script$rowid, out$package$rowid)
  expect_gt(length(key), 1000)
  s <- out$script[match(key, out$script$rowid), ]
  p <- out$package[match(key, out$package$rowid), ]
  for (col in compare_cols) {
    expect_equal(as.vector(s[[col]]), as.vector(p[[col]]), label = col)
  }
})

test_that("missing mouth size and body length: excluded from food dilution in both", {
  df <- aoc_z
  i <- which(!is.na(df$max.size.ingest.um) & df$Group != "Algae")[1:20]
  df$max.size.ingest.um[i] <- NA
  df$max.size.ingest.mm[i] <- NA
  df$body.length.cm[i] <- NA
  out <- align_both(df)
  for (o in out) {
    rows <- o$rowid %in% df$rowid[i]
    expect_true(all(is.na(o$x2M_ingest[rows])))
    expect_true(all(o$ingestible[rows] == "unknown"))
    expect_true(all(is.na(o$particles.mL.food.dilution[rows])))
    expect_true(all(o$x2M_trans[rows] == o$upper.tissue.trans.size.um[rows]))
  }
})

test_that("missing mouth size with a body length: predicted in both", {
  df <- aoc_z
  i <- which(!is.na(df$max.size.ingest.um) & df$Group != "Algae")[1:20]
  df$max.size.ingest.um[i] <- NA
  df$body.length.cm[i] <- 0.1
  out <- align_both(df)
  expected <- 1000 * 10^(0.9341 * log10(1) - 1.12)
  for (o in out) {
    rows <- o$rowid %in% df$rowid[i]
    x2m <- o$x2M_ingest[rows]
    expect_equal(x2m, rep(expected, length(x2m)))
  }
})

test_that("sampled body-length coefficients: script and package agree", {
  s <- suppressWarnings(script_env$align_data(
    aoc_z, beta_log10_body_length_input = 0.8, body_length_intercept_input = 1.5
  ))
  p <- suppressWarnings(PSSDplusplus::align_data(
    aoc_z, beta_log10_body_length_input = 0.8, body_length_intercept_input = 1.5,
    verbose = FALSE
  ))
  key <- intersect(s$rowid, p$rowid)
  s <- s[match(key, s$rowid), ]
  p <- p[match(key, p$rowid), ]
  expect_equal(s$max.size.ingest.um, p$max.size.ingest.um.used)
  for (col in compare_cols) {
    expect_equal(as.vector(s[[col]]), as.vector(p[[col]]), label = col)
  }
  # database predictions move with the coefficients; reported sizes do not
  src <- aoc_z[match(key, aoc_z$rowid), ]
  pred <- src$max.size.ingest.basis == "predicted (Jams et al. 2020)" &
    s$max.size.ingest.um < 5000
  rep <- startsWith(src$max.size.ingest.basis, "reported")
  expect_true(all(abs(s$max.size.ingest.um[pred] / src$max.size.ingest.um[pred] - 1) > 1e-6))
  expect_identical(s$max.size.ingest.um[rep], src$max.size.ingest.um[rep])
})

test_that("script MC wrappers read the matrix_function() column names", {
  mat_file <- file.path(root, "scripts/monte carlo/output/mat.rds")
  skip_if_not(file.exists(mat_file))
  pm <- readRDS(mat_file)
  expect_true(all(c("sim.beta.log10.body.length", "sim.body.length.intercept") %in% names(pm)))
  for (f in c("MC_sim_align_wrapper", "model_wrapper_sobol_parallel")) {
    body <- paste(deparse(script_env[[f]]), collapse = "
")
    expect_match(body, "sim.beta.log10.body.length", fixed = TRUE)
    expect_match(body, "sim.body.length.intercept", fixed = TRUE)
  }
})
