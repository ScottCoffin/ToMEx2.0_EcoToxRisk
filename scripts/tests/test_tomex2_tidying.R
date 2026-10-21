# End-to-end test of ToMEx2.0fxn() (scripts/utils/ToMEx2.0_Data_Tidying_functions.R):
# body length and maximum ingestible size (mouth size) follow the rules in
# body_length_functions.R (GitHub issue #6 and .claude/body_length_corrections.md),
# and the result matches data/input/aoc_z_tomex2.RDS, which is built by
# SCCWRP/ToMEx_AquaticOrganisms (ToMEx2.0_Onboarding/ToMEx2.0_Data_Tidying.R).
# Unit tests of the shared functions are in test_body_length_functions.R.
#
# Run from the repo root:
#   Rscript -e "testthat::test_file('scripts/tests/test_tomex2_tidying.R')"
#
# Needs the raw ToMEx 2.0 data, built from the validated templates in the
# ToMEx_AquaticOrganisms repo (lines 19-190 of ToMEx2.0_Data_Tidying.R) and
# saved as an RDS. Point TOMEX2_RAW_RDS at it; otherwise the test is skipped.

library(testthat)
suppressPackageStartupMessages(library(tidyverse))

root <- if (file.exists("ToMEx2.0_EcoToxRisk.Rproj")) "." else normalizePath("../..")
ref_data <- file.path(root, "scripts/monte carlo/ref data")

test_that("ToMEx2.0fxn() reproduces the published body lengths and mouth sizes", {
  raw_rds <- Sys.getenv("TOMEX2_RAW_RDS")
  skip_if(!nzchar(raw_rds) || !file.exists(raw_rds),
          "Set TOMEX2_RAW_RDS to the raw ToMEx 2.0 data (see header)")

  suppressMessages(suppressWarnings({
    source(file.path(root, "scripts/utils/functions.R"))
    source(file.path(root, "scripts/utils/ToMEx2.0_Data_Tidying_functions.R"))
  }))
  H_W_ratio <<- 0.67
  out <- suppressMessages(suppressWarnings(ToMEx2.0fxn(
    aoc_setup = readRDS(file.path(root, "data/input/aoc_setup.RDS")),
    tomex2.0 = readRDS(raw_rds),
    bodysize_addons = read_csv(file.path(ref_data, "gape_size.csv"), show_col_types = FALSE),
    R.ave.water.marine = 0.77, R.ave.water.freshwater = 0.67,
    R.ave.sediment.freshwater = 0.70, R.ave.sediment.marine = 0.75,
    beta_log10_body_length = 0.9341, body_length_intercept = 1.12,
    body_length_measured = file.path(ref_data, "body_length_measured.csv"),
    vivo_corrections = file.path(ref_data, "vivo_corrections.csv")
  )))
  pub <- readRDS(file.path(root, "data/input/aoc_z_tomex2.RDS"))
  expect_equal(nrow(out), nrow(pub))

  key <- function(d) d %>%
    transmute(source = as.character(source), doi, Species = as.character(Species),
              life_f = as.character(life_f), body.length.cm,
              mm = round(max.size.ingest.mm, 10), body.length.basis,
              max.size.ingest.basis, jams_domain, vivo_f = as.character(vivo_f),
              ssd_exclusion_reason) %>%
    count(across(everything()))
  expect_equal(arrange(key(out), across(everything())), arrange(key(pub), across(everything())))

  # reported (Koelmans et al. 2020) mouth sizes are kept
  dr <- out %>% filter(Species == "Danio rerio", life_f %in% c("adult", "early"))
  expect_equal(unique(dr$max.size.ingest.mm), 0.4)
  # row-level corrections from body_length_measured.csv
  expect_equal(unique(out$body.length.cm[out$doi == "10.1021/acs.est.8b00168"]), 5.9)
  expect_equal(unique(out$body.length.cm[out$doi == "10.1021/es404295e"]), 0.0478)
  expect_equal(sort(unique(out$body.length.cm[out$doi == "10.1002/etc.3093"])), c(0.340, 0.564))
  expect_equal(unique(out$max.size.ingest.mm[out$doi == "10.1002/etc.3093"]), 0.112) # reported wins
  # no default size anywhere: an unresolved size stays NA
  expect_true(all(is.na(out$max.size.ingest.mm[out$max.size.ingest.basis == "unresolved"])))
  expect_true(all(out$body.length.basis[is.na(out$body.length.cm)] == "unresolved"))
  # in vitro records are flagged; Richardson et al. 2021 (whole larvae) is in vivo
  expect_true(all(out$ssd_exclusion_reason[as.character(out$vivo_f) == "In Vitro"] == "in vitro"))
  expect_true(all(as.character(out$vivo_f[out$doi == "10.1016/j.marpolbul.2021.112369"]) == "In Vivo"))
  expect_true(all(as.character(out$vivo_f[out$doi == "10.1016/j.fsi.2021.04.014"]) == "In Vitro"))
  # Jams domain matches body length; < 30 mm is never "within"
  expect_false(any(out$jams_domain %in% "within" & out$body.length.cm < 3))
})
