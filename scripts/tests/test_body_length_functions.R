# Tests for body_length_functions.R (body length resolution, mouth size,
# Jams et al. 2020 model domain). Run from the repository root:
#   Rscript -e "testthat::test_file('tests/test_body_length_functions.R')"
# (in ToMEx2.0_EcoToxRisk: scripts/tests/test_body_length_functions.R)

library(testthat)
suppressPackageStartupMessages(library(dplyr))

find_file <- function(candidates) {
  for (f in candidates) if (file.exists(f)) return(f)
  stop("not found: ", paste(candidates, collapse = ", "))
}
fn_file <- find_file(c("body_length_functions.R", "../body_length_functions.R",
                       "scripts/utils/body_length_functions.R",
                       "../utils/body_length_functions.R", "../../scripts/utils/body_length_functions.R"))
source(fn_file)
measured_file <- find_file(c("body_length_measured.csv", "../body_length_measured.csv",
                             "scripts/monte carlo/ref data/body_length_measured.csv",
                             "../monte carlo/ref data/body_length_measured.csv",
                             "../../scripts/monte carlo/ref data/body_length_measured.csv"))

mm <- function(cm) 10^(0.9341 * log10(cm * 10) - 1.12)

test_that("jams_domain() uses the 9 mm and 30 mm cut-offs", {
  expect_equal(jams_domain(c(5, 3, 2.99, 0.9, 0.89, 0.0624, NA)),
               c("within", "within", "less_robust", "less_robust",
                 "outside_calibration", "outside_calibration", NA))
  # any body length below 30 mm is never "within"
  expect_false(any(jams_domain(seq(0.001, 2.999, length.out = 200)) == "within"))
})

test_that("row-level corrections match DOI (incl. alternatives), species, stage and duration", {
  measured <- read_body_length_measured(measured_file)
  df <- tibble(
    doi = c("10.1002/etc.3093", "10.1002/etc.3093", "0.1021/es404295e", "https://doi.org/10.3390/nano11030649",
            "10.1016/j.envpol.2017.10.095", "10.1016/j.envpol.2017.10.095", "10.9999/other"),
    species_f = c("Hyalella azteca", "Hyalella azteca", "Tripneustes gratilla", "Mytilus galloprovincialis",
                  "Crepidula onyx", "Crepidula onyx", "Hyalella azteca"),
    life_f = c("juvenile", "juvenile", "early", "adult", "early", "adult", "juvenile"),
    exposure.duration.d = c(10, 42, 5, 21, 95, 95, 10),
    body.length.cm = c(0.78, 0.78, 10, 15, 6, 6, 0.78),
    body.size.source = "x")
  out <- suppressWarnings(apply_body_length_measured(df, measured))
  expect_equal(out$body.length.cm, c(0.340, 0.564, 0.0478, 5.0, 0.0624, 6, 0.78))
  expect_equal(out$body.length.cm.database, df$body.length.cm)
  expect_equal(out$body.length.tier, c("B", "B", "B", "A", "A", NA, NA))
  expect_true(is.na(out$body.length.basis[6]))  # Crepidula correction is for early rows only
})

test_that("body length basis is classified and adult values on early/juvenile rows are flagged", {
  df <- tibble(
    doi = "d", species_f = c("A", "A", "A", "B", "C", "D", "E"),
    life_f = c("adult", "early", "early", "juvenile", "adult", "early", "early"),
    body.length.cm = c(15, 15, 0.0249, 7.5, 1, 2, NA),
    body.size.source = c("sealifebase.org", "sealifebase.org", "https://settlement.example",
                         "sealifebase.org", "study", "reported", NA),
    body.length.basis = NA_character_)
  out <- classify_body_length_basis(df,
                                    record_length = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
                                    tomex1 = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(out$body.length.basis,
               c("adult species reference", "adult reference (stage mismatch)",
                 "literature reference", "adult reference (stage mismatch)",
                 "study-reported (ToMEx 1.0)", "ToMEx 2.0 template", "unresolved"))
  expect_equal(out$jams_domain, jams_domain(df$body.length.cm))
})

test_that("mouth size: reported record > reported species/stage > override > prediction > unresolved", {
  df <- tibble(
    species_f = c("A", "A", "A", "B", "C", "D"),
    life_f = c("adult", "adult", "early", "adult", "adult", "adult"),
    body.length.cm = c(2, 2, 0.5, 1, 1, NA),
    max.size.ingest.mm = c(0.4, NA, NA, NA, NA, NA),
    mouth.size.reported = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  overrides <- tibble(species_f = "C", life_f = "Adult", max.size.ingest.override.um = 1000)
  out <- assign_max_size_ingest(df, overrides = overrides)
  expect_equal(out$max.size.ingest.mm, c(0.4, 0.4, mm(0.5), mm(1), 1, NA))
  expect_equal(out$max.size.ingest.um, 1000 * out$max.size.ingest.mm)
  expect_equal(out$max.size.ingest.basis,
               c("reported (record)", "reported (species, life stage)",
                 "predicted (Jams et al. 2020)", "predicted (Jams et al. 2020)",
                 "override", "unresolved"))
  expect_false(any(c("mouth.size.reported", "max.size.ingest.override.um") %in% names(out)))
})

test_that("unresolved sizes are reported, never defaulted, and fail in strict mode", {
  df <- tibble(doi = c("d1", "d2", "d3"), species_f = c("A", "B", "Alga"),
               life_f = "adult", org_f = c("Fish", "Crustacea", "Algae"),
               max.size.ingest.mm = c(NA, 0.2, NA))
  path <- tempfile(fileext = ".csv")
  expect_warning(rep <- report_unresolved_body_length(df, path = path), "1 test rows")
  expect_equal(rep$species_f, "A")  # algae are not applicable for food dilution
  expect_true(file.exists(path))
  expect_error(report_unresolved_body_length(df, strict = TRUE), "no resolved body length")
})

test_that("predicted sizes for animals < 30 mm trigger a warning", {
  df <- tibble(species_f = c("A", "B", "C"), jams_domain = c("within", "less_robust", "outside_calibration"),
               max.size.ingest.basis = c("predicted (Jams et al. 2020)", "predicted (Jams et al. 2020)",
                                         "reported (record)"))
  expect_warning(f <- warn_jams_domain(df), "1 rows \\(1 species\\)")
  expect_equal(f, c(FALSE, TRUE, FALSE))
})
