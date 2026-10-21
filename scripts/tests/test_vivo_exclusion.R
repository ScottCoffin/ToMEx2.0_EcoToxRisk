# Tests for vivo_functions.R: in vitro records never reach an SSD.
# Run from the repository root:
#   Rscript -e "testthat::test_file('tests/test_vivo_exclusion.R')"
# (in ToMEx2.0_EcoToxRisk: scripts/tests/test_vivo_exclusion.R)

library(testthat)
suppressPackageStartupMessages(library(dplyr))

find_file <- function(candidates) {
  for (f in candidates) if (file.exists(f)) return(normalizePath(f))
  stop("not found: ", paste(candidates, collapse = ", "))
}
root <- dirname(find_file(c("vivo_functions.R", "../vivo_functions.R",
                            "scripts/utils/vivo_functions.R", "../utils/vivo_functions.R",
                            "../../scripts/utils/vivo_functions.R")))
source(file.path(root, "vivo_functions.R"))
source(file.path(root, "body_length_functions.R"))
fixture <- find_file(c("tests/fixtures/vivo_fixture.csv", "fixtures/vivo_fixture.csv",
                       "scripts/tests/fixtures/vivo_fixture.csv"))
corrections <- find_file(c(file.path(root, "vivo_corrections.csv"),
                           file.path(root, "../monte carlo/ref data/vivo_corrections.csv"),
                           "scripts/monte carlo/ref data/vivo_corrections.csv"))

test_that("normalize_vivo() maps every spelling and never turns in vitro into in vivo", {
  x <- c("In Vivo", "in vivo", "invivo", " In  Vivo ", "IN-VIVO",
         "In Vitro", "in vitro", "invitro", " in   vitro", "In-Vitro")
  expect_equal(normalize_vivo(x),
               c(rep("In Vivo", 5), rep("In Vitro", 5)))
  expect_equal(normalize_vivo("in vitro"), "In Vitro")   # the old recode bug
  expect_error(normalize_vivo(c("In Vivo", "ex vivo")), "Unrecognized")
  expect_warning(out <- normalize_vivo(c("In Vivo", NA), doi = c("a", "b")), "DOIs: b")
  expect_equal(out, c("In Vivo", NA))
  expect_error(normalize_vivo(c("In Vivo", ""), strict = TRUE), "no in vitro/in vivo status")
})

test_that("ssd_eligible() keeps only in vivo rows and reasons are recorded for the rest", {
  df <- tibble(doi = c("a", "b", "c", "d"), species_f = c("A", "B", "C", "D"),
               org_f = c("Fish", "Fish", "Fish", "Crustacea"),
               vivo_f = c("In Vivo", "In Vitro", NA, "In Vivo"),
               max.size.ingest.mm = c(1, 1, 1, NA))
  out <- ssd_eligible(df)
  expect_equal(out$doi, c("a", "d"))
  expect_equal(ssd_exclusion_reason(df),
               c(NA, "in vitro", "vivo status missing",
                 "food dilution only: unresolved body length or mouth size"))
  expect_error(ssd_eligible(df, strict = TRUE), "no in vitro/in vivo status")
  expect_error(ssd_eligible(select(df, -vivo_f)), "no vivo_f column")
})

test_that("the in vitro O. mykiss study (10.1016/j.fsi.2021.04.014) never reaches an SSD", {
  fx <- readr::read_csv(fixture, show_col_types = FALSE)
  fx <- suppressWarnings(fx %>% mutate(vivo_f = normalize_vivo(vivo_f, doi = doi))) %>%
    apply_vivo_corrections(read_vivo_corrections(corrections))
  ssd <- ssd_eligible(fx)
  expect_false("10.1016/j.fsi.2021.04.014" %in% ssd$doi)
  expect_false("Species missing" %in% ssd$species_f)
  # Richardson et al. 2021 exposed whole larvae: corrected to In Vivo, so it stays
  expect_true("Pseudechinus huttoni" %in% ssd$species_f)
  expect_equal(as.character(fx$vivo_f.database[fx$species_f == "Pseudechinus huttoni"]), "In Vitro")
  rep <- report_excluded_in_vitro(fx)
  expect_setequal(rep$species_f, c("Oncorhynchus mykiss", "Species missing"))
})

test_that("unresolved report labels in vitro rows as not applicable", {
  df <- tibble(doi = c("x", "y"), species_f = c("Oncorhynchus mykiss", "B"), life_f = "adult",
               org_f = "Fish", vivo_f = c("In Vitro", "In Vivo"), max.size.ingest.mm = NA_real_)
  expect_warning(rep <- report_unresolved_body_length(df), "1 test rows")
  expect_equal(rep$reason[rep$species_f == "Oncorhynchus mykiss"], "in vitro: not applicable")
})

test_that("every SSD code path calls ssd_eligible()", {
  app <- file.path(root, "app.R")
  skip_if_not(file.exists(app), "app.R not in this repository")
  src <- readLines(app, warn = FALSE)
  chains <- grep("dplyr::filter(exp_type_f %in% exp_type_c_ssd)", src, fixed = TRUE)
  expect_gt(length(chains), 0)
  for (i in chains) expect_match(src[i - 1], "ssd_eligible()", fixed = TRUE)
  mc <- file.path(root, "monte carlo", "ssd_functions.R")
  if (file.exists(mc)) expect_equal(sum(grepl("ssd_eligible()", readLines(mc), fixed = TRUE)), 2)
})

test_that("no code path recodes in vitro to in vivo or defaults a missing mouth size to 5 mm", {
  files <- list.files(root, pattern = "\\.(R|Rmd)$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("tests/", files)]
  txt <- unlist(lapply(files, readLines, warn = FALSE))
  expect_false(any(grepl('"in vitro"\\s*~\\s*"In Vivo"', txt)))
  expect_false(any(grepl("is\\.na\\(max\\.size\\.ingest\\.um\\)\\s*~\\s*x2D_set", txt)))
  # a record with no body length and no mouth size stays unresolved (NA), never x2D_set
  df <- tibble(species_f = "A", life_f = "adult", body.length.cm = NA_real_,
               max.size.ingest.mm = NA_real_, mouth.size.reported = FALSE)
  out <- assign_max_size_ingest(df)
  expect_true(is.na(out$max.size.ingest.mm))
  expect_equal(out$max.size.ingest.basis, "unresolved")
})

test_that("built data: SSD input is in vivo only; species with only in vitro records are absent", {
  built <- c(file.path(root, "aoc_z_tomex2.RDS"), file.path(root, "../../data/input/aoc_z_tomex2.RDS"))
  built <- built[file.exists(built)]
  skip_if(length(built) == 0, "built aoc_z_tomex2.RDS not available")
  z <- readRDS(built[1])
  ssd <- ssd_eligible(z)
  expect_true(all(as.character(ssd$vivo_f) == "In Vivo"))
  sp <- if ("species_f" %in% names(z)) "species_f" else "Species"
  only_vitro <- z %>% group_by(.data[[sp]]) %>%
    summarise(all_vitro = all(as.character(vivo_f) %in% "In Vitro")) %>%
    filter(all_vitro) %>% pull(1) %>% as.character()
  expect_true("Oncorhynchus mykiss" %in% as.character(z[[sp]][z$doi == "10.1016/j.fsi.2021.04.014"]))
  expect_false(any(as.character(ssd[[sp]]) %in% only_vitro))
  expect_false("10.1016/j.fsi.2021.04.014" %in% ssd$doi)
  expect_true(all(as.character(z$vivo_f[z$doi == "10.1016/j.marpolbul.2021.112369"]) == "In Vivo"))
  if ("ssd_exclusion_reason" %in% names(z)) {
    expect_true(all(z$ssd_exclusion_reason[as.character(z$vivo_f) %in% "In Vitro"] == "in vitro"))
  }
})
