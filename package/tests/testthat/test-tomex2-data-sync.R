test_that("bundled tomex2 matches the shipped source RDS", {
  source_rds <- system.file(
    "extdata",
    "aoc_z_tomex2.RDS",
    package = "PSSDplusplus",
    mustWork = TRUE
  )
  source_data <- readRDS(source_rds)

  data("tomex2", package = "PSSDplusplus", envir = environment())

  expect_identical(tomex2, source_data)
})
