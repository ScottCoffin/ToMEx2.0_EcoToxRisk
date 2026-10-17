# Regression tests for GitHub issue #6: missing mouth size must not default to
# x2D_set for Food Dilution, and the body-length model must be used.

mouth_mm <- function(body_length_cm, beta = 0.9341, intercept = 1.12) {
  10^(beta * log10(body_length_cm * 10) - intercept)
}

test_that("records with no mouth size and no body length never enter Food Dilution", {
  tox <- PSSDplusplus::tomex2_mini
  no_mouth <- is.na(tox$max.size.ingest.um) & is.na(tox$body.length.cm)
  expect_gt(sum(no_mouth & tox$Group != "Algae"), 0)

  expect_message(
    aligned <- suppressWarnings(align_data(tox)),
    "excluded from Food Dilution"
  )
  unknown <- is.na(aligned$max.size.ingest.um) & is.na(aligned$body.length.cm)
  expect_gt(sum(unknown), 0)
  expect_true(all(is.na(aligned$x2M_ingest[unknown])))
  expect_true(all(aligned$ingestible[unknown] == "unknown"))
  expect_true(all(is.na(aligned$particles.mL.food.dilution[unknown])))
  expect_true(all(is.na(aligned$particles.kg.food.dilution[unknown])))
  # tissue translocation is unaffected
  expect_true(all(
    aligned$x2M_trans[unknown] == aligned$upper.tissue.trans.size.um[unknown]
  ))

  # ...and they never reach the Food Dilution SSD
  sim <- data.frame(
    particles_L_food_dilution = aligned$particles.mL.food.dilution * 1000,
    ingestible = aligned$ingestible
  )
  fd <- sim[!is.na(sim$particles_L_food_dilution) &
    sim$ingestible != "not ingestible", ]
  expect_false(any(fd$ingestible == "unknown"))
})

test_that("missing_mouth = 'x2D' keeps the old behaviour", {
  tox <- PSSDplusplus::tomex2_mini
  aligned <- suppressWarnings(align_data(tox, missing_mouth = "x2D"))
  unknown <- is.na(aligned$max.size.ingest.um) & is.na(aligned$body.length.cm)
  expect_true(all(aligned$x2M_ingest[unknown] == aligned$x2D_set[unknown]))
  expect_false(any(aligned$ingestible == "unknown"))
  non_algae <- unknown & aligned$Group != "Algae"
  expect_true(any(!is.na(aligned$particles.mL.food.dilution[non_algae])))
})

test_that("a body length with no mouth size gets the Jams et al. (2020) mouth size", {
  tox <- PSSDplusplus::tomex2_mini
  i <- which(is.na(tox$max.size.ingest.um) & tox$Group != "Algae")[1]
  tox$body.length.cm[i] <- 0.1 # 1 mm, e.g. a small cladoceran
  aligned <- suppressWarnings(align_data(tox, verbose = FALSE))
  row <- aligned[aligned$unique_id == i, ]
  expect_equal(row$max.size.ingest.um.used, mouth_mm(0.1) * 1000)
  expect_equal(row$x2M_ingest, mouth_mm(0.1) * 1000)
  expect_false(row$ingestible == "unknown")
})

test_that("reported mouth sizes match the body-length model at default coefficients", {
  tox <- PSSDplusplus::tomex2_mini
  aligned <- suppressWarnings(align_data(tox, verbose = FALSE))
  ok <- !is.na(aligned$max.size.ingest.um)
  expect_equal(
    aligned$max.size.ingest.um.used[ok],
    aligned$max.size.ingest.um[ok],
    tolerance = 1e-6
  )
})

test_that("sampled body-length coefficients propagate to the ingestible size", {
  tox <- PSSDplusplus::tomex2_mini
  base <- suppressWarnings(align_data(tox, verbose = FALSE))
  shifted <- suppressWarnings(align_data(
    tox,
    beta_log10_body_length_input = 0.8,
    body_length_intercept_input = 1.5,
    verbose = FALSE
  ))
  ok <- !is.na(tox$body.length.cm[base$unique_id])
  expect_equal(
    shifted$max.size.ingest.um.used[ok],
    mouth_mm(base$body.length.cm[ok], beta = 0.8, intercept = 1.5) * 1000
  )
  # "reported" ignores the coefficients when a mouth size is reported
  reported <- suppressWarnings(align_data(
    tox,
    beta_log10_body_length_input = 0.8,
    body_length_intercept_input = 1.5,
    mouth_size_source = "reported",
    verbose = FALSE
  ))
  expect_equal(
    reported$max.size.ingest.um.used[ok],
    reported$max.size.ingest.um[ok]
  )
})

test_that("matrix_function samples the intercept with its own standard error", {
  set.seed(42)
  pm <- matrix_function(
    n_sobol = 2000,
    params = param_default_values,
    x1M_set = 1,
    x2D_set = 5000
  )
  skip_if_not("sim.body.length.intercept" %in% names(pm))
  expect_equal(
    stats::sd(pm$sim.body.length.intercept),
    param_default_values$se_body_length_intercept,
    tolerance = 0.1
  )
})
