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

test_that("a reported mouth size is always used over the body-length prediction", {
  tox <- PSSDplusplus::tomex2_mini
  i <- which(!is.na(tox$max.size.ingest.um) & !is.na(tox$body.length.cm) &
    tox$Group != "Algae")[1]
  # reported gape far from what body length predicts
  tox$max.size.ingest.um[i] <- 0.2 * mouth_mm(tox$body.length.cm[i]) * 1000
  aligned <- suppressWarnings(align_data(
    tox,
    beta_log10_body_length_input = 0.8,
    body_length_intercept_input = 1.5,
    verbose = FALSE
  ))
  # every reported value is used as is, whatever the coefficients
  rep <- !is.na(aligned$max.size.ingest.um)
  expect_identical(
    aligned$max.size.ingest.um.used[rep],
    aligned$max.size.ingest.um[rep]
  )
  row <- aligned[aligned$unique_id == i, ]
  expect_equal(row$x2M_ingest, tox$max.size.ingest.um[i])
})

test_that("sampled coefficients apply only to records without a reported mouth size", {
  tox <- PSSDplusplus::tomex2_mini
  i <- which(is.na(tox$max.size.ingest.um) & tox$Group != "Algae")[1:2]
  tox$body.length.cm[i] <- c(0.1, 2)
  shifted <- suppressWarnings(align_data(
    tox,
    beta_log10_body_length_input = 0.8,
    body_length_intercept_input = 1.5,
    verbose = FALSE
  ))
  rows <- shifted[match(i, shifted$unique_id), ]
  expect_equal(
    rows$max.size.ingest.um.used,
    mouth_mm(c(0.1, 2), beta = 0.8, intercept = 1.5) * 1000
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

test_that("tomex2 keeps reported (Koelmans 2020) mouth sizes over predictions", {
  data("tomex2", package = "PSSDplusplus", envir = environment())
  reported <- data.frame(
    Species = c("Danio rerio", "Danio rerio", "Daphnia magna",
                "Daphnia magna", "Hyalella azteca", "Ceriodaphnia dubia"),
    life = c("adult", "early", "adult", "juvenile", "adult", "early"),
    mm = c(0.4, 0.4, 0.115, 0.115, 0.112, 0.0357)
  )
  for (k in seq_len(nrow(reported))) {
    rows <- as.character(tomex2$Species) == reported$Species[k] &
      as.character(tomex2$life_f) == reported$life[k]
    expect_gt(sum(rows), 0)
    expect_equal(
      unique(tomex2$max.size.ingest.mm[rows]),
      reported$mm[k],
      tolerance = 0.01
    )
    expect_equal(
      unique(tomex2$max.size.ingest.um[rows]),
      reported$mm[k] * 1000,
      tolerance = 0.01
    )
  }
})
