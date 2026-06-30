test_that("hc5_from_endpoints reproduces direct do.pSSD_mod output", {
  endpoints <- data.frame(
    Species = rep(c("sp_a", "sp_b", "sp_c", "sp_d"), each = 3),
    ec = c(
      100, 140, 190,
      250, 320, 430,
      800, 900, 1200,
      2000, 2400, 3100
    ),
    ec_sd = c(
      5, 7, 9,
      12, 16, 21,
      40, 45, 60,
      100, 120, 155
    )
  )
  direct_input <- endpoints[order(endpoints$Species), , drop = FALSE]
  direct_input$.row <- stats::ave(
    seq_len(nrow(direct_input)),
    direct_input$Species,
    FUN = seq_along
  )
  DP <- reshape2::acast(direct_input, .row ~ Species, value.var = "ec", drop = FALSE)
  DP.SD <- reshape2::acast(
    direct_input,
    .row ~ Species,
    value.var = "ec_sd",
    drop = FALSE
  )

  set.seed(20260630)
  wrapped <- hc5_from_endpoints(
    endpoints,
    sim = 40,
    cv_dp = 0.15,
    cv_uf = 0.5,
    rmore_method = "step",
    hcx = 0.05,
    apply_assessment_factors = FALSE
  )

  set.seed(20260630)
  direct <- PSSDplusplus:::do.pSSD_mod(
    DP = DP,
    DP.SD = DP.SD,
    UFt = NULL,
    UFdd = NULL,
    SIM = 40,
    CV.DP = 0.15,
    CV.UF = 0.5,
    rmore_method = "step",
    apply_assessment_factors = FALSE
  )
  direct_hc <- apply(direct, 2, stats::quantile, probs = 0.05, type = 8)

  expect_equal(wrapped$pSSD, direct)
  expect_equal(wrapped$hc, direct_hc)
})

test_that("assessment factors off leaves endpoint scale unchanged", {
  endpoints <- data.frame(
    Species = rep(c("sp_a", "sp_b", "sp_c", "sp_d"), each = 3),
    ec = c(100, 120, 150, 250, 300, 360, 800, 900, 1100, 2000, 2400, 2800)
  )
  scaled <- endpoints
  scaled$ec <- scaled$ec * 10

  set.seed(417)
  base <- hc5_from_endpoints(
    endpoints,
    sim = 60,
    cv_dp = 0.01,
    cv_uf = 0.01,
    apply_assessment_factors = FALSE
  )
  set.seed(417)
  scaled_res <- hc5_from_endpoints(
    scaled,
    sim = 60,
    cv_dp = 0.01,
    cv_uf = 0.01,
    apply_assessment_factors = FALSE
  )

  expect_equal(scaled_res$pSSD, base$pSSD * 10)
  expect_equal(stats::median(scaled_res$hc), stats::median(base$hc) * 10)
})

test_that("HC output remains in particles/L linear scale", {
  endpoints <- data.frame(
    Species = rep(c("sp_a", "sp_b", "sp_c", "sp_d"), each = 3),
    ec = c(100, 120, 150, 250, 300, 360, 800, 900, 1100, 2000, 2400, 2800)
  )

  set.seed(99)
  res <- hc5_from_endpoints(
    endpoints,
    sim = 80,
    cv_dp = 0.001,
    cv_uf = 0.001,
    hcx = 0.05,
    apply_assessment_factors = FALSE
  )

  expect_named(res, c("pSSD", "hc"))
  expect_equal(dim(res$pSSD), c(4L, 80L))
  expect_true(stats::median(res$hc) >= min(endpoints$ec))
  expect_true(stats::median(res$hc) <= stats::median(endpoints$ec))
})

test_that("wrapper matches the one-iteration internal pipeline slot", {
  endpoints <- data.frame(
    Species = rep(c("sp_a", "sp_b", "sp_c", "sp_d"), each = 3),
    ec = c(100, 140, 190, 250, 320, 430, 800, 900, 1200, 2000, 2400, 3100),
    ec_sd = c(5, 7, 9, 12, 16, 21, 40, 45, 60, 100, 120, 155)
  )
  direct_input <- endpoints[order(endpoints$Species), , drop = FALSE]
  direct_input$.row <- stats::ave(
    seq_len(nrow(direct_input)),
    direct_input$Species,
    FUN = seq_along
  )
  data_matrices <- list(
    DP = reshape2::acast(direct_input, .row ~ Species, value.var = "ec", drop = FALSE),
    DP.SD = reshape2::acast(
      direct_input,
      .row ~ Species,
      value.var = "ec_sd",
      drop = FALSE
    ),
    UFt = NULL,
    UFdd = NULL
  )

  set.seed(1234)
  wrapped <- hc5_from_endpoints(
    endpoints,
    sim = 30,
    cv_dp = 0.15,
    cv_uf = 0.5,
    rmore_method = "step",
    hcx = 0.05,
    apply_assessment_factors = FALSE
  )

  set.seed(1234)
  pipeline <- PSSDplusplus:::run_pSSD_analysis(
    data_matrices = data_matrices,
    num_iterations = 1,
    sim = 30,
    cv_dp = 0.15,
    cv_uf = 0.5,
    data_name = "test",
    rmore_method = "step",
    apply_assessment_factors = FALSE,
    silent = TRUE
  )
  pipeline_hc <- apply(pipeline$pSSD, 2, stats::quantile, probs = 0.05, type = 8)

  expect_equal(wrapped$pSSD, pipeline$pSSD)
  expect_equal(wrapped$hc, pipeline_hc)
})
