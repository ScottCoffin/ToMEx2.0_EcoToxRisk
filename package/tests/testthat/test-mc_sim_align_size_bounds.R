# Regression tests for GitHub issue #5: size bounds passed to
# MC_sim_align_parallel() must reach align_data().

test_that("MC_sim_align_wrapper passes size bounds through to align_data", {
  tox <- PSSDplusplus::tomex2_mini
  # A param row with no alignment columns, so only the size bounds differ from
  # align_data() defaults.
  params <- data.frame(simulation_id = "sim1")

  for (bounds in list(
    list(x1M = 1, x1D = 1, x2D = 500),
    list(x1M = 5, x1D = 5, x2D = 5000)
  )) {
    wrapped <- suppressWarnings(MC_sim_align_wrapper(
      tox_data = tox,
      params = params,
      simulation_id = "sim1",
      x1M_set = bounds$x1M,
      x1D_set = bounds$x1D,
      x2D_set = bounds$x2D
    ))$tox_vals
    direct <- suppressWarnings(align_data(
      tox,
      x1M_set_input = bounds$x1M,
      x1D_set_input = bounds$x1D,
      x2D_set_input = bounds$x2D
    ))
    expect_equal(
      wrapped$particles_mL_food_dilution,
      direct$particles.mL.food.dilution
    )
    expect_equal(
      wrapped$particles_mL_ox_stress,
      direct$particles.mL.ox.stress
    )
    expect_equal(wrapped$ingestible, direct$ingestible)
  }
})

test_that("x2D_set only rescales counts for records ingesting below it", {
  tox <- PSSDplusplus::tomex2_mini
  wide <- suppressWarnings(align_data(tox, x2D_set_input = 5000))
  narrow <- suppressWarnings(align_data(tox, x2D_set_input = 500))
  expect_equal(wide$unique_id, narrow$unique_id)

  below <- !is.na(wide$max.size.ingest.um) &
    wide$max.size.ingest.um < 500 &
    !is.na(wide$particles.mL.food.dilution)
  expect_gt(sum(below), 0)

  # Ingestible size limit and volume term do not depend on x2D_set here...
  expect_identical(narrow$x2M_ingest[below], wide$x2M_ingest[below])
  expect_identical(
    narrow$mu.v.poly_ingest[below],
    wide$mu.v.poly_ingest[below]
  )
  # ...so food dilution changes only by the count-range rescaling factor
  expect_equal(
    narrow$particles.mL.food.dilution[below] /
      wide$particles.mL.food.dilution[below],
    narrow$CF_bio_ingest[below] / wide$CF_bio_ingest[below]
  )

  # Records ingesting above 500 um have their volume term truncated at 500 um
  # (records with no mouth size have NA food dilution, see #6)
  above <- (is.na(wide$max.size.ingest.um) | wide$max.size.ingest.um > 500) &
    !is.na(wide$particles.mL.food.dilution)
  expect_gt(sum(above), 0)
  expect_true(all(narrow$x2M_ingest[above] == 500))
  expect_true(all(narrow$mu.v.poly_ingest[above] != wide$mu.v.poly_ingest[above]))
})

test_that("max.size.ingest.um equal to a bound is not dropped to NA", {
  tox <- PSSDplusplus::tomex2_mini
  ok <- !is.na(tox$max.size.ingest.um)
  bound <- tox$max.size.ingest.um[ok][1]
  tox$max.size.ingest.um[ok][1:2] <- bound
  aligned <- suppressWarnings(align_data(
    tox,
    x2D_set_input = bound,
    upper.tissue.trans.size.um_input = bound
  ))
  at_bound <- !is.na(aligned$max.size.ingest.um) &
    aligned$max.size.ingest.um == bound
  expect_gt(sum(at_bound), 0)
  expect_true(all(aligned$x2M_ingest[at_bound] == bound))
  expect_true(all(aligned$x2M_trans[at_bound] == bound))
})

test_that("MC_sim_align_parallel honours x2D_set", {
  skip_if_not_installed("sensobol")
  tox <- PSSDplusplus::tomex2_mini
  set.seed(1)
  pm <- matrix_function(
    n_sobol = 2,
    params = param_default_values,
    x1M_set = 1,
    x2D_set = 5000
  )
  run <- function(...) {
    set.seed(2)
    utils::capture.output(out <- suppressMessages(MC_sim_align_parallel(
      tox_data = tox,
      param_matrix = pm,
      n_sim = 1,
      num_cores = 1,
      ...
    )))
    out
  }

  default_run <- run()
  wide <- run(x1D_set = 1, x2D_set = 5000)
  narrow <- run(x1D_set = 1, x2D_set = 500)

  # Defaults are unchanged: 1-5000 um
  expect_equal(default_run, wide)

  # Narrowing the upper bound must change the food-dilution alignment
  expect_false(isTRUE(all.equal(
    wide$particles_L_food_dilution,
    narrow$particles_L_food_dilution
  )))

  # Ingestible upper limit (x2M_ingest) is capped at x2D_set, so records with
  # max.size.ingest.um > 500 or NA (mapped to x2D_set) are affected
  aligned_narrow <- suppressWarnings(align_data(
    tox,
    x2D_set_input = 500
  ))
  expect_true(all(aligned_narrow$x2M_ingest <= 500, na.rm = TRUE))
})
