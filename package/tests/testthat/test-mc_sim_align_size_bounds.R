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
