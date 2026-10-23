# Regression tests for GitHub issue #6: missing mouth size must not default to
# x2D_set for Food Dilution, and the body-length model must be used.

mouth_mm <- function(body_length_cm, beta = 0.9341, intercept = 1.12) {
  10^(beta * log10(body_length_cm * 10) - intercept)
}

test_that("records with no mouth size and no body length never enter Food Dilution", {
  tox <- PSSDplusplus::tomex2_mini
  no_mouth <- is.na(tox$max.size.ingest.um) & is.na(tox$body.length.cm)
  expect_gt(sum(no_mouth & tox$Group != "Algae"), 0)
  # treat them as in vivo records (in tomex2_mini they are the in vitro trout
  # study, which the message does not count as missing a mouth size)
  tox$vivo_f[no_mouth] <- "In Vivo"

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
  i <- which(is.na(tox$max.size.ingest.um) & tox$Group != "Algae" &
    tox$size.length.um.used.for.conversions >= 1)[1]
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
  rep <- !is.na(aligned$max.size.ingest.um) &
    !startsWith(dplyr::coalesce(aligned$max.size.ingest.basis, ""), "predicted")
  expect_gt(sum(rep), 0)
  expect_identical(
    aligned$max.size.ingest.um.used[rep],
    aligned$max.size.ingest.um[rep]
  )
  row <- aligned[aligned$unique_id == i, ]
  expect_equal(row$x2M_ingest, tox$max.size.ingest.um[i])
})

test_that("sampled coefficients apply only to records without a reported mouth size", {
  tox <- PSSDplusplus::tomex2_mini
  i <- which(is.na(tox$max.size.ingest.um) & tox$Group != "Algae" &
    tox$size.length.um.used.for.conversions >= 1)[1:2]
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

test_that("database-predicted mouth sizes are re-predicted, so sampled coefficients propagate", {
  data("tomex2", package = "PSSDplusplus", envir = environment())
  base <- suppressWarnings(align_data(tomex2, verbose = FALSE))
  shifted <- suppressWarnings(align_data(
    tomex2,
    beta_log10_body_length_input = 0.8,
    body_length_intercept_input = 1.5,
    verbose = FALSE
  ))
  expect_identical(shifted$rowid, base$rowid)
  src <- tomex2[match(base$rowid, tomex2$rowid), ]
  pred <- src$max.size.ingest.basis == "predicted (Jams et al. 2020)"
  rep <- startsWith(src$max.size.ingest.basis, "reported")
  expect_gt(sum(pred), 1000)
  expect_gt(sum(rep), 1000)

  # at the default coefficients nothing changes
  has <- !is.na(src$max.size.ingest.um)
  expect_equal(base$max.size.ingest.um.used[has], src$max.size.ingest.um[has])

  # with other coefficients, only database predictions move
  expect_equal(
    shifted$max.size.ingest.um.used[pred],
    mouth_mm(src$body.length.cm[pred], beta = 0.8, intercept = 1.5) * 1000
  )
  expect_true(all(
    shifted$max.size.ingest.basis.used[pred] == "predicted (Jams et al. 2020)"
  ))
  expect_identical(
    shifted$max.size.ingest.um.used[rep],
    src$max.size.ingest.um[rep]
  )
  expect_identical(
    shifted$max.size.ingest.basis.used[rep],
    src$max.size.ingest.basis[rep]
  )
})

test_that("reported mouth sizes are never perturbed; 'body_length' is deprecated", {
  data("tomex2", package = "PSSDplusplus", envir = environment())
  rep_src <- startsWith(tomex2$max.size.ingest.basis, "reported")
  expect_warning(
    bl <- align_data(tomex2, mouth_size_source = "body_length",
                     beta_log10_body_length_input = 0.8, verbose = FALSE),
    "deprecated"
  )
  ref <- suppressWarnings(align_data(tomex2, beta_log10_body_length_input = 0.8, verbose = FALSE))
  expect_identical(bl$max.size.ingest.um.used, ref$max.size.ingest.um.used)
  src <- tomex2[match(bl$rowid, tomex2$rowid), ]
  rep <- startsWith(src$max.size.ingest.basis, "reported")
  expect_gt(sum(rep), 1000)
  expect_identical(bl$max.size.ingest.um.used[rep], src$max.size.ingest.um[rep])
  expect_error(align_data(tomex2, mouth_size_source = "gape"))
})

test_that("MC_sim_align_parallel() propagates body-length coefficients per simulation", {
  skip_if_not_installed("sensobol") # matrix_function()
  tox <- PSSDplusplus::tomex2_mini
  set.seed(1)
  pm <- matrix_function(n_sobol = 3, params = param_default_values)
  # identical simulations except for the body-length slope
  pm <- pm[rep(1, 3), ]
  pm$simulation_id <- paste0("sim", 1:3)
  pm$sim.beta.log10.body.length <- c(0.8, 0.9341, 1.05)
  pm$sim.body.length.intercept <- 1.12
  out <- suppressWarnings(suppressMessages(utils::capture.output(
    mc <- MC_sim_align_parallel(tox, pm, n_sim = 3, num_cores = 1)
  )))
  mc <- mc[!is.na(mc$particles_L_food_dilution), ]
  n_values <- tapply(
    mc$particles_L_food_dilution, as.character(mc$unique_id),
    function(x) length(unique(signif(x, 8)))
  )
  # unique_id is the row number of the input data
  basis <- tox$max.size.ingest.basis[as.integer(names(n_values))]
  pred <- basis == "predicted (Jams et al. 2020)"
  rep <- startsWith(basis, "reported")
  expect_gt(sum(pred), 0)
  expect_gt(sum(rep), 0)
  # predicted mouth sizes move with the slope (unless capped at x2D_set in
  # every simulation); reported ones never do
  uncapped <- tox$body.length.cm[as.integer(names(n_values))] * 10 < 500
  expect_gt(sum(pred & uncapped), 0)
  expect_true(all(n_values[pred & uncapped] > 1))
  expect_true(all(n_values[rep] == 1))
  # the deprecated value is accepted (with a warning) and changes nothing
  expect_warning(utils::capture.output(suppressMessages(
    MC_sim_align_parallel(tox, pm, n_sim = 1, num_cores = 1,
                          mouth_size_source = "body_length")
  )), "deprecated")
})

test_that("matrix_function samples the intercept with its own standard error", {
  skip_if_not_installed("sensobol")
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

test_that("jams_domain() flags every body length below 30 mm", {
  expect_equal(
    jams_domain(c(3, 2.99, 0.9, 0.89, NA)),
    c("within", "less_robust", "less_robust", "outside_calibration", NA)
  )
  expect_false(any(jams_domain(seq(0.001, 2.999, length.out = 100)) == "within"))
})

test_that("an override is used before the body-length prediction, never over a reported size", {
  tox <- PSSDplusplus::tomex2_mini
  i <- which(is.na(tox$max.size.ingest.um) & tox$Group != "Algae" &
    tox$size.length.um.used.for.conversions >= 1)[1:2]
  j <- which(!is.na(tox$max.size.ingest.um) & tox$Group != "Algae" &
    tox$size.length.um.used.for.conversions >= 1)[1]
  tox$body.length.cm[i] <- 5
  tox$max.size.ingest.override.um <- NA_real_
  tox$max.size.ingest.override.um[c(i[1], j)] <- 1000
  aligned <- suppressWarnings(align_data(tox, verbose = FALSE))
  r1 <- aligned[aligned$unique_id == i[1], ]
  r2 <- aligned[aligned$unique_id == i[2], ]
  rj <- aligned[aligned$unique_id == j, ]
  expect_equal(r1$max.size.ingest.um.used, 1000)
  expect_equal(r1$max.size.ingest.basis.used, "override")
  expect_equal(r2$max.size.ingest.um.used, mouth_mm(5) * 1000)
  expect_equal(r2$max.size.ingest.basis.used, "predicted (Jams et al. 2020)")
  expect_equal(rj$max.size.ingest.um.used, tox$max.size.ingest.um[j])
})

test_that("predicted mouth sizes for animals < 30 mm raise a warning", {
  tox <- PSSDplusplus::tomex2_mini
  i <- which(is.na(tox$max.size.ingest.um) & tox$Group != "Algae" &
    tox$size.length.um.used.for.conversions >= 1)[1]
  tox$body.length.cm[i] <- 0.5
  expect_warning(
    aligned <- align_data(tox, verbose = TRUE),
    "body-length model for animals < 30 mm",
    fixed = TRUE
  )
  expect_equal(aligned$jams_domain[aligned$unique_id == i], "outside_calibration")
})

test_that("tomex2 carries the 2026 body-length corrections and never a default size", {
  data("tomex2", package = "PSSDplusplus", envir = environment())
  len <- function(doi) unique(tomex2$body.length.cm[tomex2$doi == doi])
  expect_equal(len("10.1021/acs.est.8b00168"), 5.9)        # Pinctada, Gardon 2018
  expect_equal(len("10.3390/nano11030649"), 5.0)           # Mytilus, Capolupo 2021
  expect_equal(len("10.1007/s11356-021-12506-8"), 6.536)   # Mytilus, Abidli 2021
  expect_equal(len("10.1016/j.envpol.2017.10.095"), 0.0624) # Crepidula larvae
  expect_equal(len("10.1021/es404295e"), 0.0478)           # Tripneustes larvae
  expect_equal(sort(len("10.1002/etc.3093")), c(0.340, 0.564)) # Hyalella, Au 2015
  # Hyalella keeps its reported mouth size
  expect_equal(unique(tomex2$max.size.ingest.mm[tomex2$doi == "10.1002/etc.3093"]), 0.112)
  expect_true(all(is.na(tomex2$max.size.ingest.mm[tomex2$max.size.ingest.basis == "unresolved"])))
  expect_equal(tomex2$jams_domain, jams_domain(tomex2$body.length.cm))
  # Hydra attenuata (Murphy & Quinn 2018) is freshwater
  expect_true(all(tomex2$env_f[tomex2$Species %in% "Hydra attenuata"] == "Freshwater"))
})

test_that("in vitro records never reach an SSD", {
  expect_equal(
    normalize_vivo(c("In Vivo", "in vitro", "invivo", " In-Vitro ")),
    c("In Vivo", "In Vitro", "In Vivo", "In Vitro")
  )
  expect_error(normalize_vivo("ex vivo"), "Unrecognized")
  df <- data.frame(vivo_f = c("In Vivo", "In Vitro", NA), x = 1:3)
  expect_equal(ssd_eligible(df)$x, 1L)
  expect_error(ssd_eligible(df, strict = TRUE), "no in vitro/in vivo status")
  expect_error(ssd_eligible(data.frame(x = 1)), "no vivo_f column")

  data("tomex2", package = "PSSDplusplus", envir = environment())
  # the in vitro rainbow trout study (primary kidney cell cultures)
  fsi <- tomex2$doi %in% "10.1016/j.fsi.2021.04.014"
  expect_true(all(tomex2$vivo_f[fsi] == "In Vitro"))
  expect_true(all(tomex2$ssd_exclusion_reason[fsi] == "in vitro"))
  # Richardson et al. 2021 exposed whole larvae: corrected to in vivo
  expect_true(all(tomex2$vivo_f[tomex2$doi %in% "10.1016/j.marpolbul.2021.112369"] == "In Vivo"))

  aligned <- suppressWarnings(align_data(tomex2[fsi | tomex2$Species %in% "Daphnia magna", ],
                                         verbose = FALSE))
  aligned$particles.mL.ox.stress[is.na(aligned$particles.mL.ox.stress)] <- 1
  th <- suppressWarnings(suppressMessages(
    process_environment_data(aligned, "Freshwater", nboot = 2)
  ))
  expect_true(is.data.frame(th))
  expect_false("10.1016/j.fsi.2021.04.014" %in% ssd_eligible(tomex2)$doi)
})

test_that("in vitro records are not reported as missing a mouth size", {
  df <- data.frame(max.size.ingest.um = NA_real_, body.length.cm = NA_real_,
                   Group = "Fish", Species = c("a", "b"),
                   vivo_f = c("In Vitro", "In Vivo"))
  expect_message(n <- report_missing_mouth(df), "1 record")
  expect_equal(n, 1L)
  df$vivo_f <- "In Vitro"
  expect_silent(report_missing_mouth(df))
})

test_that("process_environment_data() reports both ERMs in particles/L for surface water", {
  data("tomex2", package = "PSSDplusplus", envir = environment())
  d <- tomex2[!is.na(tomex2$effect.metric) & tomex2$tier_zero_tech_f == "Red Criteria Passed" &
                tomex2$tier_zero_risk_f == "Red Criteria Passed" & tomex2$risk.13 != 0 &
                tomex2$env_f %in% "Freshwater", ]
  al <- suppressWarnings(align_data(d, verbose = FALSE))
  set.seed(1)
  th <- suppressWarnings(suppressMessages(process_environment_data(al, "Freshwater", nboot = 11)))
  tt <- th$`Tissue Translocation (Default)`[th$Tier == "Tier2"]
  # SSD input in particles/L (particles/mL x 1000), so the HC5 lies within its range
  x <- ssd_eligible(al)
  x <- x$particles.mL.ox.stress / (x$af.time * x$af.noec) * 1000
  x <- x[is.finite(x) & x > 0]
  expect_true(is.finite(tt))
  expect_true(tt >= min(x) && tt <= max(x))
  expect_gt(tt, 100) # about 900 particles/L; was about 0.9 when left in particles/mL
})

test_that("SSD fits are seeded like the manuscript and leave the caller's RNG alone", {
  set.seed(2024)
  d <- data.frame(
    Species = paste0("sp", 1:8), Group = rep(c("Fish", "Crustacea"), 4),
    dose_new = c(3, 12, 40, 95, 210, 600, 1500, 4000),
    risk.13 = 2, bio_f = "Organism"
  )
  before <- .Random.seed
  a <- make_tiered_SSDs(d, tier = 1)
  expect_identical(.Random.seed, before)
  b <- make_tiered_SSDs(d, tier = 1)
  expect_identical(a$hc5lcl, b$hc5lcl)
  expect_identical(
    SSD_function_t1(d, hcxlcl = 5)$hc5lcl,
    { set.seed(99); dists <- ssdtools::ssd_fit_dists(
        d |> dplyr::group_by(Species, Group) |>
          dplyr::summarize(Conc = stats::quantile(dose_new, 0.25), .groups = "drop"),
        left = "Conc", dists = c("llogis", "lnorm", "gamma", "lgumbel"),
        computable = FALSE, silent = FALSE)
      stats::predict(dists, average = TRUE, nboot = 10, ci = TRUE)$lcl[5] }
  )
})
