# Regression tests for GitHub issue #8: particle geometry of polydisperse
# test materials (minimum width lost for capitalised shapes; fiber widths).

# one polydisperse surface-water record per environment, reshaped as needed
poly_template <- function(env = "Marine") {
  tox <- PSSDplusplus::tomex2
  i <- which(tox$polydispersity == "polydisperse" & tox$env_f == env &
    !grepl("sediment", tox$exposure.route, ignore.case = TRUE) &
    !is.na(tox$max.size.ingest.um) & !is.na(tox$dose.particles.mL.master))[1]
  tox[i, ]
}

make_poly <- function(shape, l_min, l_max, env = "Marine", mouth_um = 5000, ...) {
  d <- poly_template(env)
  d$shape_f <- factor(shape, levels = c("Fiber", "Fragment", "Sphere", "Not Reported"))
  d$size.length.min.um.used.for.conversions <- l_min
  d$size.length.max.um.used.for.conversions <- l_max
  d$size.length.um.used.for.conversions <- (l_min + l_max) / 2
  d$max.size.ingest.um <- mouth_um
  d$max.size.ingest.basis <- "reported (test)"
  # the database fills these with R.ave x length; they must not be used
  d$size.width.min.um.used.for.conversions <- 0.77 * l_min
  d$size.width.max.um.used.for.conversions <- 0.77 * l_max
  d$size.width.um.used.for.conversions <- NA_real_
  for (col in c("size.width.min.mm.measured", "size.width.min.mm.nominal",
                "size.width.max.mm.measured", "size.width.max.mm.nominal")) {
    d[[col]] <- NA_real_
  }
  extra <- list(...)
  for (nm in names(extra)) d[[nm]] <- extra[[nm]]
  d
}

aln <- function(d, ...) suppressWarnings(align_data(d, verbose = FALSE, ...))

test_that("a polydisperse sphere's minimum volume is a sphere's volume", {
  for (env in c("Marine", "Freshwater")) {
    a <- aln(make_poly("Sphere", 2, 50, env = env))
    expect_equal(a$size.width.min.um.used.for.conversions, 2)
    expect_equal(a$particle.volume.um3.min, (4 / 3) * pi * (2 / 2)^3)
  }
})

test_that("non-spherical minimum height uses the environment's H:W ratio", {
  a <- aln(make_poly("Fragment", 10, 50, env = "Marine"))
  expect_equal(a$size.width.min.um.used.for.conversions, a$R.ave * 10)
  expect_equal(a$size.height.min.um.used.for.conversions, a$H_W_ratio * a$R.ave * 10)
  expect_equal(
    a$particle.volume.um3.min,
    (4 / 3) * pi * (10 / 2) * (a$R.ave * 10 / 2) * (a$H_W_ratio * a$R.ave * 10 / 2)
  )
  # the sampled (or user-set) H:W ratio reaches the minimum particle
  b <- aln(make_poly("Fragment", 10, 50, env = "Marine"), H_W_ratio.marine_input = 0.5)
  expect_equal(b$size.height.min.um.used.for.conversions, 0.5 * b$R.ave * 10)
})

test_that("fiber widths: same diameter rule for min and max, never wider than long", {
  # no reported width: 15 um default for both ends, capped at the length
  a <- aln(make_poly("Fiber", 5, 2000, mouth_um = 500))
  expect_equal(a$size.width.min.um.used.for.conversions, 5) # capped at L_min
  expect_equal(a$size.width.max.um.ingest, 15)
  expect_lte(a$size.width.min.um.used.for.conversions, a$size.width.max.um.ingest)
  expect_lte(a$size.width.max.um.ingest, a$size.length.max.um.ingest)
  # partially ingestible: the maximum length is truncated at the mouth size and
  # the width does not come from the untruncated database value (0.77 x 2000)
  expect_equal(a$size.length.max.um.ingest, 500)
  expect_lt(a$size.width.max.um.ingest, 0.77 * 2000)

  # a reported width (single value) is the diameter at both ends
  b <- aln(make_poly("Fiber", 100, 2000, size.width.um.used.for.conversions = 20))
  expect_equal(b$size.width.min.um.used.for.conversions, 20)
  expect_equal(b$size.width.max.um.ingest, 20)
  expect_equal(b$size.width.max.um.trans, min(20, b$size.length.max.um.trans))

  # reported minimum/maximum widths (mm) win; measured over nominal
  c <- aln(make_poly("Fiber", 100, 2000,
    size.width.min.mm.measured = 0.010, size.width.min.mm.nominal = 0.050,
    size.width.max.mm.nominal = 0.030
  ))
  expect_equal(c$size.width.min.um.used.for.conversions, 10)
  expect_equal(c$size.width.max.um.ingest, 30)
})

test_that("fiber min width never exceeds max width on the ToMEx data", {
  a <- suppressWarnings(align_data(PSSDplusplus::tomex2, verbose = FALSE))
  f <- a[a$polydispersity == "polydisperse" & a$shape_f == "Fiber", ]
  ok <- !is.na(f$size.width.min.um.used.for.conversions) & !is.na(f$size.width.max.um.ingest)
  expect_gt(sum(ok), 100)
  expect_true(all(f$size.width.min.um.used.for.conversions[ok] <= f$size.width.max.um.ingest[ok]))
  expect_true(all(f$size.width.max.um.ingest[ok] <= f$size.length.max.um.ingest[ok]))
  p <- a[a$polydispersity == "polydisperse" & !is.na(a$size.length.min.um.used.for.conversions), ]
  expect_false(anyNA(p$size.width.min.um.used.for.conversions))
})

test_that("unrecognised shapes keep the supplied width, with a warning", {
  d <- make_poly("Sphere", 10, 50)
  d$shape_f <- "Film"
  d$size.width.min.um.used.for.conversions <- 3
  expect_warning(a <- align_data(d), "not recognised")
  expect_equal(a$size.width.min.um.used.for.conversions, 3)
})

test_that("monodisperse alignment is unchanged by the polydisperse geometry", {
  d <- make_poly("Fragment", 10, 50)
  d$polydispersity <- "monodisperse"
  d$size.length.um.used.for.conversions <- 20
  d$size.width.um.used.for.conversions <- NA_real_
  a <- aln(d)
  expect_equal(
    a$particle.volume.um3,
    (4 / 3) * pi * (20 / 2) * (a$R.ave * 20 / 2) * (a$H_W_ratio * a$R.ave * 20 / 2)
  )
})

## ---- Measured size distributions of the test material (issue #8, part B) ----

test_that("records without a measured distribution are aligned exactly as before", {
  d <- PSSDplusplus::tomex2_mini
  a0 <- suppressWarnings(align_data(d, verbose = FALSE))
  d$size.dist.family <- NA_character_
  d$size.dist.alpha <- NA_real_
  a1 <- suppressWarnings(align_data(d, verbose = FALSE))
  for (col in c("particles.mL.food.dilution", "particles.mL.ox.stress",
                "particles.kg.food.dilution", "particles.kg.ox.stress")) {
    expect_identical(a1[[col]], a0[[col]], label = col)
  }
  poly <- a0$polydispersity == "polydisperse"
  expect_true(all(a0$poly_size_method[poly] == "assumed_env_powerlaw"))
  expect_true(all(a0$poly_size_method[!poly] == "monodisperse"))
  expect_true(all(is.na(a0$poly_volume_factor_ingest[!poly])))
})

test_that("a measured power law gives the Koelmans et al. (2020) bioaccessible fraction", {
  d <- make_poly("Sphere", 1, 1000, mouth_um = 100,
    size.dist.family = "powerlaw", size.dist.alpha = 2.07)
  a <- aln(d)
  expect_equal(a$poly_size_method, "measured_powerlaw")
  expect_equal(a$CF_bioavailable_ingest, CFfnx(a = 2.07, x1D = 1, x2D = 100, x1M = 1, x2M = 1000))
  # default path for the same record uses the environment's alpha (2.07 marine)
  b <- aln(make_poly("Sphere", 1, 1000, mouth_um = 100))
  expect_equal(b$poly_size_method, "assumed_env_powerlaw")
  expect_equal(a$CF_bioavailable_ingest, b$CF_bioavailable_ingest)
})

test_that("with a measured log-normal distribution the mean ingestible volume matches numerical integration", {
  med <- 40; gsd <- 2; lmin <- 5; lmax <- 500; mouth <- 120
  d <- make_poly("Fragment", lmin, lmax, mouth_um = mouth,
    size.dist.family = "lognormal", size.dist.median.um = med, size.dist.gsd = gsd)
  a <- aln(d)
  expect_equal(a$poly_size_method, "measured_lognormal")
  R <- a$R.ave; HW <- a$H_W_ratio
  # independent check: fine grid in log(L), trapezoid rule
  u <- seq(log(lmin), log(mouth), length.out = 200001)
  L <- exp(u)
  w <- stats::dlnorm(L, log(med), log(gsd)) * L
  v <- (4 / 3) * pi * (L / 2) * (R * L / 2) * (HW * R * L / 2)
  trap <- function(y) sum(diff(u) * (head(y, -1) + tail(y, -1)) / 2)
  expect_equal(a$mu.v.mono.ingest, trap(v * w) / trap(w), tolerance = 1e-6)
  # bioaccessible number fraction from the log-normal CDF
  p <- function(x) stats::plnorm(x, log(med), log(gsd))
  expect_equal(a$CF_bioavailable_ingest, (p(mouth) - p(lmin)) / (p(lmax) - p(lmin)))
  expect_equal(a$dose.particles.mL.ingest, a$CF_bioavailable_ingest * a$dose.particles.mL.master)
  expect_equal(a$poly_volume_factor_ingest, a$mu.v.mono.ingest / a$particle.volume.um3.min)
})

test_that("D10/D50/D90 define the same log-normal as median and GSD", {
  med <- 40; gsd <- 2; z <- stats::qnorm(0.9)
  by_d <- aln(make_poly("Sphere", 5, 500, mouth_um = 120,
    size.dist.family = "lognormal", size.D10.um = med / gsd^z,
    size.D50.um = med, size.D90.um = med * gsd^z))
  by_m <- aln(make_poly("Sphere", 5, 500, mouth_um = 120,
    size.dist.family = "lognormal", size.dist.median.um = med, size.dist.gsd = gsd))
  expect_equal(by_d$mu.v.mono.ingest, by_m$mu.v.mono.ingest)
  expect_equal(by_d$mu.sa.mono.trans, by_m$mu.sa.mono.trans)
})

test_that("bad or incomplete size-distribution inputs are caught", {
  expect_error(aln(make_poly("Sphere", 5, 500, size.dist.family = "weibull")), "size.dist.family")
  expect_warning(
    a <- align_data(make_poly("Sphere", 5, 500, size.dist.family = "lognormal"), verbose = FALSE),
    "without usable parameters"
  )
  expect_equal(a$poly_size_method, "assumed_env_powerlaw")
})
