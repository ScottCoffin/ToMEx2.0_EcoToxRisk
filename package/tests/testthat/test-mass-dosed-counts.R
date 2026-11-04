# Mass-dosed tests (particle counts converted from mg/L), unit conversions and
# monodisperse fiber widths (0.4.1, ScottCoffin/ToMEx2.0_EcoToxRisk#8 review).

mass_record <- function(poly = FALSE) {
  tox <- PSSDplusplus::tomex2
  i <- which(tox$env_f == "Marine" & tox$shape_f == "Fragment" &
    tox$polydispersity == (if (poly) "polydisperse" else "monodisperse") &
    as.character(tox$dose.particles.mL.master.converted.reported) == "converted" &
    !is.na(tox$dose.mg.L.master) & !is.na(tox$density.g.cm3) &
    !is.na(tox$max.size.ingest.um) &
    !grepl("sediment", tox$exposure.route, ignore.case = TRUE))[1]
  tox[i, ]
}

test_that("a mass-dosed monodisperse test keeps its volume dose = mass / density", {
  d <- mass_record()
  for (R in c(0.5, 0.77, 1)) {
    a <- suppressWarnings(align_data(d, verbose = FALSE, R.ave.marine_input = R, H_W_ratio.marine_input = R))
    expect_equal(a$particle_count_basis, "from mass (alignment geometry)")
    # particles/mL x um3 per particle = um3/mL; mass/density: mg/L / (g/cm3) -> um3/mL
    vol_dose <- a$dose.particles.mL.master * a$particle.volume.um3
    expect_equal(vol_dose, d$dose.mg.L.master / d$density.g.cm3 * 1e6, tolerance = 1e-10)
  }
  # the database's conversion is still available for comparison
  b <- suppressWarnings(align_data(d, verbose = FALSE, particle_count_source = "database"))
  expect_equal(b$dose.particles.mL.master, d$dose.particles.mL.master)
  expect_equal(b$particle_count_basis, "from mass (database conversion)")
})

test_that("a mass-dosed polydisperse test: count from the mean particle mass over its range", {
  d <- mass_record(poly = TRUE)
  a <- suppressWarnings(align_data(d, verbose = FALSE))
  # mean volume over the whole reported range (power law a.v), full-length largest particle
  Lmax <- a$size.length.max.um.used.for.conversions
  vmax <- (4 / 3) * pi * (Lmax / 2) * (a$R.ave * Lmax / 2) * (a$H_W_ratio * a$R.ave * Lmax / 2)
  mu <- mux_polyfnx(a$a.v, vmax, a$particle.volume.um3.min)
  expect_equal(a$dose.particles.mL.master * mu, d$dose.mg.L.master / d$density.g.cm3 * 1e6, tolerance = 1e-10)
  expect_equal(a$mean.mass.per.particle.mg, massfnx(mu, d$density.g.cm3) * 1e-3)
})

test_that("reported particle counts are never changed", {
  tox <- PSSDplusplus::tomex2
  rep_rows <- which(as.character(tox$dose.particles.mL.master.converted.reported) == "reported")[1:50]
  a <- suppressWarnings(align_data(tox[rep_rows, ], verbose = FALSE))
  src <- tox[rep_rows, ][match(a$rowid, tox$rowid[rep_rows]), ]
  expect_equal(a$dose.particles.mL.master, src$dose.particles.mL.master)
  expect_true(all(a$particle_count_basis == "as reported"))
})

test_that("unit conversions use the ERM values of 1 and 5000 um particles", {
  thr <- data.frame(Tier = "Tier3", `Tissue Translocation (Default)` = 10,
                    `Food Dilution (Default)` = 2, check.names = FALSE)
  p <- param_default_values
  out <- convert_units_fxn(thr, "Marine", p)
  v <- volumefnx(length = c(1, 5000), width = c(1, 5000), height = c(1, 5000))
  sa <- SAfnx(length = c(1, 5000), width = c(1, 5000), height = c(1, 5000))
  m <- massfnx(v, p$p.ave.marine)
  expect_equal(out$food_dilution_um3.L, 2 * mux_polyfnx(p$a.v.marine, v[2], v[1]))
  expect_equal(out$tissue_translocation_um2.L, 10 * mux_polyfnx(p$a.sa.marine, sa[2], sa[1]))
  expect_equal(out$food_dilution_ug.L, 2 * mux_polyfnx(p$a.m.marine, m[2], m[1]))
  # a mean particle is far heavier than a 1 um sphere and far lighter than a 5 mm one
  mu_m <- out$food_dilution_ug.L / 2
  expect_gt(mu_m, m[1]); expect_lt(mu_m, m[2])
})

test_that("a monodisperse fiber with the default width is never wider than long", {
  tox <- PSSDplusplus::tomex2
  i <- which(tox$polydispersity == "monodisperse" & tox$shape_f == "Fiber" &
    is.na(tox$size.width.um.used.for.conversions) & tox$env_f == "Freshwater")[1]
  d <- tox[i, ]
  d$size.length.um.used.for.conversions <- 2
  a <- suppressWarnings(align_data(d, verbose = FALSE))
  expect_equal(a$size.width.um.used.for.conversions, 2)
  d$size.length.um.used.for.conversions <- 200
  a <- suppressWarnings(align_data(d, verbose = FALSE))
  expect_equal(a$size.width.um.used.for.conversions, 15)
})
