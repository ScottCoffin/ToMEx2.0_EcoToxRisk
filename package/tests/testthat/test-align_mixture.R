# align_mixture(): component-resolved alignment of co-dosed mixtures (issue #8, part C)

one_record <- function() {
  tox <- PSSDplusplus::tomex2
  i <- which(tox$polydispersity == "polydisperse" & tox$env_f == "Marine" &
    tox$shape_f == "Fragment" & !is.na(tox$max.size.ingest.um) &
    !is.na(tox$dose.particles.mL.master) & !is.na(tox$dose.mg.L.master) &
    !is.na(tox$density.g.cm3) &
    !grepl("sediment", tox$exposure.route, ignore.case = TRUE))[1]
  tox[i, ]
}

erm <- c("particles.mL.food.dilution", "particles.mL.ox.stress")

test_that("a mixture of identical halves aligns like the single material", {
  rec <- one_record()
  single <- suppressWarnings(align_data(rec, verbose = FALSE))
  comp <- rec[c(1, 1), ]
  comp$mixture_id <- "m1"
  comp$component.fraction <- c(0.5, 0.5)
  mix <- align_mixture(comp)
  expect_equal(nrow(mix), 1)
  expect_equal(mix$n_components, 2L)
  expect_equal(as.character(mix$shape_f), as.character(single$shape_f))
  for (col in erm) expect_equal(mix[[col]], single[[col]], label = col)
  expect_equal(mix$dose.particles.mL.master, rec$dose.particles.mL.master)
})

test_that("components of different shapes add up (number basis)", {
  rec <- one_record()
  comp <- rec[c(1, 1), ]
  comp$shape_f <- c("Fragment", "Sphere")
  comp$mixture_id <- "m1"
  comp$component.fraction <- c(0.3, 0.7)
  mix <- align_mixture(comp)
  # the dominant component's labels, so the record stays selectable by shape
  expect_equal(as.character(mix$shape_f), "Sphere")
  expect_equal(mix$mixture_shapes, "Fragment; Sphere")
  expect_equal(mix$n_components, 2L)
  parts <- comp
  parts$dose.particles.mL.master <- parts$dose.particles.mL.master * c(0.3, 0.7)
  a <- suppressWarnings(align_data(parts, verbose = FALSE))
  expect_true(all(a$particle_count_basis == "as reported"))
  for (col in erm) expect_equal(mix[[col]], sum(a[[col]]), label = col)
  expect_equal(mix$dose.particles.mL.master, rec$dose.particles.mL.master)
})

test_that("counts derived from the dosed mass: number and mass are both conserved", {
  rec <- one_record()
  rec$dose.particles.mL.master.converted.reported <- "converted"
  comp <- rec[c(1, 1), ]
  comp$shape_f <- c("Fragment", "Sphere")
  comp$density.g.cm3 <- c(1.05, 1.38)
  comp$mixture_id <- "m1"
  for (b in c("number", "mass")) {
    comp$component.fraction <- c(0.3, 0.7)
    comp$component.fraction.basis <- b
    mix <- align_mixture(comp)
    parts <- attr(mix, "components")
    expect_true(all(parts$particle_count_basis == "from mass (alignment geometry)"))
    n_i <- parts$.component_share * parts$dose.particles.mL.master # particles/mL per component
    mass_i <- n_i * parts$mean.mass.per.particle.mg * 1000 # mg/L per component
    # the components' masses add up to the dosed mass
    expect_equal(sum(mass_i), rec$dose.mg.L.master)
    if (b == "number") expect_equal(n_i / sum(n_i), c(0.3, 0.7))
    if (b == "mass") expect_equal(mass_i / sum(mass_i), c(0.3, 0.7))
    for (col in erm) {
      expect_equal(mix[[col]], sum(parts$.component_share * parts[[col]]), label = paste(b, col))
    }
    expect_equal(mix$dose.particles.mL.master, sum(n_i))
  }
})

test_that("mass shares of reported counts are converted to number shares", {
  rec <- one_record()
  comp <- rec[c(1, 1), ]
  comp$density.g.cm3 <- c(1, 2)
  comp$mixture_id <- "m1"
  comp$component.fraction <- c(0.5, 0.5)
  mix <- align_mixture(comp, fraction_basis = "mass")
  parts <- attr(mix, "components")
  # same sizes, twice the density: half as many particles per mg
  expect_equal(parts$.component_share, c(2, 1) / 3)
  expect_equal(mix$dose.particles.mL.master, rec$dose.particles.mL.master)
})

test_that("a component that is not ingestible contributes nothing to food dilution", {
  rec <- one_record()
  comp <- rec[c(1, 1), ]
  comp$polydispersity <- "monodisperse"
  comp$size.length.um.used.for.conversions <- c(10, 4000)
  comp$max.size.ingest.um <- 100
  comp$max.size.ingest.basis <- "reported (test)"
  comp$mixture_id <- "m1"
  comp$component.fraction <- c(0.5, 0.5)
  mix <- align_mixture(comp)
  parts <- comp[1, ]
  parts$dose.particles.mL.master <- parts$dose.particles.mL.master * 0.5
  a <- suppressWarnings(align_data(parts, verbose = FALSE))
  expect_equal(mix$particles.mL.food.dilution, a$particles.mL.food.dilution)
  expect_equal(mix$ingestible, "ingestible (mixture)")
  expect_equal(mix$size.length.um.used.for.conversions, 10)
})

test_that("mass basis: identical halves equal the whole; fractions are normalised", {
  rec <- one_record()
  whole <- rec
  whole$mixture_id <- "w"
  whole$component.fraction <- 1
  w <- align_mixture(whole, fraction_basis = "mass")
  comp <- rec[c(1, 1), ]
  comp$mixture_id <- "m1"
  comp$component.fraction <- c(1, 1)
  expect_warning(m <- align_mixture(comp, fraction_basis = "mass"), "normalising")
  for (col in erm) expect_equal(m[[col]], w[[col]], label = col)
  # denser component: fewer particles per mg, so a lower particle-based ERM dose
  comp2 <- comp
  comp2$density.g.cm3 <- c(1, 2)
  comp2$component.fraction <- c(0.5, 0.5)
  m2 <- align_mixture(comp2, fraction_basis = "mass")
  expect_true(all(is.finite(m2$particles.mL.food.dilution)))
})

test_that("align_mixture() checks its inputs", {
  rec <- one_record()
  expect_error(align_mixture(rec), "mixture_id")
  rec$mixture_id <- "a"
  rec$component.fraction <- NA_real_
  expect_error(align_mixture(rec), "fractions")
})
