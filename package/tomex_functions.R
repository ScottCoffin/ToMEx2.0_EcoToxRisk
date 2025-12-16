# Minimal ToMEx Monte Carlo + PSSD++ workflow helpers
# Source: distilled from scripts/monte carlo/EcoTox_MonteCarlo.Rmd and scripts/utils/*
# Depends on: tidyverse, truncnorm, ssdtools

## Core helpers ------------------------------------------------------------

# Correction factor from Koelmans et al.
CFfnx <- function(a, x2D, x1D, x2M, x1M) {
  (x2D^(1 - a) - x1D^(1 - a)) / (x2M^(1 - a) - x1M^(1 - a))
}

# Mu_x for polydisperse particles (Kooi et al. 2021)
mux_polyfnx <- function(a.x, x_UL, x_LL) {
  if (length(a.x) != length(x_UL) || length(a.x) != length(x_LL)) {
    stop("a.x, x_UL, and x_LL must have the same length.")
  }
  mux.poly <- numeric(length(a.x))
  for (i in seq_along(a.x)) {
    if (is.na(a.x[i]) || is.na(x_UL[i]) || is.na(x_LL[i])) {
      mux.poly[i] <- NA
    } else if (a.x[i] == 1) {
      mux.poly[i] <- (x_UL[i] - x_LL[i]) / log(x_UL[i] / x_LL[i])
    } else if (a.x[i] == 2) {
      epsilon <- 1e-10
      mux.poly[i] <- log(x_UL[i] / x_LL[i]) /
        ((x_LL[i] + epsilon)^-1 - (x_UL[i] + epsilon)^-1)
    } else {
      mux.poly[i] <- ((1 - a.x[i]) / (2 - a.x[i])) *
        ((x_UL[i]^(2 - a.x[i]) - x_LL[i]^(2 - a.x[i])) /
           (x_UL[i]^(1 - a.x[i]) - x_LL[i]^(1 - a.x[i])))
    }
  }
  mux.poly
}

## 1) Environmental parameter sampler -------------------------------------

generate_data <- function(n_sim = 10,
                          R.ave.water.marine = 0.77, R.ave.water.marine.sd = 0.29,
                          R.ave.water.freshwater = 0.67, R.ave.water.freshwater.sd = 0.28,
                          R.ave.sediment.marine = 0.75, R.ave.sediment.marine.sd = 0.30,
                          R.ave.sediment.freshwater = 0.70, R.ave.sediment.freshwater.sd = 0.33,
                          p.ave.marine = 1.10, p.ave.marine.sd = 0.14,
                          p.ave.freshwater = 1.04, p.ave.freshwater.sd = 0.12,
                          alpha.marine = 2.07, alpha.marine.sd = 0.03,
                          a.sa.marine = 1.50, a.sa.marine.sd = 0.009,
                          a.v.marine = 1.48, a.v.marine.sd = 0.063,
                          a.m.marine = 1.32, a.m.marine.sd = 0.009,
                          a.ssa.marine = 1.98, a.ssa.marine.sd = 0.297,
                          alpha.freshwater = 2.64, alpha.freshwater.sd = 0.01,
                          a.sa.freshwater = 2.00, a.sa.freshwater.sd = 0.065,
                          a.v.freshwater = 1.68, a.v.freshwater.sd = 0.081,
                          a.m.freshwater = 1.65, a.m.freshwater.sd = 0.071,
                          a.ssa.freshwater = 2.71, a.ssa.freshwater.sd = 0.009,
                          beta_log10_body_length = 0.9341, se_beta_log10_body_length = 0.1376,
                          body_length_intercept = 1.1200, se_body_length_intercept = 0.3222,
                          beta_0 = 1.308344, se_beta_0 = 0.3963612,
                          beta_1 = -0.01468148, se_beta_1 = 0.006657993) {
  R.ave.water.marine.samples <- truncnorm::rtruncnorm(n_sim, mean = R.ave.water.marine, a = 1e-4, b = 1, sd = R.ave.water.marine.sd)
  R.ave.water.freshwater.samples <- truncnorm::rtruncnorm(n_sim, mean = R.ave.water.freshwater, a = 1e-4, b = 1, sd = R.ave.water.freshwater.sd)
  R.ave.sediment.marine.samples <- truncnorm::rtruncnorm(n_sim, mean = R.ave.sediment.marine, a = 1e-4, b = 1, sd = R.ave.sediment.marine.sd)
  R.ave.sediment.freshwater.samples <- truncnorm::rtruncnorm(n_sim, mean = R.ave.sediment.freshwater, a = 1e-4, b = 1, sd = R.ave.sediment.freshwater.sd)
  p.ave.marine.samples <- rnorm(n_sim, mean = p.ave.marine, sd = p.ave.marine.sd)
  p.ave.freshwater.samples <- rnorm(n_sim, mean = p.ave.freshwater, sd = p.ave.freshwater.sd)
  alpha.marine.samples <- rnorm(n_sim, mean = alpha.marine, sd = alpha.marine.sd)
  a.sa.marine.samples <- rnorm(n_sim, mean = a.sa.marine, sd = a.sa.marine.sd)
  a.v.marine.samples <- rnorm(n_sim, mean = a.v.marine, sd = a.v.marine.sd)
  a.m.marine.samples <- rnorm(n_sim, mean = a.m.marine, sd = a.m.marine.sd)
  a.ssa.marine.samples <- rnorm(n_sim, mean = a.ssa.marine, sd = a.ssa.marine.sd)
  alpha.freshwater.samples <- rnorm(n_sim, mean = alpha.freshwater, sd = alpha.freshwater.sd)
  a.sa.freshwater.samples <- rnorm(n_sim, mean = a.sa.freshwater, sd = a.sa.freshwater.sd)
  a.v.freshwater.samples <- rnorm(n_sim, mean = a.v.freshwater, sd = a.v.freshwater.sd)
  a.m.freshwater.samples <- rnorm(n_sim, mean = a.m.freshwater, sd = a.m.freshwater.sd)
  a.ssa.freshwater.samples <- rnorm(n_sim, mean = a.ssa.freshwater, sd = a.ssa.freshwater.sd)
  sim_beta_log10_body_length_samples <- rnorm(n_sim, mean = beta_log10_body_length, sd = se_beta_log10_body_length)
  sim_body_length_intercept_samples <- rnorm(n_sim, mean = body_length_intercept, sd = se_body_length_intercept)
  sim_beta_0 <- rnorm(ceiling(n_sim * 1.2), mean = beta_0, sd = se_beta_0)
  sim_beta_1 <- rnorm(ceiling(n_sim * 1.2), mean = beta_1, sd = se_beta_1)
  sim_X50 <- -sim_beta_0 / sim_beta_1
  upper.tissue.trans.size.um.samples <- sim_X50[sim_X50 > 0][seq_len(n_sim)]
  data.frame(
    alpha.marine = alpha.marine.samples,
    a.sa.marine = a.sa.marine.samples,
    a.v.marine = a.v.marine.samples,
    a.m.marine = a.m.marine.samples,
    a.ssa.marine = a.ssa.marine.samples,
    R.ave.water.marine = R.ave.water.marine.samples,
    H_W_ratio.marine = R.ave.water.marine.samples,
    alpha.freshwater = alpha.freshwater.samples,
    a.sa.freshwater = a.sa.freshwater.samples,
    a.v.freshwater = a.v.freshwater.samples,
    a.m.freshwater = a.m.freshwater.samples,
    a.ssa.freshwater = a.ssa.freshwater.samples,
    R.ave.water.freshwater = R.ave.water.freshwater.samples,
    H_W_ratio.freshwater = R.ave.water.freshwater.samples,
    sim_beta_log10_body_length = sim_beta_log10_body_length_samples,
    sim_body_length_intercept = sim_body_length_intercept_samples,
    upper.tissue.trans.size.um = upper.tissue.trans.size.um.samples
  )
}

## 2) Threshold unit conversions -------------------------------------------

convert_units_fxn <- function(thresholds_df, environment, params, x1D_set = 1, x2D_set = 5000) {
  if (environment == "Marine") {
    a.sa.input <- params$a.sa.marine
    a.m.input <- params$a.m.marine
    a.v.input <- params$a.v.marine
  } else {
    a.sa.input <- params$a.sa.freshwater
    a.m.input <- params$a.m.freshwater
    a.v.input <- params$a.v.freshwater
  }
  thresholds_df %>%
    mutate(
      tissue_translocation_um2.L = `Tissue Translocation (Default)` *
        mux_polyfnx(a.x = a.sa.input, x_UL = x2D_set, x_LL = x1D_set),
      tissue_translocation_ug.L = `Tissue Translocation (Default)` *
        mux_polyfnx(a.x = a.m.input, x_UL = x2D_set, x_LL = x1D_set),
      food_dilution_um3.L = `Food Dilution (Default)` *
        mux_polyfnx(a.x = a.v.input, x_UL = x2D_set, x_LL = x1D_set),
      food_dilution_ug.L = `Food Dilution (Default)` *
        mux_polyfnx(a.x = a.m.input, x_UL = x2D_set, x_LL = x1D_set)
    )
}

## 3) Alignment -------------------------------------------------------------

align_data <- function(df,
                       x1M_set_input = 1,
                       x1D_set_input = 1,
                       x2D_set_input = 5000,
                       upper.tissue.trans.size.um_input = 88,
                       R.ave.marine_input = 0.77,
                       H_W_ratio.marine_input = 0.77,
                       R.ave.freshwater_input = 0.67,
                       H_W_ratio.freshwater_input = 0.67,
                       R.ave.sediment.marine_input = 0.75,
                       R.ave.sediment.freshwater_input = 0.70,
                       p.ave.marine_input = 1.10,
                       alpha.marine_input = 2.07,
                       a.sa.marine_input = 1.50,
                       a.v.marine_input = 1.48,
                       a.m.marine_input = 1.32,
                       a.ssa.marine_input = 1.98,
                       p.ave.freshwater_input = 1.04,
                       alpha.freshwater_input = 2.64,
                       a.sa.freshwater_input = 2.00,
                       a.v.freshwater_input = 1.68,
                       a.m.freshwater_input = 1.65,
                       a.ssa.freshwater_input = 2.71,
                       beta_log10_body_length_input = 0.9341,
                       body_length_intercept_input = 1.1200) {
  if (!"dose.mg.kg.sed.measured" %in% names(df)) {
    df <- df %>% mutate(dose.mg.kg.sed.measured = measured.dose.mg.kg.sediment)
  }
  if (!"dose.mg.kg.sed.nominal" %in% names(df)) {
    df <- df %>% mutate(dose.mg.kg.sed.nominal = nominal.dose.mg.kg.sediment)
  }
  if (!"dose.particles.kg.sed.nominal" %in% names(df)) {
    df <- df %>% mutate(dose.particles.kg.sed.nominal = nominal.dose.particles.kg.sediment)
  }
  if (!"environment" %in% names(df)) {
    df <- df %>% mutate(environment = env_f)
  }
  df <- df %>%
    ungroup() %>%
    mutate(
      x1M_set = x1M_set_input,
      x1D_set = x1D_set_input,
      x2D_set = x2D_set_input,
      upper.tissue.trans.size.um = upper.tissue.trans.size.um_input,
      beta_log10_body_length = beta_log10_body_length_input,
      body_length_intercept = body_length_intercept_input,
      H_W_ratio.marine = H_W_ratio.marine_input,
      H_W_ratio.freshwater = H_W_ratio.freshwater_input,
      R.ave.marine = R.ave.marine_input,
      R.ave.freshwater = R.ave.freshwater_input,
      R.ave.sediment.marine = R.ave.sediment.marine_input,
      R.ave.sediment.freshwater = R.ave.sediment.freshwater_input,
      p.ave.marine = p.ave.marine_input,
      alpha.marine = alpha.marine_input,
      a.sa.marine = a.sa.marine_input,
      a.v.marine = a.v.marine_input,
      a.m.marine = a.m.marine_input,
      a.ssa.marine = a.ssa.marine_input,
      p.ave.freshwater = p.ave.freshwater_input,
      alpha.freshwater = alpha.freshwater_input,
      a.sa.freshwater = a.sa.freshwater_input,
      a.v.freshwater = a.v.freshwater_input,
      a.m.freshwater = a.m.freshwater_input,
      a.ssa.freshwater = a.ssa.freshwater_input
    ) %>%
    mutate(
      alpha = ifelse(environment == "Marine", alpha.marine, alpha.freshwater),
      a.sa = ifelse(environment == "Marine", a.sa.marine, a.sa.freshwater),
      a.v = ifelse(environment == "Marine", a.v.marine, a.v.freshwater),
      a.m = ifelse(environment == "Marine", a.m.marine, a.m.freshwater),
      a.ssa = ifelse(environment == "Marine", a.ssa.marine, a.ssa.freshwater),
      R.ave = ifelse(environment == "Marine", R.ave.marine, R.ave.freshwater),
      H_W_ratio = ifelse(environment == "Marine", H_W_ratio.marine, H_W_ratio.freshwater),
      p.ave = ifelse(environment == "Marine", p.ave.marine, p.ave.freshwater)
    ) %>%
    mutate(
      x2M_ingest = ifelse(is.na(max.size.ingest.um), x2D_set, pmin(max.size.ingest.um, x2D_set)),
      x2M_trans = ifelse(is.na(max.size.ingest.um), upper.tissue.trans.size.um, pmin(max.size.ingest.um, upper.tissue.trans.size.um))
    ) %>%
    mutate(
      ingestible = ifelse(polydispersity == "monodisperse" & size.length.um.used.for.conversions <= x2M_ingest, "ingestible", "not ingestible"),
      translocatable = ifelse(polydispersity == "monodisperse" & size.length.um.used.for.conversions <= x2M_trans, "translocatable", "not translocatable")
    ) %>%
    mutate(
      ingestible_poly = case_when(
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions <= x2M_ingest & size.length.min.um.used.for.conversions <= x2M_ingest ~ "ingestible (all)",
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions > x2M_ingest & size.length.min.um.used.for.conversions <= x2M_ingest ~ "ingestible (some)",
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions > x2M_ingest & size.length.min.um.used.for.conversions > x2M_ingest ~ "not ingestible"
      ),
      translocatable_poly = case_when(
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions <= x2M_trans & size.length.min.um.used.for.conversions <= x2M_trans ~ "translocatable (all)",
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions > x2M_trans & size.length.min.um.used.for.conversions <= x2M_trans ~ "translocatable (some)",
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions > x2M_trans & size.length.min.um.used.for.conversions > x2M_trans ~ "not translocatable"
      )
    ) %>%
    mutate(
      CF_bioavailable_trans = ifelse(translocatable_poly == "translocatable (some)",
                                     CFfnx(a = alpha,
                                           x1D = size.length.min.um.used.for.conversions,
                                           x2D = x2M_trans,
                                           x1M = size.length.min.um.used.for.conversions,
                                           x2M = size.length.max.um.used.for.conversions), 1),
      dose.particles.mL.trans = ifelse(translocatable_poly == "translocatable (some)",
                                       CF_bioavailable_trans * dose.particles.mL.master,
                                       dose.particles.mL.master),
      CF_bioavailable_ingest = ifelse(ingestible_poly == "ingestible (some)",
                                      CFfnx(a = alpha,
                                            x1D = size.length.min.um.used.for.conversions,
                                            x2D = x2M_ingest,
                                            x1M = size.length.min.um.used.for.conversions,
                                            x2M = size.length.max.um.used.for.conversions), 1),
      dose.particles.mL.ingest = ifelse(ingestible_poly == "ingestible (some)",
                                        CF_bioavailable_ingest * dose.particles.mL.master,
                                        dose.particles.mL.master),
      size.length.max.um.trans = ifelse(translocatable_poly == "translocatable (some)", x2M_trans, size.length.max.um.used.for.conversions),
      size.length.max.um.ingest = ifelse(ingestible_poly == "ingestible (some)", x2M_ingest, size.length.max.um.used.for.conversions)
    ) %>%
    mutate(
      size.width.um.used.for.conversions = case_when(
        is.na(size.width.um.used.for.conversions) & shape_f == "Fiber" ~ 15,
        is.na(size.width.um.used.for.conversions) & shape_f == "Sphere" ~ size.length.um.used.for.conversions,
        is.na(size.width.um.used.for.conversions) & shape_f == "Fragment" ~ size.length.um.used.for.conversions * R.ave,
        TRUE ~ size.width.um.used.for.conversions
      ),
      size.height.um.used.for.conversions = case_when(
        is.na(size.height.um.used.for.conversions) & shape_f == "Fiber" ~ 15,
        is.na(size.height.um.used.for.conversions) & shape_f == "Sphere" ~ size.length.um.used.for.conversions,
        is.na(size.height.um.used.for.conversions) & shape_f == "Fragment" ~ size.width.um.used.for.conversions * H_W_ratio,
        TRUE ~ size.height.um.used.for.conversions
      )
    ) %>%
    mutate(
      mu.v.poly.ingest = mux_polyfnx(a.x = a.v, x_UL = size.length.max.um.ingest, x_LL = size.length.min.um.used.for.conversions),
      mu.m.poly.ingest = mux_polyfnx(a.x = a.m, x_UL = size.length.max.um.ingest, x_LL = size.length.min.um.used.for.conversions),
      mu.sa.poly.ingest = mux_polyfnx(a.x = a.sa, x_UL = size.length.max.um.ingest, x_LL = size.length.min.um.used.for.conversions),
      mu.ssa.poly.ingest = mux_polyfnx(a.x = a.ssa, x_UL = size.length.max.um.ingest, x_LL = size.length.min.um.used.for.conversions),
      mu.v.poly.trans = mux_polyfnx(a.x = a.v, x_UL = size.length.max.um.trans, x_LL = size.length.min.um.used.for.conversions),
      mu.m.poly.trans = mux_polyfnx(a.x = a.m, x_UL = size.length.max.um.trans, x_LL = size.length.min.um.used.for.conversions),
      mu.sa.poly.trans = mux_polyfnx(a.x = a.sa, x_UL = size.length.max.um.trans, x_LL = size.length.min.um.used.for.conversions),
      mu.ssa.poly.trans = mux_polyfnx(a.x = a.ssa, x_UL = size.length.max.um.trans, x_LL = size.length.min.um.used.for.conversions)
    ) %>%
    mutate(
      CF_aligned_trans = mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set) /
        mux_polyfnx(a.x = a.sa, x_UL = size.length.max.um.trans, x_LL = size.length.min.um.used.for.conversions),
      CF_aligned_ingest = mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set) /
        mux_polyfnx(a.x = a.v, x_UL = size.length.max.um.ingest, x_LL = size.length.min.um.used.for.conversions)
    ) %>%
    mutate(
      EC_poly_p.particles.mL_trans = case_when(
        polydispersity == "polydisperse" ~ dose.particles.mL.trans / mu.sa.poly.trans,
        polydispersity == "monodisperse" ~ dose.particles.mL.trans
      ),
      EC_poly_p.particles.mL_ingest = case_when(
        polydispersity == "polydisperse" ~ dose.particles.mL.ingest / mu.v.poly.ingest,
        polydispersity == "monodisperse" ~ dose.particles.mL.ingest
      ),
      EC_env_sa.particles.mL_trans = EC_poly_p.particles.mL_trans * CF_aligned_trans,
      EC_env_v.particles.mL_ingest = EC_poly_p.particles.mL_ingest * CF_aligned_ingest
    ) %>%
    mutate(
      particles.mL.ox.stress = EC_env_sa.particles.mL_trans,
      particles.mL.food.dilution = ifelse(Group == "Algae", NA, EC_env_v.particles.mL_ingest)
    ) %>%
    mutate(unique_id = row_number()) %>%
    filter(size.length.um.used.for.conversions >= x1D_set_input)
  df
}

## 3) SSD helpers ----------------------------------------------------------

SSD_function_t1 <- function(filtered.data, hcxlcl, nboot = 10) {
  collapsed <- filtered.data %>%
    group_by(Species, Group) %>%
    summarize(Conc = quantile(dose_new, 0.25), .groups = "drop")
  metadata <- filtered.data %>%
    summarize(n_species = n_distinct(Species),
              n_groups = n_distinct(Group),
              n_datapoints = n(),
              .groups = "drop")
  if (nrow(collapsed) < 6) {
    return(list(hc5lcl = NA, n_species = metadata$n_species, n_groups = metadata$n_groups, n_datapoints = metadata$n_datapoints))
  }
  dists <- ssdtools::ssd_fit_dists(collapsed, left = "Conc",
                                   dists = c("llogis", "lnorm", "gamma", "lgumbel"),
                                   computable = FALSE, silent = FALSE)
  preds <- predict(dists, average = TRUE, nboot = nboot, ci = TRUE)
  list(hc5lcl = c(preds$lcl[hcxlcl]),
       n_species = metadata$n_species,
       n_groups = metadata$n_groups,
       n_datapoints = metadata$n_datapoints)
}

SSD_function_t2 <- function(filtered.data, hcx, nboot = 10) {
  collapsed <- filtered.data %>%
    group_by(Species, Group) %>%
    summarize(Conc = quantile(dose_new, 0.25), .groups = "drop")
  metadata <- filtered.data %>%
    summarize(n_species = n_distinct(Species),
              n_groups = n_distinct(Group),
              n_datapoints = n(),
              .groups = "drop")
  if (nrow(collapsed) < 6) {
    return(list(hcx_est = NA, hcx05cl = NA, hcx95cl = NA,
                n_species = metadata$n_species, n_groups = metadata$n_groups, n_datapoints = metadata$n_datapoints))
  }
  dists <- ssdtools::ssd_fit_dists(collapsed, left = "Conc",
                                   dists = c("llogis", "lnorm", "gamma", "lgumbel"),
                                   computable = FALSE, silent = FALSE)
  preds <- predict(dists, average = TRUE, nboot = nboot, ci = TRUE)
  list(hcx_est = c(preds$est[hcx]),
       hcx05cl = c(preds$lcl[hcx]),
       hcx95cl = c(preds$ucl[hcx]),
       n_species = metadata$n_species,
       n_groups = metadata$n_groups,
       n_datapoints = metadata$n_datapoints)
}

SSD_function_t3_4 <- function(filtered.data, hcx, nboot = 10) {
  collapsed <- filtered.data %>%
    filter(risk.13 != 1, bio_f %in% c("Organism", "Population")) %>%
    group_by(Species, Group) %>%
    summarize(Conc = quantile(dose_new, 0.50), .groups = "drop")
  metadata <- filtered.data %>%
    filter(risk.13 != 1, bio_f %in% c("Organism", "Population")) %>%
    summarize(n_species = n_distinct(Species),
              n_groups = n_distinct(Group),
              n_datapoints = n(),
              .groups = "drop")
  if (nrow(collapsed) < 6) {
    return(list(hcx_est = NA, hcx05cl = NA, hcx95cl = NA,
                n_species = metadata$n_species, n_groups = metadata$n_groups, n_datapoints = metadata$n_datapoints))
  }
  dists <- ssdtools::ssd_fit_dists(collapsed, left = "Conc",
                                   dists = c("llogis", "lnorm", "gamma", "lgumbel"),
                                   computable = FALSE, silent = FALSE)
  preds <- predict(dists, average = TRUE, nboot = nboot, ci = TRUE)
  list(hcx_est = c(preds$est[hcx]),
       hcx05cl = c(preds$lcl[hcx]),
       hcx95cl = c(preds$ucl[hcx]),
       n_species = metadata$n_species,
       n_groups = metadata$n_groups,
       n_datapoints = metadata$n_datapoints)
}

process_environment_data <- function(data, env_filter, upper.tissue.trans.size.um, x1D_set, x2D_set, nboot = 10) {
  filtered_data_small_default_t1_2 <- data %>%
    ungroup() %>%
    mutate(dose_new = particles.mL.ox.stress / (af.time * af.noec)) %>%
    drop_na(dose_new) %>%
    filter(dose_new > 0) %>%
    mutate(dose_new = dose_new * 1000) %>%
    filter(between(size.length.um.used.for.conversions, x1D_set, upper.tissue.trans.size.um),
           shape_f != "Not Reported",
           poly_f != "Not Reported",
           environment %in% env_filter,
           Group != "Bacterium",
           Group != "Plant",
           effect.metric != "HONEC",
           translocatable != "not translocatable")
  
  filtered_data_small_default_t3_4 <- filtered_data_small_default_t1_2 %>%
    filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
  filtered_data_large_default_t1_2 <- data %>%
    filter(Group != "Algae") %>%
    mutate(dose_new = particles.mL.food.dilution / (af.time * af.noec)) %>%
    drop_na(dose_new) %>%
    filter(dose_new > 0) %>%
    mutate(dose_new = dose_new * 1000) %>%
    filter(between(size.length.um.used.for.conversions, x1D_set, x2D_set),
           ingestible != "not ingestible",
           poly_f != "Not Reported",
           shape_f != "Not Reported",
           environment %in% env_filter,
           Group != "Bacterium",
           Group != "Plant",
           effect.metric != "HONEC")
  filtered_data_large_default_t3_4 <- filtered_data_large_default_t1_2 %>%
    filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
  small_default_t1 <- SSD_function_t1(filtered_data_small_default_t1_2, hcxlcl = 5, nboot = nboot)
  small_default_t2 <- SSD_function_t2(filtered_data_small_default_t1_2, hcx = 5, nboot = nboot)
  small_default_t3 <- SSD_function_t3_4(filtered_data_small_default_t3_4, hcx = 5, nboot = nboot)
  small_default_t4 <- SSD_function_t3_4(filtered_data_small_default_t3_4, hcx = 10, nboot = nboot)
  large_default_t1 <- SSD_function_t1(filtered_data_large_default_t1_2, hcxlcl = 5, nboot = nboot)
  large_default_t2 <- SSD_function_t2(filtered_data_large_default_t1_2, hcx = 5, nboot = nboot)
  large_default_t3 <- SSD_function_t3_4(filtered_data_large_default_t3_4, hcx = 5, nboot = nboot)
  large_default_t4 <- SSD_function_t3_4(filtered_data_large_default_t3_4, hcx = 10, nboot = nboot)
  tibble::tibble(
    Tier = c("Tier1", "Tier2", "Tier3", "Tier4"),
    `Tissue Translocation (Default)` = c(small_default_t1$hc5lcl, small_default_t2$hcx_est, small_default_t3$hcx_est, small_default_t4$hcx_est),
    `Food Dilution (Default)` = c(large_default_t1$hc5lcl, large_default_t2$hcx_est, large_default_t3$hcx_est, large_default_t4$hcx_est),
    `Tissue Translocation (5th %)` = c(NA, small_default_t2$hcx05cl, small_default_t3$hcx05cl, small_default_t4$hcx05cl),
    `Tissue Translocation (95th %)` = c(NA, small_default_t2$hcx95cl, small_default_t3$hcx95cl, small_default_t4$hcx95cl),
    `Food Dilution (5th %)` = c(NA, large_default_t2$hcx05cl, large_default_t3$hcx05cl, large_default_t4$hcx05cl),
    `Food Dilution (95th %)` = c(NA, large_default_t2$hcx95cl, large_default_t3$hcx95cl, large_default_t4$hcx95cl),
    `N Species (tissue trans)` = c(small_default_t1$n_species, small_default_t2$n_species, small_default_t3$n_species, small_default_t4$n_species),
    `N Species (food dilution)` = c(large_default_t1$n_species, large_default_t2$n_species, large_default_t3$n_species, large_default_t4$n_species),
    `N Groups (tissue trans)` = c(small_default_t1$n_groups, small_default_t2$n_groups, small_default_t3$n_groups, small_default_t4$n_groups),
    `N Groups (food dilution)` = c(large_default_t1$n_groups, large_default_t2$n_groups, large_default_t3$n_groups, large_default_t4$n_groups),
    `N Datapoints (tissue trans)` = c(small_default_t1$n_datapoints, small_default_t2$n_datapoints, small_default_t3$n_datapoints, small_default_t4$n_datapoints),
    `N Datapoints (food dilution)` = c(large_default_t1$n_datapoints, large_default_t2$n_datapoints, large_default_t3$n_datapoints, large_default_t4$n_datapoints)
  )
}

## 4) Wrapper for Sobol/MC iteration ---------------------------------------

model_wrapper_sobol_parallel <- function(data, params_row, simulation_id,
                                         x1D_set = 1, x2D_set = 5000, nboot = 10) {
  aligned <- align_data(
    data,
    upper.tissue.trans.size.um_input = params_row$upper.tissue.trans.size.um,
    body_length_intercept_input = params_row$sim_body_length_intercept,
    beta_log10_body_length_input = params_row$sim_beta_log10_body_length,
    R.ave.marine_input = params_row$R.ave.water.marine,
    H_W_ratio.marine_input = params_row$H_W_ratio.marine,
    alpha.marine_input = params_row$alpha.marine,
    a.sa.marine_input = params_row$a.sa.marine,
    a.v.marine_input = params_row$a.v.marine,
    a.m.marine_input = params_row$a.m.marine,
    R.ave.freshwater_input = params_row$R.ave.water.freshwater,
    H_W_ratio.freshwater_input = params_row$H_W_ratio.freshwater,
    alpha.freshwater_input = params_row$alpha.freshwater,
    a.sa.freshwater_input = params_row$a.sa.freshwater,
    a.v.freshwater_input = params_row$a.v.freshwater,
    a.m.freshwater_input = params_row$a.m.freshwater
  )
  aligned <- aligned %>%
    drop_na(effect.metric) %>%
    filter(tier_zero_tech_f == "Red Criteria Passed")
  marine_thresholds <- process_environment_data(
    aligned, "Marine",
    upper.tissue.trans.size.um = params_row$upper.tissue.trans.size.um,
    x1D_set = x1D_set, x2D_set = x2D_set, nboot = nboot
  ) %>%
    convert_units_fxn(environment = "Marine", params = params_row, x1D_set = x1D_set, x2D_set = x2D_set)
  freshwater_thresholds <- process_environment_data(
    aligned, "Freshwater",
    upper.tissue.trans.size.um = params_row$upper.tissue.trans.size.um,
    x1D_set = x1D_set, x2D_set = x2D_set, nboot = nboot
  ) %>%
    convert_units_fxn(environment = "Freshwater", params = params_row, x1D_set = x1D_set, x2D_set = x2D_set)
  list(
    simulation_id = simulation_id,
    base_thresholds = list(marine = marine_thresholds, freshwater = freshwater_thresholds),
    tox_vals = aligned
  )
}

## Reprex pointers ----------------------------------------------------------
# generate_data(n_sim = 5)
# aligned_det <- align_data(readRDS("data/input/aoc_z_tomex2.RDS"), x2D_set_input = 5000, upper.tissue.trans.size.um_input = 88)
# thresholds_demo <- process_environment_data(aligned_det, c("Marine","Freshwater"), upper.tissue.trans.size.um = 88, x1D_set = 1, x2D_set = 5000, nboot = 2)
