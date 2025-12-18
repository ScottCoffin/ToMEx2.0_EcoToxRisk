# SSD helper functions ----------------------------------------------------

SSD_function_t1 <- function(filtered.data, hcxlcl, nboot = 10) {
  collapsed <- filtered.data %>%
    dplyr::group_by(Species, Group) %>%
    dplyr::summarize(Conc = stats::quantile(dose_new, 0.25), .groups = "drop")

  metadata <- filtered.data %>%
    dplyr::summarize(
      n_species = dplyr::n_distinct(Species),
      n_groups = dplyr::n_distinct(Group),
      n_datapoints = dplyr::n(),
      .groups = "drop"
    )

  if (nrow(collapsed) < 6) {
    return(list(
      hc5lcl = NA,
      n_species = metadata$n_species,
      n_groups = metadata$n_groups,
      n_datapoints = metadata$n_datapoints
    ))
  }

  dists <- ssdtools::ssd_fit_dists(
    collapsed,
    left = "Conc",
    dists = c("llogis", "lnorm", "gamma", "lgumbel"),
    computable = FALSE,
    silent = FALSE
  )
  preds <- stats::predict(dists, average = TRUE, nboot = nboot, ci = TRUE)

  list(
    hc5lcl = c(preds$lcl[hcxlcl]),
    n_species = metadata$n_species,
    n_groups = metadata$n_groups,
    n_datapoints = metadata$n_datapoints
  )
}

SSD_function_t2 <- function(filtered.data, hcx, nboot = 10) {
  collapsed <- filtered.data %>%
    dplyr::group_by(Species, Group) %>%
    dplyr::summarize(Conc = stats::quantile(dose_new, 0.25), .groups = "drop")

  metadata <- filtered.data %>%
    dplyr::summarize(
      n_species = dplyr::n_distinct(Species),
      n_groups = dplyr::n_distinct(Group),
      n_datapoints = dplyr::n(),
      .groups = "drop"
    )

  if (nrow(collapsed) < 6) {
    return(list(
      hcx_est = NA,
      hcx05cl = NA,
      hcx95cl = NA,
      n_species = metadata$n_species,
      n_groups = metadata$n_groups,
      n_datapoints = metadata$n_datapoints
    ))
  }

  dists <- ssdtools::ssd_fit_dists(
    collapsed,
    left = "Conc",
    dists = c("llogis", "lnorm", "gamma", "lgumbel"),
    computable = FALSE,
    silent = FALSE
  )
  preds <- stats::predict(dists, average = TRUE, nboot = nboot, ci = TRUE)

  list(
    hcx_est = c(preds$est[hcx]),
    hcx05cl = c(preds$lcl[hcx]),
    hcx95cl = c(preds$ucl[hcx]),
    n_species = metadata$n_species,
    n_groups = metadata$n_groups,
    n_datapoints = metadata$n_datapoints
  )
}

SSD_function_t3_4 <- function(filtered.data, hcx, nboot = 10) {
  collapsed <- filtered.data %>%
    dplyr::filter(risk.13 != 1, bio_f %in% c("Organism", "Population")) %>%
    dplyr::group_by(Species, Group) %>%
    dplyr::summarize(Conc = stats::quantile(dose_new, 0.50), .groups = "drop")

  metadata <- filtered.data %>%
    dplyr::filter(risk.13 != 1, bio_f %in% c("Organism", "Population")) %>%
    dplyr::summarize(
      n_species = dplyr::n_distinct(Species),
      n_groups = dplyr::n_distinct(Group),
      n_datapoints = dplyr::n(),
      .groups = "drop"
    )

  if (nrow(collapsed) < 6) {
    return(list(
      hcx_est = NA,
      hcx05cl = NA,
      hcx95cl = NA,
      n_species = metadata$n_species,
      n_groups = metadata$n_groups,
      n_datapoints = metadata$n_datapoints
    ))
  }

  dists <- ssdtools::ssd_fit_dists(
    collapsed,
    left = "Conc",
    dists = c("llogis", "lnorm", "gamma", "lgumbel"),
    computable = FALSE,
    silent = FALSE
  )
  preds <- stats::predict(dists, average = TRUE, nboot = nboot, ci = TRUE)

  list(
    hcx_est = c(preds$est[hcx]),
    hcx05cl = c(preds$lcl[hcx]),
    hcx95cl = c(preds$ucl[hcx]),
    n_species = metadata$n_species,
    n_groups = metadata$n_groups,
    n_datapoints = metadata$n_datapoints
  )
}

# Threshold aggregation used for manuscript-style outputs -----------------

process_environment_data <- function(data,
                                     env_filter = "Marine",
                                     upper.tissue.trans.size.um = 88,
                                     x1D_set = 1,
                                     x2D_set = 5000,
                                     nboot = 10) {
  if (!"particles.mL.ox.stress" %in% names(data)) {
    data <- data %>%
      dplyr::mutate(particles.mL.ox.stress = EC_env_sa.particles.mL_trans)
  }
  if (!"particles.mL.food.dilution" %in% names(data)) {
    data <- data %>%
      dplyr::mutate(particles.mL.food.dilution = EC_env_v.particles.mL_ingest)
  }

  filtered_data_small_default_t1_2 <- data %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dose_new = particles.mL.ox.stress / (af.time * af.noec)) %>%
    tidyr::drop_na(dose_new) %>%
    dplyr::filter(
      dose_new > 0,
      dplyr::between(size.length.um.used.for.conversions, x1D_set, upper.tissue.trans.size.um),
      translocatable != "not translocatable",
      shape_f != "Not Reported",
      poly_f != "Not Reported",
      environment %in% env_filter,
      Group != "Bacterium",
      Group != "Plant",
      effect.metric != "HONEC"
    )
  filtered_data_small_default_t3_4 <- filtered_data_small_default_t1_2 %>%
    dplyr::filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))

  filtered_data_large_default_t1_2 <- data %>%
    dplyr::filter(Group != "Algae") %>%
    dplyr::mutate(dose_new = particles.mL.food.dilution / (af.time * af.noec)) %>%
    tidyr::drop_na(dose_new) %>%
    dplyr::filter(dose_new > 0) %>%
    dplyr::mutate(dose_new = dose_new * 1000) %>%
    dplyr::filter(
      dplyr::between(size.length.um.used.for.conversions, x1D_set, x2D_set),
      ingestible != "not ingestible",
      poly_f != "Not Reported",
      shape_f != "Not Reported",
      environment %in% env_filter,
      Group != "Bacterium",
      Group != "Plant",
      effect.metric != "HONEC",
      dose_new > 0
    )
  filtered_data_large_default_t3_4 <- filtered_data_large_default_t1_2 %>%
    dplyr::filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))

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
    `Tissue Translocation (Default)` = c(
      small_default_t1$hc5lcl,
      small_default_t2$hcx_est,
      small_default_t3$hcx_est,
      small_default_t4$hcx_est
    ),
    `Food Dilution (Default)` = c(
      large_default_t1$hc5lcl,
      large_default_t2$hcx_est,
      large_default_t3$hcx_est,
      large_default_t4$hcx_est
    ),
    `Tissue Translocation (5th %)` = c(NA, small_default_t2$hcx05cl, small_default_t3$hcx05cl, small_default_t4$hcx05cl),
    `Tissue Translocation (95th %)` = c(NA, small_default_t2$hcx95cl, small_default_t3$hcx95cl, small_default_t4$hcx95cl),
    `Food Dilution (5th %)` = c(NA, large_default_t2$hcx05cl, large_default_t3$hcx05cl, large_default_t4$hcx05cl),
    `Food Dilution (95th %)` = c(NA, large_default_t2$hcx95cl, large_default_t3$hcx95cl, large_default_t4$hcx95cl),
    `N Species (tissue trans)` = c(
      small_default_t1$n_species,
      small_default_t2$n_species,
      small_default_t3$n_species,
      small_default_t4$n_species
    ),
    `N Species (food dilution)` = c(
      large_default_t1$n_species,
      large_default_t2$n_species,
      large_default_t3$n_species,
      large_default_t4$n_species
    ),
    `N Groups (tissue trans)` = c(
      small_default_t1$n_groups,
      small_default_t2$n_groups,
      small_default_t3$n_groups,
      small_default_t4$n_groups
    ),
    `N Groups (food dilution)` = c(
      large_default_t1$n_groups,
      large_default_t2$n_groups,
      large_default_t3$n_groups,
      large_default_t4$n_groups
    ),
    `N Datapoints (tissue trans)` = c(
      small_default_t1$n_datapoints,
      small_default_t2$n_datapoints,
      small_default_t3$n_datapoints,
      small_default_t4$n_datapoints
    ),
    `N Datapoints (food dilution)` = c(
      large_default_t1$n_datapoints,
      large_default_t2$n_datapoints,
      large_default_t3$n_datapoints,
      large_default_t4$n_datapoints
    )
  )
}
