#### Special Plotting Functions #####

##### visualize with histrogram fxn

parameter_histograms_function <- function(mat) {
  # Define colors for marine and freshwater parameters
  param_colors <- c(
    "Marine" = "#1f78b4", # Dark Blue
    "Freshwater" = "#a6cee3", # Light Blue
    "Body Length" = "#66c2a5", # Green shades for body length
    "Translocation" = "#fdae61" # Orange for translocation
  )

  # Prepare data for each group
  alpha_data <- mat %>%
    select(alpha.marine, alpha.freshwater) %>%
    rename(
      "Length Alpha (Marine)" = "alpha.marine",
      "Length Alpha (Freshwater)" = "alpha.freshwater"
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(Environment = if_else(grepl("Marine", name), "Marine", "Freshwater"))

  mass_data <- mat %>%
    select(a.m.marine, a.m.freshwater) %>%
    rename(
      "Mass Alpha (Marine)" = "a.m.marine",
      "Mass Alpha (Freshwater)" = "a.m.freshwater"
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(Environment = if_else(grepl("Marine", name), "Marine", "Freshwater"))

  volume_data <- mat %>%
    select(a.v.marine, a.v.freshwater) %>%
    rename(
      "Volume Alpha (Marine)" = "a.v.marine",
      "Volume Alpha (Freshwater)" = "a.v.freshwater"
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(Environment = if_else(grepl("Marine", name), "Marine", "Freshwater"))

  surface_area_data <- mat %>%
    select(a.sa.marine, a.sa.freshwater) %>%
    rename(
      "Surface Area Alpha (Marine)" = "a.sa.marine",
      "Surface Area Alpha (Freshwater)" = "a.sa.freshwater"
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(Environment = if_else(grepl("Marine", name), "Marine", "Freshwater"))

  ratio_data <- mat %>%
    select(R.ave.water.marine, R.ave.water.freshwater) %>%
    rename(
      "Length-Width Ratio (Marine)" = "R.ave.water.marine",
      "Length-Width Ratio (Freshwater)" = "R.ave.water.freshwater"
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(Environment = if_else(grepl("Marine", name), "Marine", "Freshwater"))

  H_W_ratio_data <- mat %>%
    select(H_W_ratio.marine, H_W_ratio.freshwater) %>%
    rename(
      "Height-Width Ratio (Marine)" = "H_W_ratio.marine",
      "Height-Width Ratio (Freshwater)" = "H_W_ratio.freshwater"
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(Environment = if_else(grepl("Marine", name), "Marine", "Freshwater"))

  body_length_data <- mat %>%
    select(sim.beta.log10.body.length, sim.body.length.intercept) %>%
    rename(
      "log10(Body Length, mm) - beta" = "sim.beta.log10.body.length",
      "log10(Body Length, mm) - intercept" = "sim.body.length.intercept"
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(
      Environment = case_when(
        grepl("beta", name) ~ "Beta (mm)",
        grepl("intercept", name) ~ "Intercept (mm)"
      )
    )

  translocation_data <- mat %>%
    select(upper.tissue.trans.size.um) %>%
    rename(
      "Maximum Translocatable Particle Length (µm)" = "upper.tissue.trans.size.um"
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(Environment = "Translocation")

  # Calculate common limits for alphas
  alpha_min <- min(c(
    alpha_data$value,
    mass_data$value,
    volume_data$value,
    surface_area_data$value
  )) *
    0.9
  alpha_max <- max(c(
    alpha_data$value,
    mass_data$value,
    volume_data$value,
    surface_area_data$value
  )) *
    1.1

  # alpha length
  alpha_length_plot <- alpha_data %>%
    ggplot(aes(x = value, fill = Environment)) +
    geom_density(
      aes(y = after_stat(density)),
      size = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    scale_fill_manual(values = param_colors) +
    scale_color_manual(values = param_colors) +
    theme_minimal(base_size = 12) +
    labs(x = "Value", y = "Relative Frequency") +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    ggtitle("Length Alpha")

  alpha_mass_plot <- mass_data %>%
    ggplot(aes(x = value, fill = Environment)) +
    geom_density(
      aes(y = after_stat(density)),
      size = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    scale_fill_manual(values = param_colors) +
    scale_color_manual(values = param_colors) +
    theme_minimal(base_size = 12) +
    labs(x = "Value", y = "Relative Frequency") +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    ggtitle("Mass Alpha")

  # alpha length
  alpha_volume_plot <- volume_data %>%
    ggplot(aes(x = value, fill = Environment)) +
    geom_density(
      aes(y = after_stat(density)),
      size = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    scale_fill_manual(values = param_colors) +
    scale_color_manual(values = param_colors) +
    theme_minimal(base_size = 12) +
    labs(x = "Value", y = "Relative Frequency") +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    ggtitle("Volume Alpha")

  # alpha length
  alpha_surface_area_plot <- surface_area_data %>%
    ggplot(aes(x = value, fill = Environment)) +
    geom_density(
      aes(y = after_stat(density)),
      size = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    scale_fill_manual(values = param_colors) +
    scale_color_manual(values = param_colors) +
    theme_minimal(base_size = 12) +
    labs(x = "Value", y = "Relative Frequency") +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 13),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
    ) +
    ggtitle("Surface Area Alpha")

  ## length to width
  length_width_plot <- ratio_data %>%
    ggplot(aes(x = value, fill = Environment)) +
    geom_density(
      aes(y = after_stat(density)),
      size = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    scale_fill_manual(values = param_colors) +
    scale_color_manual(values = param_colors) +
    theme_minimal(base_size = 12) +
    labs(x = "Value", y = "Relative Frequency") +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 13),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    ggtitle("Length to Width Ratio")

  ## height to width
  height_width_plot <- H_W_ratio_data %>%
    ggplot(aes(x = value, fill = Environment)) +
    geom_density(
      aes(y = after_stat(density)),
      size = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    scale_fill_manual(values = param_colors) +
    scale_color_manual(values = param_colors) +
    theme_minimal(base_size = 12) +
    labs(x = "Value", y = "Relative Frequency") +
    theme(
      # legend.position = "none",
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 13),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    ggtitle("Height to Width Ratio")

  ##### Body length
  body_length_plot <- body_length_data %>%
    ggplot(aes(x = value, fill = Environment)) +
    geom_density(
      aes(y = after_stat(density)),
      size = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    scale_fill_manual(values = c("#66c2a5", "#66c282"), name = "Parameter") +
    scale_color_manual(values = c("#66c2a5", "#66c282"), name = "Parameter") +
    theme_minimal(base_size = 12) +
    labs(x = "Value", y = "Relative Frequency") +
    theme(
      legend.position = c(0.2, 0.7),
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 13),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
    ) +
    ggtitle("Body Length to Mouth Size Estimation")

  tissue_trans_plot <- translocation_data %>%
    ggplot(aes(x = value)) +
    geom_density(
      aes(y = after_stat(density)),
      size = 0.8,
      adjust = 1.5,
      alpha = 0.6,
      color = "#fdae61",
      fill = "#fdae61"
    ) +
    theme_minimal(base_size = 12) +
    labs(x = "Particle Length (µm)", y = "Relative Frequency") +
    theme(
      legend.title = element_blank(),
      legend.position = c(0.2, 0.7),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 13),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
    ) +
    ggtitle("Tissue Translocation Size Limit (µm)")

  alpha_combined_plot <- ggarrange(
    alpha_length_plot,
    alpha_mass_plot,
    alpha_volume_plot,
    alpha_surface_area_plot,
    length_width_plot,
    height_width_plot,
    common.legend = T,
    nrow = 3,
    ncol = 2,
    legend = "right"
  )

  tissue_body <- ggarrange(body_length_plot, tissue_trans_plot, ncol = 2)

  return(list(
    alpha_combined_plot = alpha_combined_plot,
    tissue_body = tissue_body
  ))
}
