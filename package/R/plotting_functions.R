#### Special Plotting Functions #####

#' Visualize sampled parameter distributions
#'
#' Generates histograms for the Latin hypercube samples created by
#' \link[=matrix_function]{matrix_function()}, returning ggplot objects for reuse.
#'
#' @param mat Parameter matrix produced by \link[=matrix_function]{matrix_function()}.
#'
#' @return A list of ggplot objects summarizing parameter draws.
#' @export
parameter_histograms_function <- function(mat) {
  # Define colors for marine and freshwater parameters
  env_colors <- c(
    "Marine" = "#1f78b4", # Dark Blue
    "Freshwater" = "#a6cee3" # Light Blue
  )
  body_colors <- c(
    "Beta (mm)" = "#66c2a5",
    "Intercept (mm)" = "#66c282"
  )

  # Prepare data for each group
  alpha_data <- mat %>%
    dplyr::select(alpha.marine, alpha.freshwater) %>%
    dplyr::rename(
      "Length Alpha (Marine)" = "alpha.marine",
      "Length Alpha (Freshwater)" = "alpha.freshwater"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      Environment = factor(
        dplyr::if_else(
          grepl("Marine", name),
          "Marine",
          "Freshwater"
        ),
        levels = names(env_colors)
      )
    )

  mass_data <- mat %>%
    dplyr::select(a.m.marine, a.m.freshwater) %>%
    dplyr::rename(
      "Mass Alpha (Marine)" = "a.m.marine",
      "Mass Alpha (Freshwater)" = "a.m.freshwater"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      Environment = factor(
        dplyr::if_else(
          grepl("Marine", name),
          "Marine",
          "Freshwater"
        ),
        levels = names(env_colors)
      )
    )

  volume_data <- mat %>%
    dplyr::select(a.v.marine, a.v.freshwater) %>%
    dplyr::rename(
      "Volume Alpha (Marine)" = "a.v.marine",
      "Volume Alpha (Freshwater)" = "a.v.freshwater"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      Environment = factor(
        dplyr::if_else(
          grepl("Marine", name),
          "Marine",
          "Freshwater"
        ),
        levels = names(env_colors)
      )
    )

  surface_area_data <- mat %>%
    dplyr::select(a.sa.marine, a.sa.freshwater) %>%
    dplyr::rename(
      "Surface Area Alpha (Marine)" = "a.sa.marine",
      "Surface Area Alpha (Freshwater)" = "a.sa.freshwater"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      Environment = factor(
        dplyr::if_else(
          grepl("Marine", name),
          "Marine",
          "Freshwater"
        ),
        levels = names(env_colors)
      )
    )

  ratio_data <- mat %>%
    dplyr::select(R.ave.water.marine, R.ave.water.freshwater) %>%
    dplyr::rename(
      "Length-Width Ratio (Marine)" = "R.ave.water.marine",
      "Length-Width Ratio (Freshwater)" = "R.ave.water.freshwater"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      Environment = factor(
        dplyr::if_else(
          grepl("Marine", name),
          "Marine",
          "Freshwater"
        ),
        levels = names(env_colors)
      )
    )

  H_W_ratio_data <- mat %>%
    dplyr::select(H_W_ratio.marine, H_W_ratio.freshwater) %>%
    dplyr::rename(
      "Height-Width Ratio (Marine)" = "H_W_ratio.marine",
      "Height-Width Ratio (Freshwater)" = "H_W_ratio.freshwater"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      Environment = factor(
        dplyr::if_else(
          grepl("Marine", name),
          "Marine",
          "Freshwater"
        ),
        levels = names(env_colors)
      )
    )

  body_length_data <- mat %>%
    dplyr::select(sim.beta.log10.body.length, sim.body.length.intercept) %>%
    dplyr::rename(
      "log10(Body Length, mm) - beta" = "sim.beta.log10.body.length",
      "log10(Body Length, mm) - intercept" = "sim.body.length.intercept"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      Environment = factor(
        dplyr::case_when(
          grepl("beta", name) ~ "Beta (mm)",
          grepl("intercept", name) ~ "Intercept (mm)"
        ),
        levels = names(body_colors)
      )
    )

  translocation_data <- mat %>%
    dplyr::select(upper.tissue.trans.size.um) %>%
    dplyr::rename(
      "Maximum Translocatable Particle Length (um)" = "upper.tissue.trans.size.um"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(Environment = "Translocation")

  # alpha length
  alpha_length_plot <- alpha_data %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = Environment)) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      linewidth = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = unname(env_colors),
      breaks = c("Marine", "Freshwater")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Value", y = "Relative Frequency") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold")
    ) +
    ggplot2::ggtitle("Length Alpha")

  alpha_mass_plot <- mass_data %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = Environment)) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      linewidth = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = unname(env_colors),
      breaks = c("Marine", "Freshwater")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Value", y = "Relative Frequency") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold")
    ) +
    ggplot2::ggtitle("Mass Alpha")

  alpha_volume_plot <- volume_data %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = Environment)) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      linewidth = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = unname(env_colors),
      breaks = c("Marine", "Freshwater")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Value", y = "Relative Frequency") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold")
    ) +
    ggplot2::ggtitle("Volume Alpha")

  alpha_surface_area_plot <- surface_area_data %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = Environment)) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      linewidth = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = unname(env_colors),
      breaks = c("Marine", "Freshwater")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Value", y = "Relative Frequency") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle("Surface Area Alpha")

  length_width_plot <- ratio_data %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = Environment)) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      linewidth = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = unname(env_colors),
      breaks = c("Marine", "Freshwater")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Value", y = "Relative Frequency") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold")
    ) +
    ggplot2::ggtitle("Length to Width Ratio")

  height_width_plot <- H_W_ratio_data %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = Environment)) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      linewidth = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = unname(body_colors),
      breaks = c("Beta (mm)", "Intercept (mm)")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Value", y = "Relative Frequency") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold")
    ) +
    ggplot2::ggtitle("Height to Width Ratio")

  body_length_plot <- body_length_data %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = Environment)) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      linewidth = 0.8,
      adjust = 1.5,
      alpha = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = unname(body_colors),
      breaks = c("Beta (mm)", "Intercept (mm)"),
      name = "Parameter"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Value", y = "Relative Frequency") +
    ggplot2::theme(
      legend.position = c(0.2, 0.7),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle("Body Length to Mouth Size Estimation")

  tissue_trans_plot <- translocation_data %>%
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      linewidth = 0.8,
      adjust = 1.5,
      alpha = 0.6,
      color = "#fdae61",
      fill = "#fdae61"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Particle Length (um)", y = "Relative Frequency") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = c(0.2, 0.7),
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle("Tissue Translocation Size Limit (um)")

  alpha_combined_plot <- ggpubr::ggarrange(
    alpha_length_plot,
    alpha_mass_plot,
    alpha_volume_plot,
    alpha_surface_area_plot,
    length_width_plot,
    height_width_plot,
    common.legend = TRUE,
    nrow = 3,
    ncol = 2,
    legend = "right"
  )

  tissue_body <- ggpubr::ggarrange(
    body_length_plot,
    tissue_trans_plot,
    ncol = 2
  )

  list(
    alpha_combined_plot = alpha_combined_plot,
    tissue_body = tissue_body
  )
}
