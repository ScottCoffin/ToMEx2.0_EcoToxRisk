#### Special Plotting Functions #####

#' Visualize sampled parameter distributions
#'
#' Generates histograms for the Latin hypercube samples created by
#' \link[=matrix_function]{matrix_function()}, returning ggplot objects for reuse.
#'
#' @param mat Parameter matrix produced by \link[=matrix_function]{matrix_function()}.
#' @param metrics Optional character vector of metrics to include. Defaults to
#'   all available metrics in `mat`. Supported values include
#'   `"length_alpha"`, `"mass_alpha"`, `"volume_alpha"`, `"surface_area_alpha"`,
#'   `"length_width_ratio"`, `"height_width_ratio"`, `"body_length"`,
#'   `"translocation"`. Shorthand aliases such as `"mass"`, `"volume"`, and
#'   `"surface_area"` are also accepted.
#' @param compartments Optional character vector of compartments to include.
#'   Defaults to surface water (`"Marine"`, `"Freshwater"`). Use
#'   `"Marine Sediment"` and `"Freshwater Sediment"` to include sediment draws.
#'
#' @return A list of ggplot objects summarizing parameter draws.
#' @export
parameter_histograms_function <- function(mat, metrics = NULL, compartments = NULL) {
  # Define colors for marine and freshwater parameters
  env_colors <- c(
    "Marine" = "#1f78b4", # Dark Blue
    "Freshwater" = "#a6cee3", # Light Blue
    "Marine Sediment" = "#33a02c", # Dark Green
    "Freshwater Sediment" = "#b2df8a" # Light Green
  )
  body_colors <- c(
    "Beta (mm)" = "#66c2a5",
    "Intercept (mm)" = "#66c282"
  )

  mat <- as.data.frame(mat)

  metric_registry <- list(
    length_alpha = list(
      type = "env",
      title = "Length Alpha",
      columns = c(
        "alpha.marine" = "Length Alpha (Marine)",
        "alpha.freshwater" = "Length Alpha (Freshwater)",
        "alpha.sediment.marine" = "Length Alpha (Marine Sediment)",
        "alpha.sediment.freshwater" = "Length Alpha (Freshwater Sediment)"
      ),
      show_x_title = FALSE
    ),
    mass_alpha = list(
      type = "env",
      title = "Mass Alpha",
      columns = c(
        "a.m.marine" = "Mass Alpha (Marine)",
        "a.m.freshwater" = "Mass Alpha (Freshwater)",
        "a.m.sediment.marine" = "Mass Alpha (Marine Sediment)",
        "a.m.sediment.freshwater" = "Mass Alpha (Freshwater Sediment)"
      ),
      show_x_title = FALSE
    ),
    volume_alpha = list(
      type = "env",
      title = "Volume Alpha",
      columns = c(
        "a.v.marine" = "Volume Alpha (Marine)",
        "a.v.freshwater" = "Volume Alpha (Freshwater)",
        "a.v.sediment.marine" = "Volume Alpha (Marine Sediment)",
        "a.v.sediment.freshwater" = "Volume Alpha (Freshwater Sediment)"
      ),
      show_x_title = FALSE
    ),
    surface_area_alpha = list(
      type = "env",
      title = "Surface Area Alpha",
      columns = c(
        "a.sa.marine" = "Surface Area Alpha (Marine)",
        "a.sa.freshwater" = "Surface Area Alpha (Freshwater)",
        "a.sa.sediment.marine" = "Surface Area Alpha (Marine Sediment)",
        "a.sa.sediment.freshwater" = "Surface Area Alpha (Freshwater Sediment)"
      ),
      show_x_title = TRUE
    ),
    length_width_ratio = list(
      type = "env",
      title = "Length to Width Ratio",
      columns = c(
        "R.ave.water.marine" = "Length-Width Ratio (Marine)",
        "R.ave.water.freshwater" = "Length-Width Ratio (Freshwater)",
        "R.ave.sediment.marine" = "Length-Width Ratio (Marine Sediment)",
        "R.ave.sediment.freshwater" = "Length-Width Ratio (Freshwater Sediment)"
      ),
      show_x_title = TRUE
    ),
    height_width_ratio = list(
      type = "env",
      title = "Height to Width Ratio",
      columns = c(
        "H_W_ratio.marine" = "Height-Width Ratio (Marine)",
        "H_W_ratio.freshwater" = "Height-Width Ratio (Freshwater)",
        "H_W_ratio.sediment.marine" = "Height-Width Ratio (Marine Sediment)",
        "H_W_ratio.sediment.freshwater" = "Height-Width Ratio (Freshwater Sediment)"
      ),
      show_x_title = TRUE
    ),
    body_length = list(
      type = "body",
      title = "Body Length to Mouth Size Estimation",
      columns = c(
        "sim.beta.log10.body.length" = "log10(Body Length, mm) - beta",
        "sim.body.length.intercept" = "log10(Body Length, mm) - intercept"
      ),
      group_map = c(
        "log10(Body Length, mm) - beta" = "Beta (mm)",
        "log10(Body Length, mm) - intercept" = "Intercept (mm)"
      ),
      show_x_title = TRUE,
      legend_position = c(0.2, 0.7)
    ),
    translocation = list(
      type = "single",
      title = "Tissue Translocation Size Limit (um)",
      columns = c(
        "upper.tissue.trans.size.um" = "Maximum Translocatable Particle Length (um)"
      ),
      show_x_title = TRUE,
      color = "#fdae61",
      fill = "#fdae61",
      x_label = "Particle Length (um)"
    )
  )

  metric_aliases <- c(
    "length_alpha" = "length_alpha",
    "alpha_length" = "length_alpha",
    "mass_alpha" = "mass_alpha",
    "mass" = "mass_alpha",
    "volume_alpha" = "volume_alpha",
    "volume" = "volume_alpha",
    "surface_area_alpha" = "surface_area_alpha",
    "surface_area" = "surface_area_alpha",
    "surfacearea" = "surface_area_alpha",
    "length_width_ratio" = "length_width_ratio",
    "length_width" = "length_width_ratio",
    "height_width_ratio" = "height_width_ratio",
    "height_width" = "height_width_ratio",
    "body_length" = "body_length",
    "body" = "body_length",
    "translocation" = "translocation",
    "tissue_translocation" = "translocation"
  )

  normalize_key <- function(value) {
    gsub("[^a-z0-9]+", "_", tolower(value))
  }

  resolve_metrics <- function(metrics, registry, aliases) {
    if (is.null(metrics)) {
      return(names(registry))
    }
    resolved <- vapply(metrics, function(metric) {
      key <- normalize_key(metric)
      unname(aliases[key])
    }, character(1))

    if (any(is.na(resolved))) {
      invalid <- metrics[is.na(resolved)]
      stop(
        "Unknown metric(s): ",
        paste(invalid, collapse = ", "),
        ". Supported values include: ",
        paste(names(registry), collapse = ", "),
        call. = FALSE
      )
    }

    unique(resolved)
  }

  resolve_compartments <- function(compartments, allowed, default_compartments = allowed) {
    if (is.null(compartments)) {
      return(default_compartments)
    }
    if (!length(compartments)) {
      return(character(0))
    }
    compartment_map <- c(
      "marine" = "Marine",
      "freshwater" = "Freshwater",
      "marine sediment" = "Marine Sediment",
      "marine_sediment" = "Marine Sediment",
      "sediment marine" = "Marine Sediment",
      "sediment_marine" = "Marine Sediment",
      "freshwater sediment" = "Freshwater Sediment",
      "freshwater_sediment" = "Freshwater Sediment",
      "sediment freshwater" = "Freshwater Sediment",
      "sediment_freshwater" = "Freshwater Sediment"
    )
    resolved <- unname(compartment_map[tolower(compartments)])
    if (any(is.na(resolved))) {
      invalid <- compartments[is.na(resolved)]
      stop(
        "Unknown compartment(s): ",
        paste(invalid, collapse = ", "),
        ". Supported values include: ",
        paste(allowed, collapse = ", "),
        call. = FALSE
      )
    }
    unique(resolved)
  }

  detect_compartment <- function(column) {
    if (grepl("sediment", column, ignore.case = TRUE) &&
      grepl("marine", column, ignore.case = TRUE)) {
      return("Marine Sediment")
    }
    if (grepl("sediment", column, ignore.case = TRUE) &&
      grepl("fresh", column, ignore.case = TRUE)) {
      return("Freshwater Sediment")
    }
    if (grepl("marine", column, ignore.case = TRUE)) {
      return("Marine")
    }
    if (grepl("fresh", column, ignore.case = TRUE)) {
      return("Freshwater")
    }
    NA_character_
  }

  build_env_data <- function(mat, columns, compartments) {
    cols <- intersect(names(columns), names(mat))
    if (!length(cols)) {
      return(NULL)
    }
    compartment_map <- vapply(cols, detect_compartment, character(1))
    keep <- cols[!is.na(compartment_map) & compartment_map %in% compartments]
    if (!length(keep)) {
      return(NULL)
    }
    labels <- columns[keep]
    rename_map <- stats::setNames(keep, labels)
    group_map <- stats::setNames(compartment_map[keep], labels)
    data <- mat %>%
      dplyr::select(dplyr::all_of(keep)) %>%
      dplyr::rename(!!!rename_map) %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "name",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        Group = group_map[name]
      ) %>%
      dplyr::filter(!is.na(Group))
    if (!nrow(data)) {
      return(NULL)
    }
    data$Group <- factor(data$Group, levels = unique(data$Group))
    data
  }

  build_group_data <- function(mat, columns, group_map) {
    cols <- intersect(names(columns), names(mat))
    if (!length(cols)) {
      return(NULL)
    }
    labels <- columns[cols]
    rename_map <- stats::setNames(cols, labels)
    group_labels <- group_map[labels]
    data <- mat %>%
      dplyr::select(dplyr::all_of(cols)) %>%
      dplyr::rename(!!!rename_map) %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "name",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        Group = group_labels[name]
      ) %>%
      dplyr::filter(!is.na(Group))
    if (!nrow(data)) {
      return(NULL)
    }
    data$Group <- factor(data$Group, levels = unique(data$Group))
    data
  }

  build_single_data <- function(mat, columns) {
    cols <- intersect(names(columns), names(mat))
    if (!length(cols)) {
      return(NULL)
    }
    labels <- columns[cols]
    rename_map <- stats::setNames(cols, labels)
    data <- mat %>%
      dplyr::select(dplyr::all_of(cols)) %>%
      dplyr::rename(!!!rename_map) %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "name",
        values_to = "value"
      )
    if (!nrow(data)) {
      return(NULL)
    }
    data
  }

  base_density_plot <- function(data, title, fill_colors, show_x_title, legend_position = "right") {
    x_title <- if (isTRUE(show_x_title)) {
      ggplot2::element_text(size = 13)
    } else {
      ggplot2::element_blank()
    }
    ggplot2::ggplot(data, ggplot2::aes(x = value, fill = Group)) +
      ggplot2::geom_density(
        ggplot2::aes(y = ggplot2::after_stat(density)),
        linewidth = 0.8,
        adjust = 1.5,
        alpha = 0.6
      ) +
      ggplot2::scale_fill_manual(
        values = fill_colors,
        breaks = names(fill_colors)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = "Value", y = "Relative Frequency") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
        axis.title.y = ggplot2::element_blank(),
        axis.title.x = x_title,
        axis.text.x = ggplot2::element_text(size = 10),
        axis.text.y = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 12, face = "bold"),
        legend.position = legend_position
      ) +
      ggplot2::ggtitle(title)
  }

  compartments <- resolve_compartments(
    compartments,
    names(env_colors),
    default_compartments = c("Marine", "Freshwater")
  )
  resolved_metrics <- resolve_metrics(metrics, metric_registry, metric_aliases)

  metric_available <- function(definition) {
    cols <- intersect(names(definition$columns), names(mat))
    if (!length(cols)) {
      return(FALSE)
    }
    if (definition$type == "env") {
      compartment_map <- vapply(cols, detect_compartment, character(1))
      return(any(!is.na(compartment_map) & compartment_map %in% compartments))
    }
    TRUE
  }

  available_metrics <- names(metric_registry)[vapply(metric_registry, metric_available, logical(1))]
  selected_metrics <- if (is.null(metrics)) {
    available_metrics
  } else {
    intersect(resolved_metrics, available_metrics)
  }

  build_plot <- function(metric_key) {
    definition <- metric_registry[[metric_key]]
    if (definition$type == "env") {
      data <- build_env_data(mat, definition$columns, compartments)
      if (is.null(data)) {
        return(NULL)
      }
      fill_colors <- env_colors[levels(data$Group)]
      return(
        base_density_plot(
          data = data,
          title = definition$title,
          fill_colors = fill_colors,
          show_x_title = definition$show_x_title
        )
      )
    }
    if (definition$type == "body") {
      data <- build_group_data(mat, definition$columns, definition$group_map)
      if (is.null(data)) {
        return(NULL)
      }
      fill_colors <- body_colors[levels(data$Group)]
      plot <- base_density_plot(
        data = data,
        title = definition$title,
        fill_colors = fill_colors,
        show_x_title = definition$show_x_title,
        legend_position = definition$legend_position
      )
      return(
        plot +
          ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 11)
          )
      )
    }
    if (definition$type == "single") {
      data <- build_single_data(mat, definition$columns)
      if (is.null(data)) {
        return(NULL)
      }
      return(
        ggplot2::ggplot(data, ggplot2::aes(x = value)) +
          ggplot2::geom_density(
            ggplot2::aes(y = ggplot2::after_stat(density)),
            linewidth = 0.8,
            adjust = 1.5,
            alpha = 0.6,
            color = definition$color,
            fill = definition$fill
          ) +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::labs(x = definition$x_label, y = "Relative Frequency") +
          ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.position = c(0.2, 0.7),
            plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
            axis.title.y = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_text(size = 13),
            axis.text.x = ggplot2::element_text(size = 10),
            axis.text.y = ggplot2::element_blank()
          ) +
          ggplot2::ggtitle(definition$title)
      )
    }
    NULL
  }

  plots <- lapply(selected_metrics, build_plot)
  names(plots) <- selected_metrics
  plots <- plots[!vapply(plots, is.null, logical(1))]

  alpha_keys <- intersect(
    c(
      "length_alpha",
      "mass_alpha",
      "volume_alpha",
      "surface_area_alpha",
      "length_width_ratio",
      "height_width_ratio"
    ),
    names(plots)
  )

  alpha_combined_plot <- NULL
  if (length(alpha_keys)) {
    alpha_plots <- plots[alpha_keys]
    alpha_cols <- min(2, length(alpha_plots))
    alpha_rows <- ceiling(length(alpha_plots) / alpha_cols)
    alpha_combined_plot <- ggpubr::ggarrange(
      plotlist = alpha_plots,
      common.legend = TRUE,
      nrow = alpha_rows,
      ncol = alpha_cols,
      legend = "right"
    )
  }

  tissue_keys <- intersect(c("body_length", "translocation"), names(plots))
  tissue_body <- NULL
  if (length(tissue_keys)) {
    tissue_plots <- plots[tissue_keys]
    tissue_body <- ggpubr::ggarrange(
      plotlist = tissue_plots,
      ncol = min(2, length(tissue_plots))
    )
  }

  list(
    alpha_combined_plot = alpha_combined_plot,
    tissue_body = tissue_body
  )
}
