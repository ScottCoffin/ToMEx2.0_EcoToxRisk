# SSD helper functions ----------------------------------------------------

#' Fit a Tier 1 SSD and return the HC5 lower confidence limit
#'
#' Collapses endpoints to the 25th percentile per species-group combination,
#' fits a model-averaged SSD (log-logistic, log-normal, gamma, log-Gumbel),
#' and returns the lower confidence limit of the HC5 — the most protective
#' estimate suitable for Tier 1 risk assessment.
#'
#' @param filtered.data A data frame containing aligned toxicity data with
#'   columns `dose_new` (aligned dose), `Species`, and `Group`.
#' @param hcxlcl Integer row index into the bootstrapped prediction table
#'   corresponding to the desired hazard concentration percentile (typically 5
#'   for HC5).
#' @param nboot Number of bootstrap iterations for confidence interval
#'   estimation. Default `10`; increase to 1000+ for publication-quality CIs.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{hc5lcl}{Lower confidence limit of the HC5 (NA if < 6 species).}
#'     \item{n_species}{Number of distinct species in the input data.}
#'     \item{n_groups}{Number of distinct taxonomic groups.}
#'     \item{n_datapoints}{Total number of data points.}
#'   }
#'
#' @seealso [SSD_function_t2_3()], [SSD_function_t3_4()], [make_tiered_SSDs()]
#' @noRd
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

#' Fit a Tier 2 SSD and return the HC5 estimate with 95% CI
#'
#' Collapses endpoints to the 25th percentile per species-group combination,
#' fits a model-averaged SSD (log-logistic, log-normal, gamma, log-Gumbel),
#' and returns the HC5 point estimate with 95% bootstrap confidence interval.
#' Equivalent to [SSD_function_t2_3()].
#'
#' @param filtered.data A data frame with columns `dose_new`, `Species`,
#'   `Group`.
#' @param hcx Integer row index into the prediction table for the desired
#'   hazard concentration percentile (typically 5 for HC5).
#' @param nboot Number of bootstrap iterations. Default `10`.
#'
#' @return A named list with elements `hcx_est`, `hcx05cl`, `hcx95cl`,
#'   `n_species`, `n_groups`, `n_datapoints` (all NA if < 6 species).
#'
#' @seealso [SSD_function_t1()], [SSD_function_t2_3()], [SSD_function_t3_4()],
#'   [make_tiered_SSDs()]
#' @noRd
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

#' Fit a Tier 3/4 SSD and return the HCx estimate with 95% CI
#'
#' Filters to organism- and population-level endpoints (`bio_f` in
#' `c("Organism", "Population")`) with `risk.13 != 1`, collapses to the median
#' (50th percentile) per species-group, fits a model-averaged SSD, and returns
#' the HCx estimate with 95% bootstrap CI. Use `hcx = 5` for Tier 3 (HC5) and
#' `hcx = 10` for Tier 4 (HC10).
#'
#' @param filtered.data A data frame with columns `dose_new`, `Species`,
#'   `Group`, `bio_f`, and `risk.13`.
#' @param hcx Integer row index for the hazard concentration (5 for HC5, 10
#'   for HC10).
#' @param nboot Number of bootstrap iterations. Default `10`.
#'
#' @return A named list with elements `hcx_est`, `hcx05cl`, `hcx95cl`,
#'   `n_species`, `n_groups`, `n_datapoints` (all NA if < 6 species after
#'   filtering).
#'
#' @seealso [SSD_function_t1()], [SSD_function_t2_3()], [make_tiered_SSDs()]
#' @noRd
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

#' Fit a Tier 2 SSD (alias for SSD_function_t2)
#'
#' Identical to [SSD_function_t2()]. Provided as a named alias for clarity
#' when working with tiers 2–3 style data (25th-percentile collapse, all
#' endpoints). Use `hcx = 5` for HC5.
#'
#' @inheritParams SSD_function_t2
#' @inherit SSD_function_t2 return
#'
#' @seealso [SSD_function_t1()], [SSD_function_t2()], [SSD_function_t3_4()],
#'   [make_tiered_SSDs()]
#' @noRd
SSD_function_t2_3 <- function(filtered.data, hcx, nboot = 10) {
  SSD_function_t2(filtered.data = filtered.data, hcx = hcx, nboot = nboot)
}

#' Fit a tiered SSD with the tier supplied as a parameter
#'
#' A unified interface to the tier-specific SSD functions. Data collapsing,
#' filtering, and HCx selection all follow the tier definitions from
#' [Mehinto et al. (2022)](https://link.springer.com/article/10.1186/s43591-022-00033-3)
#' and Coffin et al. (2026):
#'
#' | Tier | Data filter | Collapse | Returns |
#' |------|-------------|----------|---------|
#' | 1 | All endpoints | 25th %ile | HC5 lower CL |
#' | 2 | All endpoints | 25th %ile | HC5 est + 95% CI |
#' | 3 | Organism/Population, risk.13 != 1 | Median | HC5 est + 95% CI |
#' | 4 | Organism/Population, risk.13 != 1 | Median | HC10 est + 95% CI |
#'
#' The return value is an S3 object of class `"tiered_SSD"`. In addition to the
#' scalar HCx summaries it carries all data needed to reproduce the SSD plot via
#' `ggplot2::autoplot()`. The internal per-tier helpers (`SSD_function_t1`,
#' `SSD_function_t2_3`, `SSD_function_t3_4`) are preserved in the source for
#' backward compatibility with existing scripts but are no longer exported.
#'
#' @param filtered.data A data frame of aligned toxicity data. Must contain
#'   `dose_new`, `Species`, `Group`. Tiers 3–4 additionally require `bio_f`
#'   and `risk.13`.
#' @param tier Integer 1–4 specifying the tier to apply.
#' @param hcx Integer percentile for the hazard concentration. Defaults to 5
#'   (HC5) for tiers 1–3 and 10 (HC10) for tier 4.
#' @param nboot Number of bootstrap iterations for confidence intervals.
#'   Default `10`; use 1000+ for publication results.
#'
#' @return An object of class `"tiered_SSD"` (a list) containing:
#'   \describe{
#'     \item{hc5lcl}{(Tier 1 only) HC5 lower confidence limit.}
#'     \item{hcx_est, hcx05cl, hcx95cl}{(Tiers 2–4) HC5/HC10 estimate and 95% CI.}
#'     \item{n_species, n_groups, n_datapoints}{Sample-size metadata.}
#'     \item{tier, hcx, nboot}{Parameters used.}
#'     \item{predictions}{Full 100-row `ssdtools::predict()` output (one row per
#'       percent 1–100). Serves as `aoc_pred` for plotting.}
#'     \item{collapsed_data}{Species-level collapsed data frame with `Conc`,
#'       `frac` (Hazen plotting position), `Species`, and `Group` columns.
#'       Serves as `aoc_ssd` for plotting.}
#'     \item{hc_data}{Single-row subset of `predictions` at the chosen HCx
#'       percentile with an added `est_format` column (scientific notation).
#'       Serves as `aoc_hc` for annotation.}
#'     \item{fit_dists}{The `ssdtools` distribution-fit object. Useful for
#'       re-predicting with different settings.}
#'   }
#'   All plot-data fields are `NULL` when fewer than 6 species remain after
#'   tier-specific filtering.
#'
#' @examples
#' \dontrun{
#' aligned <- align_data(tomex2)
#' fd_water <- aligned |>
#'   dplyr::filter(
#'     environment %in% c("Freshwater", "Marine"),
#'     ingestible != "not ingestible",
#'     EC_env_v.particles.mL_ingest > 0, Group != "Algae",
#'     shape_f != "Not Reported", poly_f != "Not Reported"
#'   ) |>
#'   dplyr::mutate(
#'     dose_new = (EC_env_v.particles.mL_ingest / (af.time * af.noec)) * 1000
#'   ) |>
#'   tidyr::drop_na(dose_new) |>
#'   dplyr::filter(dose_new > 0)
#'
#' result <- make_tiered_SSDs(fd_water, tier = 3)
#' print(result)
#' ggplot2::autoplot(result, x_label = "particles/L", erm_label = "Food Dilution")
#' }
#'
#' @seealso [align_data()], [autoplot.tiered_SSD()], [print.tiered_SSD()]
#' @importFrom ggplot2 autoplot
#' @export
make_tiered_SSDs <- function(filtered.data, tier, hcx = NULL, nboot = 10) {
  tier <- as.integer(tier)
  if (!tier %in% 1L:4L) stop("`tier` must be 1, 2, 3, or 4.")
  if (is.null(hcx)) hcx <- if (tier == 4L) 10L else 5L

  if (tier %in% 3L:4L) {
    work <- filtered.data %>%
      dplyr::filter(.data$risk.13 != 1, .data$bio_f %in% c("Organism", "Population"))
    collapsed <- work %>%
      dplyr::group_by(.data$Species, .data$Group) %>%
      dplyr::summarize(Conc = stats::quantile(.data$dose_new, 0.50), .groups = "drop")
    metadata <- work %>%
      dplyr::summarize(
        n_species    = dplyr::n_distinct(.data$Species),
        n_groups     = dplyr::n_distinct(.data$Group),
        n_datapoints = dplyr::n(),
        .groups = "drop"
      )
  } else {
    collapsed <- filtered.data %>%
      dplyr::group_by(.data$Species, .data$Group) %>%
      dplyr::summarize(Conc = stats::quantile(.data$dose_new, 0.25), .groups = "drop")
    metadata <- filtered.data %>%
      dplyr::summarize(
        n_species    = dplyr::n_distinct(.data$Species),
        n_groups     = dplyr::n_distinct(.data$Group),
        n_datapoints = dplyr::n(),
        .groups = "drop"
      )
  }

  meta <- list(
    n_species    = metadata$n_species,
    n_groups     = metadata$n_groups,
    n_datapoints = metadata$n_datapoints,
    tier         = tier,
    hcx          = hcx,
    nboot        = nboot
  )

  if (nrow(collapsed) < 6) {
    scalar <- if (tier == 1L) list(hc5lcl = NA) else list(hcx_est = NA, hcx05cl = NA, hcx95cl = NA)
    result <- c(scalar, meta, list(predictions = NULL, collapsed_data = NULL, hc_data = NULL, fit_dists = NULL))
    class(result) <- "tiered_SSD"
    return(result)
  }

  # Empirical CDF (Hazen plotting positions) on sorted collapsed data
  collapsed_ssd <- collapsed[order(collapsed$Conc), ]
  n_ssd <- nrow(collapsed_ssd)
  collapsed_ssd$frac <- (seq_len(n_ssd) - 0.5) / n_ssd

  dists <- ssdtools::ssd_fit_dists(
    collapsed,
    left       = "Conc",
    dists      = c("llogis", "lnorm", "gamma", "lgumbel"),
    computable = FALSE,
    silent     = FALSE
  )
  preds <- stats::predict(dists, average = TRUE, nboot = nboot, ci = TRUE)
  if (!"dist" %in% names(preds)) preds$dist <- "average"

  # HC row for plot annotation (row hcx = hcx-th percentile; proportion col = hcx/100)
  hc_data <- preds[hcx, , drop = FALSE]
  hc_data$est_format <- formatC(hc_data$est, format = "e", digits = 2)

  if (tier == 1L) {
    scalar <- list(hc5lcl = c(preds$lcl[hcx]))
  } else {
    scalar <- list(
      hcx_est = c(preds$est[hcx]),
      hcx05cl = c(preds$lcl[hcx]),
      hcx95cl = c(preds$ucl[hcx])
    )
  }

  result <- c(scalar, meta, list(
    predictions    = preds,
    collapsed_data = collapsed_ssd,
    hc_data        = hc_data,
    fit_dists      = dists
  ))
  class(result) <- "tiered_SSD"
  result
}

#' Print method for tiered_SSD objects
#'
#' Displays a compact summary of the SSD fit: tier, HCx, species count, and the
#' hazard concentration estimate(s).
#'
#' @param x A `tiered_SSD` object returned by [make_tiered_SSDs()].
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
print.tiered_SSD <- function(x, ...) {
  cat(sprintf("tiered_SSD  |  Tier %d  |  HC%d\n", x$tier, x$hcx))
  cat(sprintf("Species: %d   Groups: %d   Datapoints: %d\n",
              x$n_species, x$n_groups, x$n_datapoints))
  if (x$tier == 1L) {
    cat(sprintf("HC%d lower CL: %s\n", x$hcx,
                if (is.na(x$hc5lcl)) "NA (< 6 species)" else formatC(x$hc5lcl, format = "e", digits = 3)))
  } else {
    if (is.na(x$hcx_est)) {
      cat(sprintf("HC%d: NA (< 6 species after tier filtering)\n", x$hcx))
    } else {
      cat(sprintf("HC%d: %s  [95%% CI: %s, %s]\n", x$hcx,
                  formatC(x$hcx_est,  format = "e", digits = 3),
                  formatC(x$hcx05cl, format = "e", digits = 3),
                  formatC(x$hcx95cl, format = "e", digits = 3)))
    }
  }
  invisible(x)
}

#' ggplot2 autoplot method for tiered_SSD objects
#'
#' Produces a publication-ready Species Sensitivity Distribution plot from a
#' `tiered_SSD` object. The plot combines the model-averaged SSD curve and
#' confidence ribbon (`predictions`), empirical species points with repelled
#' labels (`collapsed_data`), and crosshair lines with a text annotation
#' pointing to the actual HCx position on the curve (`hc_data`).
#'
#' @param object A `tiered_SSD` object returned by [make_tiered_SSDs()].
#' @param x_label Character string for the x-axis label (dose metric and
#'   units, e.g. `"particles/L"`). Default `"Concentration"`.
#' @param erm_label Optional character string naming the ERM (e.g.
#'   `"Food Dilution"`). Appended to the plot subtitle. Default `NULL`.
#' @param base_size Base font size (pt) passed to `theme_bw()`. Controls
#'   overall text size and is the primary way to scale the plot for different
#'   output dimensions. Default `12`.
#' @param point_size Size of the species data points. Default `2`.
#' @param label_size Size of the repelled species labels. Default `3.5`.
#' @param ... Additional arguments passed to [ggplot2::theme()].
#'
#' @return A `ggplot` object. Returns `NULL` (with a warning) when fewer than 6
#'   species were available and no SSD could be fit.
#'
#' @examples
#' \dontrun{
#' result <- make_tiered_SSDs(fd_water, tier = 3)
#' ggplot2::autoplot(result, x_label = "particles/L", erm_label = "Food Dilution")
#' # Larger base size for a wider figure
#' ggplot2::autoplot(result, x_label = "particles/L", base_size = 14)
#' }
#'
#' @seealso [plotly_ssd()] for an interactive version, [make_tiered_SSDs()]
#' @importFrom ggplot2 autoplot ggplot aes geom_line geom_point geom_vline
#'   geom_hline annotate coord_trans scale_x_continuous scale_y_continuous
#'   labs xlab theme element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales trans_breaks trans_format math_format percent_format
#' @export
autoplot.tiered_SSD <- function(object,
                                x_label    = "Concentration",
                                erm_label  = NULL,
                                base_size  = 12,
                                point_size = 2,
                                label_size = 3.5,
                                ...) {
  if (is.null(object$predictions)) {
    warning("Not enough species (< 6) to produce SSD plot.")
    return(NULL)
  }

  aoc_pred <- object$predictions
  aoc_ssd  <- object$collapsed_data
  aoc_hc   <- object$hc_data

  subtitle <- if (!is.null(erm_label)) {
    paste0("(ERM = ", erm_label, ")  |  Tier ", object$tier)
  } else {
    paste0("Tier ", object$tier)
  }

  dist_label <- paste0("Distribution: ", aoc_hc$dist[1])

  # ssdtools 0.3.7 uses 'proportion' (0.01-0.99); older versions used 'percent' (1-99)
  y_col       <- if ("proportion" %in% names(aoc_pred)) "proportion" else "percent"
  pct_divisor <- if (y_col == "percent") 100 else 1

  # HC y-position (fraction, [0,1]) and annotation label
  hc_y      <- aoc_hc[[y_col]][1] / pct_divisor
  hc_x      <- aoc_hc$est[1]
  hc_annot  <- paste0("HC", object$hcx, ": ", aoc_hc$est_format, "\n", x_label)

  ggplot2::ggplot(aoc_pred, ggplot2::aes(x = .data$est)) +
    ssdtools::geom_xribbon(
      ggplot2::aes(xmin = .data$lcl, xmax = .data$ucl,
                   y = .data[[y_col]] / pct_divisor),
      alpha = 0.2, color = "grey"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[[y_col]] / pct_divisor), color = "gray"
    ) +
    ggplot2::geom_point(
      data = aoc_ssd,
      ggplot2::aes(x = .data$Conc, y = .data$frac, color = .data$Group),
      size = point_size
    ) +
    ggrepel::geom_text_repel(
      data          = aoc_ssd,
      ggplot2::aes(x = .data$Conc, y = .data$frac,
                   label = .data$Species, color = .data$Group),
      nudge_x       = 0.2,
      size          = label_size,
      segment.alpha = 0.5,
      max.overlaps  = Inf
    ) +
    # Crosshairs pointing to the actual HCx position on the curve
    ggplot2::geom_vline(xintercept = hc_x, linetype = "dashed",
                        color = "red", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = hc_y, linetype = "dashed",
                        color = "red", linewidth = 0.5) +
    ggplot2::annotate("text",
      x = hc_x, y = hc_y,
      hjust = -0.12, vjust = -0.4,
      label = hc_annot,
      color = "red", size = label_size * 0.85
    ) +
    ggplot2::scale_y_continuous(
      "Species Affected (%)",
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::xlab(x_label) +
    ggplot2::labs(
      title    = "Microplastics Species Sensitivity Distribution",
      subtitle = subtitle,
      caption  = dist_label
    ) +
    ggplot2::coord_trans(x = "log10") +
    ggplot2::scale_x_continuous(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 15),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      ...
    )
}

#' Interactive plotly SSD for tiered_SSD objects
#'
#' Produces an interactive Plotly version of the Species Sensitivity
#' Distribution plot. The layout mirrors [autoplot.tiered_SSD()] but uses a
#' log-scale x-axis with hover tooltips and an arrow annotation that points
#' exactly to the HCx position on the SSD curve.
#'
#' Requires the `plotly` package (`install.packages("plotly")`).
#'
#' @param object A `tiered_SSD` object returned by [make_tiered_SSDs()].
#' @param x_label Character string for the x-axis and hover label (dose metric
#'   and units, e.g. `"particles/L"`). Default `"Concentration"`.
#' @param erm_label Optional ERM name shown in the plot title (e.g.
#'   `"Food Dilution"`). Default `NULL`.
#' @param lower_size Optional lower particle-size limit (µm) for the x-axis
#'   title (e.g. `1`). Only used when `upper_size` is also supplied.
#' @param upper_size Optional upper particle-size limit (µm) for the x-axis
#'   title (e.g. `5000`). Formatted with comma separator.
#'
#' @return A `plotly` htmlwidget. Returns `NULL` (with a warning) when fewer
#'   than 6 species were available.
#'
#' @examples
#' \dontrun{
#' result <- make_tiered_SSDs(fd_water, tier = 3)
#' plotly_ssd(result, x_label = "particles/L", erm_label = "Food Dilution")
#' plotly_ssd(result, x_label = "particles/L",
#'            lower_size = 1, upper_size = 5000)
#' }
#'
#' @seealso [autoplot.tiered_SSD()] for the static ggplot2 version,
#'   [make_tiered_SSDs()]
#' @export
plotly_ssd <- function(object,
                       x_label    = "Concentration",
                       erm_label  = NULL,
                       lower_size = NULL,
                       upper_size = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  }
  if (is.null(object$predictions)) {
    warning("Not enough species (< 6) to produce SSD plot.")
    return(NULL)
  }

  aoc_pred <- object$predictions
  aoc_ssd  <- object$collapsed_data
  aoc_hc   <- object$hc_data

  # Normalise to proportion [0,1] (ssdtools 0.3.7 = 'proportion'; older = 'percent')
  y_col       <- if ("proportion" %in% names(aoc_pred)) "proportion" else "percent"
  pct_divisor <- if (y_col == "percent") 100 else 1

  aoc_pred$y_prop <- aoc_pred[[y_col]] / pct_divisor
  hc_y            <- aoc_hc[[y_col]][1] / pct_divisor
  hc_x            <- aoc_hc$est[1]

  # CI ribbon as a closed polygon (sorted ascending, then reversed for UCL leg)
  sorted   <- aoc_pred[order(aoc_pred$y_prop), ]
  ribbon_x <- c(sorted$lcl, rev(sorted$ucl))
  ribbon_y <- c(sorted$y_prop, rev(sorted$y_prop))

  # Pre-compute a named color vector that scales to any number of groups.
  # Passing a named hex vector to plotly's `colors` argument avoids the
  # internal RColorBrewer::brewer.pal call, which warns when N > 8 (Set2 max).
  groups       <- sort(unique(aoc_ssd$Group))
  n_groups     <- length(groups)
  group_colors <- stats::setNames(
    grDevices::colorRampPalette(c(
      "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
      "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
    ))(n_groups),
    groups
  )
  # Per-point color vector for the text-labels trace.
  # Using `color = ~Group` on a text-mode trace still triggers brewer.pal
  # internally even when `colors` is a named vector; bypassed via textfont.
  aoc_ssd$label_color <- unname(group_colors[as.character(aoc_ssd$Group)])

  # Build hover text strings (pre-computed to avoid formula-environment issues)
  aoc_ssd$hover_text <- paste0(
    "Species: ", aoc_ssd$Species, "<br>",
    "Group: ", aoc_ssd$Group, "<br>",
    "Concentration: ", signif(aoc_ssd$Conc, 3), " ", x_label, "<br>",
    "% Species Affected: ", signif(aoc_ssd$frac * 100, 2), "%"
  )
  aoc_pred$hover_text <- paste0(
    "% Species Affected: ", signif(aoc_pred$y_prop * 100, 2), "%<br>",
    "Predicted Concentration: ", signif(aoc_pred$est, 3), " ", x_label
  )

  # x-axis title, optionally including particle size range
  x_title <- if (!is.null(lower_size) && !is.null(upper_size)) {
    paste0(x_label, " (", lower_size, " to ",
           formatC(upper_size, big.mark = ",", format = "fg", digits = 3), " µm)")
  } else {
    x_label
  }

  title_text <- paste0(
    "Microplastics Species Sensitivity Distribution",
    if (!is.null(erm_label))
      paste0("<br><span style='font-size:14px'>(ERM = ", erm_label, ")</span>")
    else ""
  )

  hc_annot_text <- paste0(
    "HC", object$hcx, " (", object$hcx, "% Hazard Confidence Level)<br>",
    aoc_hc$est_format, " ", x_label
  )

  # Build figure with explicit sequential assignments — avoids pipe/Reduce
  # lazy-eval issues where plotly captures variable *names* and tries to
  # re-evaluate them after the function returns, producing blank JS.
  fig <- plotly::plot_ly()

  # CI ribbon (closed polygon)
  fig <- plotly::add_trace(fig,
    name      = "95% CI",
    showlegend = FALSE,
    type      = "scatter",
    mode      = "lines",
    x         = ribbon_x,
    y         = ribbon_y,
    fill      = "toself",
    fillcolor = "rgba(180,180,180,0.3)",
    line      = list(color = "rgba(180,180,180,0)"),
    hoverinfo = "none"
  )

  # Model-averaged SSD curve
  fig <- plotly::add_trace(fig,
    name      = "Model-averaged SSD",
    showlegend = FALSE,
    type      = "scatter",
    mode      = "lines",
    x         = aoc_pred$est,
    y         = aoc_pred$y_prop,
    line      = list(color = "black"),
    text      = aoc_pred$hover_text,
    hoverinfo = "text"
  )

  # One trace per Group — bypasses plotly's brewer.pal color machinery
  for (grp in groups) {
    d   <- aoc_ssd[aoc_ssd$Group == grp, ]
    col <- unname(group_colors[grp])
    fig <- plotly::add_trace(fig,
      type          = "scatter",
      mode          = "markers+text",
      name          = grp,
      x             = d$Conc,
      y             = d$frac,
      marker        = list(color = col, size = 8),
      text          = d$Species,
      textposition  = "top right",
      textfont      = list(color = col, size = 11),
      customdata    = d$hover_text,
      hovertemplate = "%{customdata}<extra></extra>"
    )
  }

  # plotly uses log10 coordinate space for annotations on log-type axes
  # (i.e. x must be log10(actual_value), not the actual value itself)
  fig <- plotly::layout(fig,
    title  = list(text = title_text, x = 0.5, xanchor = "center"),
    margin = list(t = 110),
    xaxis  = list(
      type           = "log",
      exponentformat = "e",
      title          = x_title
    ),
    yaxis  = list(
      title      = "Species Affected (%)",
      range      = c(0, 1.05),
      tickformat = ".0%"
    ),
    annotations = list(list(
      x          = log10(hc_x),
      y          = hc_y,
      xref       = "x",
      yref       = "y",
      text       = hc_annot_text,
      showarrow  = TRUE,
      arrowhead  = 7,
      arrowcolor = "red",
      ay         = -80,
      font       = list(color = "red", size = 12)
    ))
  )

  fig
}

#' Apply all four SSD tiers to one environment/ERM subset and collect thresholds.
#' @noRd
process_environment_data <- function(data,
                                     env_filter = "Marine",
                                     upper.tissue.trans.size.um = 88,
                                     x1D_set = 1,
                                     x2D_set = 5000,
                                     nboot = 10) {
  is_sediment <- any(grepl("sediment", env_filter, ignore.case = TRUE))
  particles_ox_col <- if (is_sediment) {
    "particles.kg.ox.stress"
  } else {
    "particles.mL.ox.stress"
  }
  particles_food_col <- if (is_sediment) {
    "particles.kg.food.dilution"
  } else {
    "particles.mL.food.dilution"
  }
  unit_multiplier <- if (is_sediment) 1 else 1000

  if (!particles_ox_col %in% names(data)) {
    if (is_sediment && "EC_env_sa.particles.kg_trans" %in% names(data)) {
      data <- data %>%
        dplyr::mutate(particles.kg.ox.stress = EC_env_sa.particles.kg_trans)
    } else if (!is_sediment && "EC_env_sa.particles.mL_trans" %in% names(data)) {
      data <- data %>%
        dplyr::mutate(particles.mL.ox.stress = EC_env_sa.particles.mL_trans)
    }
  }
  if (!particles_food_col %in% names(data)) {
    if (is_sediment && "EC_env_v.particles.kg_ingest" %in% names(data)) {
      data <- data %>%
        dplyr::mutate(particles.kg.food.dilution = EC_env_v.particles.kg_ingest)
    } else if (!is_sediment && "EC_env_v.particles.mL_ingest" %in% names(data)) {
      data <- data %>%
        dplyr::mutate(particles.mL.food.dilution = EC_env_v.particles.mL_ingest)
    }
  }

  filtered_data_small_default_t1_2 <- data %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dose_new = .data[[particles_ox_col]] / (af.time * af.noec)
    ) %>%
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
    dplyr::mutate(
      dose_new = .data[[particles_food_col]] / (af.time * af.noec)
    ) %>%
    tidyr::drop_na(dose_new) %>%
    dplyr::filter(dose_new > 0) %>%
    dplyr::mutate(dose_new = dose_new * unit_multiplier) %>%
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
