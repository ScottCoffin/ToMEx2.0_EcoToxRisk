##### Measured (laboratory) size distributions of polydisperse test materials #####
#
# align_data() assumes, by default, that a polydisperse test material follows
# the *environmental* power-law length distribution (Kooi et al. 2021) between
# its reported minimum and maximum lengths: no study in ToMEx reports the size
# distribution of its test material (see the SI of Coffin et al. 2026). These
# helpers let a record carry its own measured distribution instead. They are
# used only for records that supply one, so they never change the alignment
# of records that do not (and therefore never change thresholds derived from
# ToMEx as it stands).
#
# Optional per-record input columns (all may be absent or NA):
#   size.dist.family     "powerlaw" or "lognormal" (number-based length distribution)
#   size.dist.alpha      power-law exponent of the length distribution ("powerlaw")
#   size.dist.median.um  median length, um ("lognormal")
#   size.dist.gsd        geometric standard deviation ("lognormal")
#   size.D10.um, size.D50.um, size.D90.um
#                        number-based percentiles; used for "lognormal" when the
#                        median/GSD are not given (D50 = median; GSD from D10/D90)
# The reported minimum and maximum lengths
# (size.length.{min,max}.um.used.for.conversions) truncate the distribution.

#' Number-based length distribution families accepted by align_data()
#' @noRd
.poly_dist_families <- c("powerlaw", "lognormal")

#' Log-normal parameters (meanlog, sdlog) from a median and GSD, or from
#' D10/D50/D90 (vectorised). Percentiles are number-based.
#' @noRd
.poly_dist_lognormal_params <- function(median = NA_real_, gsd = NA_real_,
                                        d10 = NA_real_, d50 = NA_real_,
                                        d90 = NA_real_) {
  n <- max(length(median), length(gsd), length(d10), length(d50), length(d90))
  median <- rep_len(median, n)
  gsd <- rep_len(gsd, n)
  d10 <- rep_len(d10, n)
  d50 <- rep_len(d50, n)
  d90 <- rep_len(d90, n)
  z90 <- stats::qnorm(0.9)
  sd_from_d <- dplyr::case_when(
    !is.na(d10) & !is.na(d90) ~ (log(d90) - log(d10)) / (2 * z90),
    !is.na(d50) & !is.na(d90) ~ (log(d90) - log(d50)) / z90,
    !is.na(d10) & !is.na(d50) ~ (log(d50) - log(d10)) / z90,
    TRUE ~ NA_real_
  )
  meanlog <- dplyr::coalesce(
    log(median), log(d50),
    ifelse(!is.na(d10) & !is.na(d90), (log(d10) + log(d90)) / 2, NA_real_)
  )
  sdlog <- dplyr::coalesce(log(gsd), sd_from_d)
  list(meanlog = meanlog, sdlog = sdlog)
}

#' Unnormalised CDF of the number-based length distribution (scalar inputs).
#' @noRd
.poly_dist_cdf <- function(x, family, alpha, meanlog, sdlog) {
  switch(family,
    powerlaw = if (isTRUE(all.equal(alpha, 1))) log(x) else x^(1 - alpha) / (1 - alpha),
    lognormal = stats::plnorm(x, meanlog, sdlog)
  )
}

#' Number fraction of particles between `lmin` and `lmax` that are no longer
#' than `limit` (vectorised); NA where the distribution is not specified.
#' For a power law this equals CFfnx() (Koelmans et al. 2020, eq. 2).
#' @noRd
.poly_dist_fraction <- function(family, alpha, meanlog, sdlog, lmin, lmax, limit) {
  n <- length(family)
  out <- rep(NA_real_, n)
  if (all(is.na(family))) return(out)
  for (i in seq_len(n)) {
    if (is.na(family[i]) || anyNA(c(lmin[i], lmax[i], limit[i]))) next
    up <- min(lmax[i], limit[i])
    if (up <= lmin[i]) {
      out[i] <- if (limit[i] >= lmin[i]) 1 else 0
      next
    }
    cdf <- function(x) .poly_dist_cdf(x, family[i], alpha[i], meanlog[i], sdlog[i])
    out[i] <- (cdf(up) - cdf(lmin[i])) / (cdf(lmax[i]) - cdf(lmin[i]))
  }
  out
}

#' Width and height (um) of a particle of length L under align_data()'s
#' geometry rules (scalar shape; vectorised over L).
#' @noRd
.poly_dist_dims <- function(L, shape, R, H_W_ratio, fiber_width) {
  switch(shape,
    Sphere = list(width = L, height = L),
    Fiber = {
      w <- pmin(fiber_width, L)
      list(width = w, height = w * H_W_ratio)
    },
    # Fragment, Not Reported
    list(width = R * L, height = R * L * H_W_ratio)
  )
}

#' Mean particle volume ("volume") or surface area ("surface_area") of the
#' particles between `lower` and `upper` length (vectorised), integrating the
#' number-based length distribution with align_data()'s shape geometry.
#' NA where no distribution is given or the shape is not one of
#' "Sphere", "Fiber", "Fragment", "Not Reported".
#' @noRd
.poly_dist_mean <- function(family, alpha, meanlog, sdlog, lower, upper,
                            shape, R, H_W_ratio, fiber_width,
                            metric = c("volume", "surface_area")) {
  metric <- match.arg(metric)
  shape <- as.character(shape)
  n <- length(family)
  out <- rep(NA_real_, n)
  if (all(is.na(family))) return(out)
  known <- c("Sphere", "Fiber", "Fragment", "Not Reported")
  for (i in seq_len(n)) {
    if (is.na(family[i]) || !shape[i] %in% known ||
      anyNA(c(lower[i], upper[i], R[i], H_W_ratio[i]))) {
      next
    }
    g <- function(L) {
      d <- .poly_dist_dims(L, shape[i], R[i], H_W_ratio[i], fiber_width[i])
      if (metric == "volume") {
        volumefnx(length = L, width = d$width, height = d$height)
      } else {
        SAfnx(length = L, width = d$width, height = d$height)
      }
    }
    if (upper[i] <= lower[i]) {
      out[i] <- g(lower[i])
      next
    }
    dens <- switch(family[i],
      powerlaw = function(L) L^(-alpha[i]),
      lognormal = function(L) stats::dlnorm(L, meanlog[i], sdlog[i])
    )
    # integrate on log(L): smooth for both families over several decades
    num <- stats::integrate(function(u) g(exp(u)) * dens(exp(u)) * exp(u),
      log(lower[i]), log(upper[i]), rel.tol = 1e-8, stop.on.error = FALSE)
    den <- stats::integrate(function(u) dens(exp(u)) * exp(u),
      log(lower[i]), log(upper[i]), rel.tol = 1e-8, stop.on.error = FALSE)
    out[i] <- num$value / den$value
  }
  out
}

#' Validate and resolve the optional size-distribution columns of `df`.
#' Returns `df` with resolved columns `.sd_family`, `.sd_alpha`,
#' `.sd_meanlog`, `.sd_sdlog` (NA when no usable distribution).
#' @noRd
.poly_dist_resolve <- function(df) {
  cols <- c(
    "size.dist.family", "size.dist.alpha", "size.dist.median.um",
    "size.dist.gsd", "size.D10.um", "size.D50.um", "size.D90.um"
  )
  for (col in cols) {
    if (!col %in% names(df)) {
      df[[col]] <- if (col == "size.dist.family") NA_character_ else NA_real_
    }
  }
  family <- tolower(trimws(as.character(df$size.dist.family)))
  family[family %in% c("", "na")] <- NA_character_
  bad <- setdiff(unique(family[!is.na(family)]), .poly_dist_families)
  if (length(bad) > 0) {
    stop("size.dist.family must be one of ",
      paste(dQuote(.poly_dist_families, FALSE), collapse = ", "),
      " (or NA); got ", paste(dQuote(bad, FALSE), collapse = ", "),
      call. = FALSE
    )
  }
  ln <- .poly_dist_lognormal_params(
    df$size.dist.median.um, df$size.dist.gsd,
    df$size.D10.um, df$size.D50.um, df$size.D90.um
  )
  poly <- if ("polydispersity" %in% names(df)) {
    as.character(df$polydispersity) %in% "polydisperse"
  } else {
    rep(FALSE, nrow(df))
  }
  usable <- !is.na(family) & poly &
    dplyr::case_when(
      family == "powerlaw" ~ !is.na(df$size.dist.alpha),
      family == "lognormal" ~ !is.na(ln$meanlog) & !is.na(ln$sdlog) & ln$sdlog > 0,
      TRUE ~ FALSE
    )
  incomplete <- !is.na(family) & !usable & poly
  if (any(incomplete)) {
    warning(sum(incomplete), " polydisperse record(s) name a size.dist.family ",
      "without usable parameters; they use the environmental power law.",
      call. = FALSE
    )
  }
  df$.sd_family <- ifelse(usable, family, NA_character_)
  df$.sd_alpha <- df$size.dist.alpha
  df$.sd_meanlog <- ln$meanlog
  df$.sd_sdlog <- ln$sdlog
  df
}
