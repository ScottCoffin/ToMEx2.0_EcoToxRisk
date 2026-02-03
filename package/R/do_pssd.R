# -------------------------------------------------------------------------------------------------
# New PSSD++ Function from Coffin et al. (2025), "A Probabilistic Risk Framework for Microplastics Integrating
# Uncertainty Across Toxicological and Environmental Variability: Development and Application to Marine and Freshwater Ecosystems"
# # Arguments:   - DP : matrix of data points
#              - DP.SD: matrix of data point-level standard deviations based on alignments.
#              - UFt : matrix of uncertainty factors for the exposure time
#              - UFdd : matrix of uncertainty factors for the dose-descriptor
#              - SIM : number of iterations in the simulaation
#              - CV.DP : coefficient of variation for the interlaboratory variation
#              - CV.UF : coefficient of variation for the use of non-substance-specific
#                uncertainty factors
#              - rmore_method: method for generating random samples ("step" i.e., trapezoidal; or "lognormal" i.e., short-cut method)
#              - rmore() function is trapezoidal distribution function in rmore.r from the original Wigger et al. (2019) publication

# Helper to align unnamed args to expected positions
resolve_pssd_args <- function(args, arg_order) {
  if (length(args) == 0) {
    return(args)
  }

  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  unnamed_idx <- which(arg_names == "")
  if (length(unnamed_idx) == 0) {
    return(args)
  }

  unnamed <- args[unnamed_idx]
  args <- args[-unnamed_idx]

  for (value in unnamed) {
    current_names <- names(args)
    if (is.null(current_names)) {
      current_names <- character(0)
    }
    next_name <- setdiff(arg_order, current_names)[1]
    if (is.na(next_name)) {
      break
    }
    args[[next_name]] <- value
  }

  args
}

#' Generate probabilistic species sensitivity distributions (PSSD++)
#'
#' Propagates intra- and inter-laboratory uncertainty to sample NOEC
#' distributions per species with either stepwise or log-normal
#' bootstrapping across endpoints.
#'
#' @param .data Matrix of endpoint values (rows = endpoints, cols = species). When
#'   provided, it is treated as `DP`.
#' @param ... Additional arguments used by the sampler. Expected arguments include
#'   `DP`, `DP.SD`, `UFt`, `UFdd`, `SIM`, `CV.DP`, `CV.UF`, and `rmore_method`.
#'
#' @return Matrix of simulated NOEC draws with species in rows and iterations in columns.
#' @export
do.pSSD_mod <- function(.data = NULL, ...) {
  args <- list(...)
  if (!missing(.data) && (is.null(args$DP) || !"DP" %in% names(args))) {
    args <- c(list(DP = .data), args)
  }

  arg_order <- c("DP", "DP.SD", "UFt", "UFdd", "SIM", "CV.DP", "CV.UF", "rmore_method")
  args <- resolve_pssd_args(args, arg_order)

  if (any(names(args) == "")) {
    stop("Too many unnamed arguments supplied to do.pSSD_mod().", call. = FALSE)
  }

  allowed <- arg_order
  unknown <- setdiff(names(args), allowed)
  if (length(unknown) > 0) {
    stop(
      "Unknown arguments in do.pSSD_mod(): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(args$rmore_method)) {
    args$rmore_method <- "step"
  }

  required <- c("DP", "DP.SD", "UFt", "UFdd", "SIM", "CV.DP", "CV.UF")
  missing_args <- setdiff(required, names(args))
  if (length(missing_args) > 0) {
    stop(
      "Missing required arguments in do.pSSD_mod(): ",
      paste(missing_args, collapse = ", "),
      call. = FALSE
    )
  }

  do.call(do_pSSD_mod_internal, args)
}

do_pSSD_mod_internal <- function(
  DP,
  DP.SD,
  UFt,
  UFdd,
  SIM,
  CV.DP,
  CV.UF,
  rmore_method = "step"
) {
  # Coerce to 2-D matrices to avoid apply() errors on vectors/1-D objects
  if (!is.null(DP) && (is.null(dim(DP)) || length(dim(DP)) < 2)) {
    DP <- matrix(DP, ncol = 1)
  }
  if (!is.null(DP.SD) && (is.null(dim(DP.SD)) || length(dim(DP.SD)) < 2)) {
    DP.SD <- matrix(DP.SD, ncol = 1)
  }
  if (!is.null(UFt) && (is.null(dim(UFt)) || length(dim(UFt)) < 2)) {
    UFt <- matrix(UFt, ncol = 1)
  }
  if (!is.null(UFdd) && (is.null(dim(UFdd)) || length(dim(UFdd)) < 2)) {
    UFdd <- matrix(UFdd, ncol = 1)
  }

  if (is.null(DP) || nrow(DP) == 0 || ncol(DP) == 0) {
    warning("No species with data after filtering; returning empty PSSD.")
    return(matrix(NA_real_, nrow = 0, ncol = SIM))
  }

  # Ensure trapezoid helpers are in scope for mc2d::rtrunc()
  rtrapezoid <- trapezoid::rtrapezoid
  ptrapezoid <- trapezoid::ptrapezoid
  qtrapezoid <- trapezoid::qtrapezoid
  dtrapezoid <- trapezoid::dtrapezoid

  #step (slow) or lognormal (fast)
  # Check for species with no data
  if (any(apply(DP, 2, function(x) length(which(!is.na(x)))) == 0)) {
    warning(
      "No data is available for one or more species, it/they won't contribute to the PSSD calculation."
    )
    ind.sp.rem <- which(apply(DP, 2, function(x) length(which(!is.na(x)))) == 0)
    DP <- DP[, -ind.sp.rem, drop = FALSE]
    if (!is.null(DP.SD)) {
      DP.SD <- DP.SD[, -ind.sp.rem, drop = FALSE]
    }
    UFt <- UFt[, -ind.sp.rem, drop = FALSE]
    UFdd <- UFdd[, -ind.sp.rem, drop = FALSE]
  }

  if (is.null(dim(DP)) || ncol(DP) == 0) {
    warning("No species with data after filtering; returning empty PSSD.")
    return(matrix(NA_real_, nrow = 0, ncol = SIM))
  }

  # Calculate corrected endpoints
  corr.endpoints <- DP / (UFdd * UFt)
  sort.endpoints <- apply(corr.endpoints, 2, sort)

  # Initialize matrix for results
  NOEC_comb <- matrix(NA, ncol(DP), SIM, dimnames = list(colnames(DP), NULL))

  for (sp in colnames(DP)) {
    # Identify min and max indices
    ind.min <- which.min(corr.endpoints[, sp])
    ind.max <- which.max(corr.endpoints[, sp])

    # Handle CV.DP as a matrix or scalar
    CV.DP2 <- if (is.matrix(CV.DP)) CV.DP[, sp] else CV.DP

    # Calculate sp.min and sp.max with explicit non-negative checks
    sp.min <- max(
      0,
      corr.endpoints[ind.min, sp] *
        (1 - (sqrt(sum((CV.DP2 / 2.45)^2) + 2 * (CV.UF / 2.45)^2) * 2.45))
    )
    sp.max <- max(
      0,
      corr.endpoints[ind.max, sp] *
        (1 + (sqrt(sum((CV.DP2 / 2.45)^2) + 2 * (CV.UF / 2.45)^2) * 2.45))
    )

    # Handle cases based on the number of unique endpoints
    if (length(unique(sort.endpoints[[sp]])) == 1) {
      # Single endpoint: Use truncated triangular distribution
      NOEC_comb[sp, ] <- mc2d::rtrunc(
        "rtriang",
        min = sp.min,
        mode = sort.endpoints[[sp]][1],
        max = sp.max,
        n = SIM,
        linf = 0
      )
    } else if (length(sort.endpoints[[sp]]) == 2) {
      # Two endpoints: Use truncated trapezoidal distribution
      NOEC_comb[sp, ] <- mc2d::rtrunc(
        "rtrapezoid",
        SIM,
        mode1 = sort.endpoints[[sp]][1],
        mode2 = sort.endpoints[[sp]][2],
        min = sp.min,
        max = sp.max,
        linf = 0
      )

      ### Three endpoints: Use rmore function based on the specified method  ##
      # (step-wise trapezoidal - i.e., original Wigger et al. (2020))
    } else if (length(sort.endpoints[[sp]]) > 2 & rmore_method == "step") {
      # Three or more endpoints: Use step distribution
      NOEC_comb[sp, ] <- rmore(
        values = sort.endpoints[[sp]],
        max = sp.max,
        min = sp.min,
        N = SIM,
        linf = 0
      )

      # Log10 -> normal distribution bootstrapping shortcut method developed here
    } else if (length(sort.endpoints[[sp]]) > 2 & rmore_method == "lognormal") {
      # Three or more endpoints: Log-Normal distribution bootstrapping with uncertainty factor
      #Removing the DP.SD[, sp] from the groups to test
      high <- (sqrt(sum(
        c(mean(DP.SD[, sp] / DP[, sp], na.rm = T), 1 + CV.DP2, 1 + CV.UF)^2,
        na.rm = T
      )))
      low <- 1 / high

      uncertainty_factor <- stats::runif(min = low, max = high, n = SIM)

      data <- 10^stats::rnorm(
        mean = mean(log10(sort.endpoints[[sp]])),
        sd = stats::sd(log10(sort.endpoints[[sp]])),
        n = SIM
      )

      # Ensure no negative values in the final result
      NOEC_comb[sp, ] <- data * uncertainty_factor
    }
  }

  # Return the final matrix
  return(NOEC_comb)
}


########################## Original Function from https://zenodo.org/records/3267194 #########
# Function for generating NOEC distributions for each species
#
# Yields a matrix with dimensions number of species * number of iterations in the simulation
#
# Arguments:   - DP : matrix of data points
#              - DP.SD: matrix of data point-level standard deviations based on alignments.
#              - UFt : matrix of uncertainty factors for the exposure time
#              - UFdd : matrix of uncertainty factors for the dose-descriptor
#              - SIM : number of iterations in the simulation
#              - CV.DP : coefficient of variation for the interlaboratory variation
#              - CV.UF : coefficient of variation for the use of non-substance-specific
#                uncertainty factors
#
# Date of last modification: 03.07.2019
#
# Associated publication: Systematic consideration of parameter uncertainty and variability in
#                         probabilistic species sensitivity distributions
#
# Authors: Henning Wigger, Delphine Kawecki, Bernd Nowack and Vronique Adam
#
# Institute: Empa, Swiss Federal Laboratories for Materials Science and Technology,
#            Technology and Society Laboratory, Lerchenfeldstrasse 5, 9014 St. Gallen, Switzerland
#
# submitted to Integrated Environmental Assessment and Management in July 2019
#' Probabilistic SSD generator (original implementation)
#'
#' Samples NOEC distributions per species using triangular, trapezoidal, or
#' stepwise sampling depending on available endpoints.
#'
#' @param .data Matrix of endpoint values (rows = endpoints, cols = species). When
#'   provided, it is treated as `DP`.
#' @param ... Additional arguments used by the sampler. Expected arguments include
#'   `DP`, `UFt`, `UFdd`, `SIM`, `CV.DP`, and `CV.UF`.
#'
#' @return Matrix of simulated NOEC draws with species in rows and iterations in columns.
#' @export
do.pSSD <- function(.data = NULL, ...) {
  args <- list(...)
  if (!missing(.data) && (is.null(args$DP) || !"DP" %in% names(args))) {
    args <- c(list(DP = .data), args)
  }

  arg_order <- c("DP", "UFt", "UFdd", "SIM", "CV.DP", "CV.UF")
  args <- resolve_pssd_args(args, arg_order)

  if (any(names(args) == "")) {
    stop("Too many unnamed arguments supplied to do.pSSD().", call. = FALSE)
  }

  allowed <- arg_order
  unknown <- setdiff(names(args), allowed)
  if (length(unknown) > 0) {
    stop(
      "Unknown arguments in do.pSSD(): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  required <- c("DP", "UFt", "UFdd", "SIM", "CV.DP", "CV.UF")
  missing_args <- setdiff(required, names(args))
  if (length(missing_args) > 0) {
    stop(
      "Missing required arguments in do.pSSD(): ",
      paste(missing_args, collapse = ", "),
      call. = FALSE
    )
  }

  do.call(do_pSSD_internal, args)
}

do_pSSD_internal <- function(DP, UFt, UFdd, SIM, CV.DP, CV.UF) {
  if (!is.null(DP) && (is.null(dim(DP)) || length(dim(DP)) < 2)) {
    DP <- matrix(DP, ncol = 1)
  }
  if (!is.null(UFt) && (is.null(dim(UFt)) || length(dim(UFt)) < 2)) {
    UFt <- matrix(UFt, ncol = 1)
  }
  if (!is.null(UFdd) && (is.null(dim(UFdd)) || length(dim(UFdd)) < 2)) {
    UFdd <- matrix(UFdd, ncol = 1)
  }
  if (is.null(DP) || nrow(DP) == 0 || ncol(DP) == 0) {
    warning("No species with data after filtering; returning empty PSSD.")
    return(matrix(NA_real_, nrow = 0, ncol = SIM))
  }

  # Ensure trapezoid helpers are in scope for mc2d::rtrunc()
  rtrapezoid <- trapezoid::rtrapezoid
  ptrapezoid <- trapezoid::ptrapezoid
  qtrapezoid <- trapezoid::qtrapezoid
  dtrapezoid <- trapezoid::dtrapezoid

  # test if there is no data available for one species
  if (any(apply(DP, 2, function(x) length(which(!is.na(x)))) == 0)) {
    warning(
      "No data is available for one or more species, it/they won't contribute to the PSSD calculation."
    )
    # find which species has no data
    ind.sp.rem <- which(apply(DP, 2, function(x) length(which(!is.na(x)))) == 0)
    # remove those columns
    DP <- DP[, -ind.sp.rem, drop = FALSE]
    UFt <- UFt[, -ind.sp.rem, drop = FALSE]
    UFdd <- UFdd[, -ind.sp.rem, drop = FALSE]
  }

  if (is.null(dim(DP)) || ncol(DP) == 0) {
    warning("No species with data after filtering; returning empty PSSD.")
    return(matrix(NA_real_, nrow = 0, ncol = SIM))
  }

  # Create the step distributions (or triangular or trapezoidal) for each species
  # Create an empty matrix in which step distributions will be compiled
  NOEC_comb <- matrix(NA, ncol(DP), SIM, dimnames = list(colnames(DP), NULL))

  # Fill in the matrix. If there is only one data point, NOEC stays the same. If  there are
  # 2 endpoints, a uniform distribution is produced. If  there are more than 2 endpoints, a step
  # distribution is produced. One line is for one species.
  # store the corrected endpoints
  corr.endpoints <- DP / (UFdd * UFt)
  sort.endpoints <- apply(corr.endpoints, 2, sort)

  for (sp in colnames(DP)) {
    # store the indices of the minimal and maximal data point
    ind.min <- which.min(corr.endpoints[, sp])
    ind.max <- which.max(corr.endpoints[, sp])

    # calculate the theoretical minimum and maximum of the distribution we are looking for
    sp.min <- corr.endpoints[ind.min, sp] *
      (1 -
        (sqrt((CV.DP / 2.45)^2 + (CV.UF / 2.45)^2 + (CV.UF / 2.45)^2) * 2.45))
    sp.max <- corr.endpoints[ind.max, sp] *
      (1 +
        (sqrt((CV.DP / 2.45)^2 + (CV.UF / 2.45)^2 + (CV.UF / 2.45)^2) * 2.45))

    # For species with one unique data point, NOEC stays the same:
    if (length(unique(sort.endpoints[[sp]])) == 1) {
      NOEC_comb[sp, ] <- mc2d::rtrunc(
        "rtriang",
        min = sp.min,
        mode = sort.endpoints[[sp]][1],
        max = sp.max,
        n = SIM,
        linf = 0
      )

      # For species with two endpoints:
    } else if (length(sort.endpoints[[sp]]) == 2) {
      # Create a trapezoidal distribution including both endpoints
      NOEC_comb[sp, ] <- mc2d::rtrunc(
        "rtrapezoid",
        SIM,
        mode1 = sort.endpoints[[sp]][1],
        mode2 = sort.endpoints[[sp]][2],
        min = sp.min,
        max = sp.max,
        linf = 0
      )

      # For species with three endpoints or more:
    } else {
      # Sample from this step distribution for each species
      NOEC_comb[sp, ] <- rmore(
        values = sort.endpoints[[sp]],
        max = sp.max,
        min = sp.min,
        N = SIM,
        linf = 0
      )
    }
  }
  # return the whole matrix
  return(NOEC_comb)
}
