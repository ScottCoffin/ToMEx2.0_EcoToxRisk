#' Compute an HC distribution from per-species endpoints
#'
#' Thin, documented wrapper around the package's internal probabilistic-SSD
#' sampler, so externally generated, for example ML-predicted, per-species
#' effect concentrations can be run through the same engine as
#' [make_all_pSSDs()].
#'
#' @param endpoints A data frame with columns `Species` and `ec`. `Species`
#'   identifies species, and `ec` is a positive numeric effect concentration on
#'   the linear scale, in particles/L. Optional column `ec_sd` supplies
#'   endpoint-level standard deviations on the same scale for `DP.SD`.
#' @param sim Integer number of Monte Carlo iterations.
#' @param cv_dp Numeric coefficient of variation for endpoint uncertainty.
#' @param cv_uf Numeric coefficient of variation for uncertainty factors.
#' @param rmore_method Sampling method for species with three or more endpoints;
#'   either `"step"` or `"lognormal"`.
#' @param hcx Numeric hazard concentration quantile, for example `0.05` for HC5.
#' @param apply_assessment_factors Logical. Defaults to `FALSE`, meaning inputs
#'   are already NOEC-equivalent and assessment factors are treated as unit
#'   multipliers.
#'
#' @return A list with `pSSD`, a species x `sim` matrix of sampled endpoint
#'   concentrations, and `hc`, a numeric vector of HCx values, one per
#'   simulation.
#' @export
#'
#' @examples
#' endpoints <- data.frame(
#'   Species = rep(c("sp_a", "sp_b", "sp_c"), each = 3),
#'   ec = c(20, 25, 40, 80, 95, 120, 300, 330, 380)
#' )
#' res <- hc5_from_endpoints(endpoints, sim = 50)
#' str(res)
hc5_from_endpoints <- function(
  endpoints,
  sim = 1000,
  cv_dp = 0.15,
  cv_uf = 0.5,
  rmore_method = c("step", "lognormal"),
  hcx = 0.05,
  apply_assessment_factors = FALSE
) {
  rmore_method <- match.arg(rmore_method)

  if (!is.data.frame(endpoints)) {
    stop("`endpoints` must be a data frame.", call. = FALSE)
  }
  required_cols <- c("Species", "ec")
  missing_cols <- setdiff(required_cols, names(endpoints))
  if (length(missing_cols) > 0) {
    stop(
      "`endpoints` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.numeric(endpoints$ec) || any(!is.finite(endpoints$ec)) ||
      any(endpoints$ec <= 0)) {
    stop("`endpoints$ec` must contain positive finite numeric values.", call. = FALSE)
  }
  if ("ec_sd" %in% names(endpoints) &&
      (!is.numeric(endpoints$ec_sd) || any(!is.finite(endpoints$ec_sd)) ||
       any(endpoints$ec_sd < 0))) {
    stop("`endpoints$ec_sd` must contain finite non-negative numeric values.", call. = FALSE)
  }
  if (!is.numeric(sim) || length(sim) != 1L || !is.finite(sim) || sim < 1) {
    stop("`sim` must be a positive integer.", call. = FALSE)
  }
  sim <- as.integer(sim)
  if (!is.numeric(hcx) || length(hcx) != 1L || !is.finite(hcx) ||
      hcx <= 0 || hcx >= 1) {
    stop("`hcx` must be a single numeric value between 0 and 1.", call. = FALSE)
  }
  if (!is.logical(apply_assessment_factors) ||
      length(apply_assessment_factors) != 1L ||
      is.na(apply_assessment_factors)) {
    stop("`apply_assessment_factors` must be TRUE or FALSE.", call. = FALSE)
  }

  endpoints <- endpoints[order(endpoints$Species), , drop = FALSE]
  endpoints$.row <- stats::ave(
    seq_len(nrow(endpoints)),
    endpoints$Species,
    FUN = seq_along
  )

  DP <- reshape2::acast(endpoints, .row ~ Species, value.var = "ec", drop = FALSE)
  DP.SD <- if ("ec_sd" %in% names(endpoints)) {
    reshape2::acast(endpoints, .row ~ Species, value.var = "ec_sd", drop = FALSE)
  } else {
    NULL
  }

  pSSD <- do.pSSD_mod(
    DP = DP,
    DP.SD = DP.SD,
    UFt = NULL,
    UFdd = NULL,
    SIM = sim,
    CV.DP = cv_dp,
    CV.UF = cv_uf,
    rmore_method = rmore_method,
    apply_assessment_factors = apply_assessment_factors
  )

  list(
    pSSD = pSSD,
    hc = apply(pSSD, 2, stats::quantile, probs = hcx, type = 8)
  )
}
