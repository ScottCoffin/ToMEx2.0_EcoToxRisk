#' ToMEx 2.0 ecotoxicity dataset
#'
#' Microplastics ecotoxicity records used throughout the PSSD++ workflow.
#' Data are bundled as a plain data frame identical to the source RDS shipped
#' in `inst/extdata/aoc_z_tomex2.RDS`.
#'
#' @format A data frame with the fields described in the ToMEx 2.0 metadata.
#' @usage data(tomex2)
"tomex2"

#' Default alignment parameters (means and SDs)
#'
#' Vectorized defaults for the Monte Carlo alignment sampler
#' (\link[=matrix_function]{matrix_function()}), providing means and standard deviations for length/width
#' ratios, densities, alpha parameters, height:width ratios, translocation size
#' limit, and body-length model coefficients.
#'
#' @format A data frame of parameter means and standard deviations.
#' @usage data(param_default_values)
#' @export
"param_default_values"
