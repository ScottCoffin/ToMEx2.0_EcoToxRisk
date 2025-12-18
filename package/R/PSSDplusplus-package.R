#' PSSDplusplus: Probabilistic Microplastic Ecotoxicity Workflow
#'
#' Implements the Monte Carlo alignment and PSSD++ analysis used for the ToMEx 2.0
#' microplastic ecotoxicity assessment. The package bundles helper functions,
#' plotting utilities, and workflow wrappers to reproduce the manuscript
#' thresholds with the bundled dataset.
#'
#' @keywords internal
#' @import dplyr tidyr purrr tibble stringr ggplot2 sensobol truncnorm doParallel doSNOW tictoc reshape2 future future.apply progressr ssdtools parallel utils mc2d trapezoid
#' @importFrom magrittr %>%
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix Matrix
#' @importFrom crayon blue green red yellow
#' @importFrom data.table as.data.table
#' @importFrom foreach %dopar% foreach
#' @importFrom scales math_format
#' @importFrom stats quantile rnorm predict
"_PACKAGE"
