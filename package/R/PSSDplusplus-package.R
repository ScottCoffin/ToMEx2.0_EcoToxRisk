#' PSSDplusplus: Probabilistic Microplastic Ecotoxicity Workflow (PSSD++)
#'
#' Implements the PSSD++ probabilistic ecotoxicological risk framework for
#' microplastics, as described in Coffin et al. (2026,
#' \doi{10.1016/j.jhazmat.2025.141021}). The package provides Monte Carlo
#' alignment of laboratory toxicity data to environmentally relevant metrics
#' (ERMs) -- food dilution and tissue translocation -- via Latin Hypercube
#' Sampling of particle-trait parameters. Aligned data are used to fit
#' probabilistic Species Sensitivity Distributions (pSSDs) and derive Predicted
#' No-Effect Concentrations (PNECs) across surface-water and sediment compartments
#' for freshwater and marine environments. Bundles the preprocessed ToMEx 2.0
#' microplastic ecotoxicity database (~13,000 data points) and a compact 10-species
#' subset for rapid testing. Supports parallel computation, deterministic caching,
#' and structured diagnostic output.
#'
#' @references Coffin, S. et al. (2026). \emph{Journal of Hazardous Materials}
#' 503, 141021. \doi{10.1016/j.jhazmat.2025.141021}
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
