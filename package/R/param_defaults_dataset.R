#' Default alignment parameters
#'
#' Backing data frame for `PSSDplusplus::param_default_values`, matching the
#' defaults used in the alignment sampler.
#'
#' @name param_default_values
#' @docType data
#' @format A data frame with one row of default means/standard deviations for alignment parameters.
#' @keywords internal
param_default_values <- data.frame(
  # width to length ratios
  R.ave.water.marine = 0.77,
  R.ave.water.marine.sd = 0.29,
  R.ave.water.freshwater = 0.67,
  R.ave.water.freshwater.sd = 0.28,
  R.ave.sediment.marine = 0.75,
  R.ave.sediment.marine.sd = 0.30,
  R.ave.sediment.freshwater = 0.70,
  R.ave.sediment.freshwater.sd = 0.33,

  # density (g/cm^3)
  p.ave.marine = 1.10,
  p.ave.marine.sd = 0.14,
  p.ave.freshwater = 1.04,
  p.ave.freshwater.sd = 0.12,

  # alpha values
  alpha.marine = 2.07,
  alpha.marine.sd = 0.03, # length
  a.sa.marine = 1.50,
  a.sa.marine.sd = 0.009,
  a.v.marine = 1.48,
  a.v.marine.sd = 0.06,
  a.m.marine = 1.45,
  a.m.marine.sd = 0.08,
  alpha.freshwater = 2.25,
  alpha.freshwater.sd = 0.19,
  a.sa.freshwater = 1.55,
  a.sa.freshwater.sd = 0.028,
  a.v.freshwater = 1.68,
  a.v.freshwater.sd = 0.081,
  a.m.freshwater = 1.65,
  a.m.freshwater.sd = 0.071,
  a.ssa.freshwater = 2.71,
  a.ssa.freshwater.sd = 0.009,
  beta_log10_body_length = 0.9341,
  se_beta_log10_body_length = 0.1376,
  body_length_intercept = 1.1200,
  se_body_length_intercept = 0.3222,
  beta_0 = 1.308344,
  se_beta_0 = 0.3963612,
  beta_1 = -0.01468148,
  se_beta_1 = 0.006657993
)
