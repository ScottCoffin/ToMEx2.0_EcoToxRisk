# Function for generating probability distributions based on more than 2 modal values and their 
# uncertainties
# 
# Yields a vector of length: number of iterations in the simulation
# 
# Arguments:   - values : vector of modal values
#              - coef : vector of uncertainty coefficients. Defaults to NULL. Either coef or max
#                and min need to be provided
#              - max : maximum of the probability distribution
#              - min : minimum of the probability distribution
#              - N : number of iterations in the simulation
#              - linf : lower truncation limit (useful for example to avoid negative concentrations).
#                Defaults to -Inf for no lower truncation
#              - lsup : upper truncation limit. Defaults to Inf for no upper truncation
# 
# Date of last modification: 03.07.2019
# 
# Associated publication: Systematic consideration of parameter uncertainty and variability in
#                         probabilistic species sensitivity distributions
# 
# Authors: Henning Wigger, Delphine Kawecki, Bernd Nowack and V?ronique Adam
# 
# Institute: Empa, Swiss Federal Laboratories for Materials Science and Technology,
#            Technology and Society Laboratory, Lerchenfeldstrasse 5, 9014 St. Gallen, Switzerland
# 
# submitted to Integrated Environmental Assessment and Management in July 2019

# -------------------------------------------------------------------------------------------------

require(crayon)

# Helper function for parameter validation
validate_rtrunc_args <- function(dist, minv, maxv, mode = NULL, lsup, linf) {
  # Check for NA, NaN, or non-numeric
  vals <- c(minv, maxv, lsup, linf)
  if (!is.numeric(vals) || any(is.na(vals) | is.nan(vals))) {
    stop("All rtrunc arguments must be numeric and not NA/NaN.")
  }
  if (lsup <= linf) stop("lsup must be greater than linf.")
  if (lsup < minv) stop("lsup must be >= min.")
  if (linf > maxv) stop("linf must be <= max.")
  if (dist == "rtriang") {
    if (is.null(mode)) stop("mode must be provided for rtriang.")
    if (!is.numeric(mode) || is.na(mode) || is.nan(mode)) stop("mode must be numeric and not NA/NaN.")
    if (minv > mode || mode > maxv) stop("min <= mode <= max must hold for rtriang.")
    if (minv >= maxv) stop("min < max must hold for rtriang.")
  } else if (dist == "runif") {
    if (minv >= maxv) stop("min < max must hold for runif.")
  }
}

rmore <- function(values,
                  coef = NULL,
                  max = NULL,
                  min = NULL,
                  N,
                  linf = -Inf,
                  lsup = Inf) {
  # Input checks
  stopifnot(any(!is.null(max), !is.null(min), !is.null(coef)))
  if (!is.null(max) | !is.null(min)) {
    stopifnot(!is.null(max),
              !is.null(min),
              !any(values > max),
              !any(values < min))
  }
  require(trapezoid)
  require(mc2d)
  
  # Calculate min/max
  if (!is.null(coef)) {
    dist.min <- values[which.min(values)] * (1 - coef[which.min(values)])
    dist.max <- values[which.max(values)] * (1 + coef[which.max(values)])
  } else {
    dist.min <- min
    dist.max <- max
  }
  
  sort.values <- sort(values)
  freq.uni.values <- table(sort.values)
  uni.values <- as.numeric(names(freq.uni.values))
  val.min <- sort.values[1]
  val.max <- sort.values[length(sort.values)]
  
  # Heights
  height.values <- rep(NA, length(uni.values))
  for (i in 2:(length(uni.values)-1)) {
    height.values[i] <- max(freq.uni.values[i]*N/(uni.values[i+1]-uni.values[i]),
                            freq.uni.values[i]*N/(uni.values[i]-uni.values[i-1]))
  }
  height.values[1] <- freq.uni.values[1]*N/(uni.values[2]-uni.values[1])
  height.values[length(height.values)] <-
    freq.uni.values[length(height.values)]*N/(uni.values[length(height.values)]-
                                                uni.values[length(height.values)-1])
  
  # LEFT TRIANGLE
  if ((length(unique(sort.values)) < length(sort.values)) &
      (val.min == sort.values[2])) {
    n <- length(which(sort.values == val.min))
    val.minP1 <- sort.values[n+1]
    h <- height.values[1]
    A <- (val.min - dist.min)*h/2
    if (dist.min < linf) {
      left <- NULL
    } else {
      # Validation
      validate_rtrunc_args("rtriang", dist.min, val.min, val.min, lsup, linf)
      left <- rtrunc("rtriang",
                     n = N,
                     min = dist.min,
                     mode = val.min,
                     max = val.min,
                     lsup = lsup,
                     linf = linf)
    }
  } else {
    A <- N*(val.min-dist.min)/(2*(sort.values[2]-val.min))
    validate_rtrunc_args("rtriang", dist.min, val.min, val.min, lsup, linf)
    left <- rtrunc("rtriang",
                   n = N,
                   min = dist.min,
                   mode = val.min,
                   max = val.min,
                   lsup = lsup,
                   linf = linf)
  }
  
  # MID SEGMENTS
  mid <- list()
  for (i in 1:(length(uni.values)-1)) {
    if (height.values[i] == 1 & height.values[i+1] == 1) {
      validate_rtrunc_args("runif", uni.values[i], uni.values[i+1], lsup = lsup, linf = linf)
      mid[[i]] <- rtrunc("runif",
                         n = N,
                         min = uni.values[i],
                         max = uni.values[i+1],
                         lsup = lsup,
                         linf = linf)
    } else if (isTRUE(all.equal(height.values[i], height.values[i+1]))) {
      validate_rtrunc_args("runif", uni.values[i], uni.values[i+1], lsup = lsup, linf = linf)
      mid[[i]] <- rtrunc("runif",
                         n = N*freq.uni.values[i],
                         min = uni.values[i],
                         max = uni.values[i+1],
                         lsup = lsup,
                         linf = linf)
    } else if (height.values[i] > height.values[i+1]) {
      h1 <- height.values[i]
      h2 <- height.values[i+1]
      max.trunc.distr <- (h1*uni.values[i+1]-h2*uni.values[i])/(h1-h2)
      A <- ((max.trunc.distr - uni.values[i])*h1 - (max.trunc.distr - uni.values[i+1])*h2)/2
      validate_rtrunc_args("rtriang", uni.values[i], max.trunc.distr, uni.values[i], lsup = min(uni.values[i+1], lsup), linf = linf)
      mid[[i]] <- rtrunc("rtriang",
                         n = N,
                         min = uni.values[i],
                         mode = uni.values[i],
                         max = max.trunc.distr,
                         lsup = min(uni.values[i+1], lsup),
                         linf = linf)
    } else if (height.values[i] < height.values[i+1]) {
      h1 <- height.values[i]
      h2 <- height.values[i+1]
      min.trunc.distr <- (h2*uni.values[i] - h1*uni.values[i+1])/(h2 - h1)
      A <- ((uni.values[i+1] - min.trunc.distr)*h2 - (uni.values[i] - min.trunc.distr)*h1)/2
      validate_rtrunc_args("rtriang", min.trunc.distr, uni.values[i+1], uni.values[i+1], lsup = lsup, linf = max(uni.values[i], linf))
      mid[[i]] <- rtrunc("rtriang",
                         n = N,
                         min = min.trunc.distr,
                         mode = uni.values[i+1],
                         max = uni.values[i+1],
                         lsup = lsup,
                         linf = max(uni.values[i], linf))
    }
  }
  
  # RIGHT TRIANGLE
  if ((length(unique(sort.values)) < length(sort.values)) &
      (val.max == sort.values[length(sort.values)-1])) {
    n <- length(which(sort.values == val.max))
    val.maxM1 <- sort.values[length(sort.values)-n]
    h <- height.values[length(height.values)]
    A <- (dist.max - val.max)*h /2
    if (dist.max > lsup) {
      right <- NULL
    } else {
      validate_rtrunc_args("rtriang", val.max, dist.max, val.max, lsup, linf)
      right <- rtrunc("rtriang",
                      n = N,
                      min = val.max,
                      mode = val.max,
                      max = dist.max,
                      lsup = lsup,
                      linf = linf)
    }
  } else {
    A <- N*(dist.max-val.max)/(2*(val.max-sort.values[length(sort.values)-1]))
    validate_rtrunc_args("rtriang", sort.values[length(sort.values)], dist.max, sort.values[length(sort.values)], lsup, linf)
    right <- rtrunc("rtriang",
                    n = N, 
                    min = sort.values[length(sort.values)], 
                    mode = sort.values[length(sort.values)],
                    max = dist.max,
                    lsup = lsup,
                    linf = linf)
  }
  
  # Combine
  step_distr <- c(left, do.call("c", mid), right)
  return(sample(step_distr, N))
}

