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

# -------------------------------------------------------------------------------------------------


#################################################################
#### PSSD++ Function #####
#################################################################
do.pSSD_mod <- function(DP, DP.SD, UFt, UFdd, SIM, CV.DP, CV.UF, rmore_method = "step") {
  # Check for species with no data
  if (any(apply(DP, 2, function(x) length(which(!is.na(x)))) == 0)) {
    warning("No data is available for one or more species, it/they won't contribute to the PSSD calculation.")
    ind.sp.rem <- which(apply(DP, 2, function(x) length(which(!is.na(x)))) == 0)
    DP <- DP[, -ind.sp.rem]
    DP.SD <- DP.SD[, -ind.sp.rem]
    UFt <- UFt[, -ind.sp.rem]
    UFdd <- UFdd[, -ind.sp.rem]
  }
  
  require(trapezoid)
  require(mc2d)
  
  
  
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
    sp.min <- max(0, corr.endpoints[ind.min, sp] * (1 - (sqrt(sum((CV.DP2 / 2.45)^2) + 2 * (CV.UF / 2.45)^2) * 2.45)))
    sp.max <- max(0, corr.endpoints[ind.max, sp] * (1 + (sqrt(sum((CV.DP2 / 2.45)^2) + 2 * (CV.UF / 2.45)^2) * 2.45)))
    
    # Handle cases based on the number of unique endpoints
    if (length(unique(sort.endpoints[[sp]])) == 1) {
      # Single endpoint: Use truncated triangular distribution
      NOEC_comb[sp, ] <- rtrunc("rtriang", min = sp.min, 
                                mode = sort.endpoints[[sp]][1], 
                                max = sp.max, 
                                n = SIM, linf = 0)
      
    } else if (length(sort.endpoints[[sp]]) == 2) {
      # Two endpoints: Use truncated trapezoidal distribution
      NOEC_comb[sp, ] <- rtrunc("rtrapezoid", SIM, 
                                mode1 = sort.endpoints[[sp]][1], 
                                mode2 = sort.endpoints[[sp]][2], 
                                min = sp.min, max = sp.max, linf = 0)
      
      ### Three endpoints: Use rmore function based on the specified method  ##
      # (step-wise trapezoidal - i.e., original Wigger et al. (2020))
    } else if (length(sort.endpoints[[sp]]) > 2 & rmore_method == "step") {
      # Three or more endpoints: Use step distribution
      NOEC_comb[sp,] <- rmore(values = sort.endpoints[[sp]], max = sp.max, min = sp.min, N = SIM, linf = 0)
      
      # Log10 -> normal distribution bootstrapping shortcut method developed here
    } else if (length(sort.endpoints[[sp]]) > 2 & rmore_method == "lognormal") {
    
    # Three or more endpoints: Log-Normal distribution bootstrapping with uncertainty factor
    #Removing the DP.SD[, sp] from the groups to test
    high <- (sqrt(sum(c(mean(DP.SD[, sp]/DP[, sp], na.rm = T), 1+CV.DP2, 1+CV.UF)^2, na.rm = T)))
    low <- 1/high
    
    uncertainty_factor <- runif(min = low, max = high, n = SIM)
    
    data <- 10^rnorm(mean = mean(log10(sort.endpoints[[sp]])), 
                     sd = sd(log10(sort.endpoints[[sp]])),  
                     n = SIM)
    
    # Ensure no negative values in the final result
    NOEC_comb[sp, ] <- data * uncertainty_factor
    
  }
  }
  
  # Return the final matrix
  return(NOEC_comb)
}
