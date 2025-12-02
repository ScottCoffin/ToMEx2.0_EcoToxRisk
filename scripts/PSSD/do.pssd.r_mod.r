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
#### Scott's attempt at fixing the negative number issues #####
#################################################################
do.pSSD_mod <- function(DP, DP.SD, UFt, UFdd, SIM, CV.DP, CV.UF) {
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
    
    # Debugging: Print sp.min and sp.max to ensure they are non-negative
    # print(paste("Species:", sp, "sp.min:", sp.min, "sp.max:", sp.max))
    
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
      
    } else {
      # Three or more endpoints: Bootstrap with uncertainty factor
      #Removing the DP.SD[, sp] from the groups to test
      high <- (sqrt(sum(c(mean(DP.SD[, sp]/DP[, sp], na.rm = T), 1+CV.DP2, 1+CV.UF)^2, na.rm = T)))
      low <- 1/high
      
      # Debugging: Print low and high to ensure they are valid
      #  print(paste("Species:", sp, "low:", low, "high:", high))
      
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


###################################
############### PSSD++ Method #######################
#######################################################
#Simulates species-specific effect concentrations using triangular (n =1) /trapezoidal (n = 2) or bootstrap distributions (n >2),
# improving on PSSD+ by supporting lognormal or empirical bootstrapping and handling uncertainty more efficiently.

do.pSSD_mod_multiMode <- function(DP, DP.SD, UFt, UFdd, SIM, CV.DP, CV.UF, method = "empirical") {
  
  
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
    
    # Calculate sp.min and sp.max with non-negativity enforced
    sp.min <- max(0, corr.endpoints[ind.min, sp] * (1 - (sqrt(sum((CV.DP2 / 2.45)^2) + 2 * (CV.UF / 2.45)^2) * 2.45)))
    sp.max <- max(0, corr.endpoints[ind.max, sp] * (1 + (sqrt(sum((CV.DP2 / 2.45)^2) + 2 * (CV.UF / 2.45)^2) * 2.45)))
    
    if (length(unique(sort.endpoints[[sp]])) == 1) {
      # Single endpoint: Use truncated triangular distribution
      NOEC_comb[sp, ] <- rtrunc("rtriang", min = sp.min,
                                mode = sort.endpoints[[sp]][1],
                                max = sp.max, n = SIM, linf = 0)
      
    } else if (length(sort.endpoints[[sp]]) == 2) {
      # Two endpoints: Use truncated trapezoidal distribution
      NOEC_comb[sp, ] <- rtrunc("rtrapezoid", SIM,
                                mode1 = sort.endpoints[[sp]][1],
                                mode2 = sort.endpoints[[sp]][2],
                                min = sp.min, max = sp.max, linf = 0)
      
    } else {
      # Three or more endpoints: lognormal or empirical bootstrapping
      observed <- sort.endpoints[[sp]]
      
      # Uncertainty factor (lognormal-like range)
      high <- sqrt(sum(c(mean(DP.SD[, sp] / DP[, sp], na.rm = TRUE), 1 + CV.DP2, 1 + CV.UF)^2, na.rm = TRUE))
      low <- 1 / high
      uncertainty_factor <- runif(SIM, min = low, max = high)
      
      # use lognormal dist
      if (tolower(method) == "lognormal") {
        data <- 10^rnorm(mean = mean(log10(observed)),
                         sd = sd(log10(observed)),
                         n = SIM)
        #NOEC_comb[sp, ] <- pmax(data * uncertainty_factor, 1e-6)
        NOEC_comb[sp, ] <- data * uncertainty_factor
        
        # Use empirical bootstrapping
      } else if (tolower(method) == "empirical") {
        boot_sample <- sample(observed, size = SIM, replace = TRUE)
        jitter_width <- max(0.05 * IQR(observed), 1e-6)
        jittered <- boot_sample + runif(SIM, min = 0, max = jitter_width)
        NOEC_comb[sp, ] <- pmax(jittered * uncertainty_factor, 1e-6)
      
        # use CFD method from PSSD+
      } else if (tolower(method) == "rmore") {
        coef <- rep(sqrt(sum(c(mean(DP.SD[, sp] / DP[, sp], na.rm = TRUE), 1 + CV.DP2, 1 + CV.UF)^2, na.rm = TRUE)), length(observed))
        
        NOEC_comb[sp, ] <- pmax(rmore(values = observed, coef = coef, N = SIM, linf = 1e-6), 1e-6)
        
        # Sample from this step distribution for each species
        #NOEC_comb[sp,] <- rmore(values = sort.endpoints[[sp]], max = sp.max, min = sp.min, N = SIM, linf = 0)
      
      } else {
        stop("Invalid method. Use 'lognormal' or 'empirical'.")
      }
    }
  
  return(NOEC_comb)
  }
}


########################## Original Function from https://zenodo.org/records/3267194 #########
do.pSSD <- function(DP,
                    UFt,
                    UFdd,
                    SIM,
                    CV.DP,
                    CV.UF){
  
  # test if there is no data available for one species
  if(any(apply(DP,2,function(x) length(which(!is.na(x)))) == 0)){
    warning("No data is available for one or more species, it/they won't contribute to the PSSD calculation.")
    # find which species has no data
    ind.sp.rem <- which(apply(DP,2,function(x) length(which(!is.na(x)))) == 0)
    # remove those columns
    DP <- DP[,-ind.sp.rem]
    UFt <- UFt[,-ind.sp.rem]
    UFdd <- UFdd[,-ind.sp.rem]
  }
  
  # Create the step distributions (or triangular or trapezoidal) for each species
  # Create an empty matrix in which step distributions will be compiled
  NOEC_comb <- matrix(NA, ncol(DP), SIM,
                      dimnames = list(colnames(DP), NULL))
  
  # Fill in the matrix. If there is only one data point, NOEC stays the same. If  there are
  # 2 endpoints, a uniform distribution is produced. If  there are more than 2 endpoints, a step
  # distribution is produced. One line is for one species.
  require(trapezoid)
  require(mc2d)
  
  # store the corrected endpoints
  corr.endpoints <- DP/(UFdd*UFt)
  sort.endpoints <- apply(corr.endpoints, 2, sort)
  
  for (sp in colnames(DP)){
    
    # store the indices of the minimal and maximal data point
    ind.min <- which.min(corr.endpoints[,sp])
    ind.max <- which.max(corr.endpoints[,sp])
    
    # calculate the theoretical minimum and maximum of the distribution we are looking for
    sp.min <- corr.endpoints[ind.min,sp]*(1-(sqrt((CV.DP/2.45)^2 + (CV.UF/2.45)^2 + (CV.UF/2.45)^2)*2.45))
    sp.max <- corr.endpoints[ind.max,sp]*(1+(sqrt((CV.DP/2.45)^2 + (CV.UF/2.45)^2 + (CV.UF/2.45)^2)*2.45))
    
    # For species with one unique data point, NOEC stays the same:
    if(length(unique(sort.endpoints[[sp]])) == 1){
      NOEC_comb[sp,] <- rtrunc("rtriang", min = sp.min, 
                               mode = sort.endpoints[[sp]][1],
                               max = sp.max,
                               n = SIM, linf = 0)
      
      
      # For species with two endpoints:
    } else if(length(sort.endpoints[[sp]]) == 2){
      # Create a trapezoidal distribution including both endpoints
      NOEC_comb[sp,] <- rtrunc("rtrapezoid", SIM,
                               mode1 = sort.endpoints[[sp]][1],
                               mode2 = sort.endpoints[[sp]][2],
                               min = sp.min, max = sp.max,
                               linf = 0)
      
      
      # For species with three endpoints or more:
    } else {
      
      
      # Sample from this step distribution for each species
      NOEC_comb[sp,] <- rmore(values = sort.endpoints[[sp]], max = sp.max, min = sp.min, N = SIM, linf = 0)
      
    } 
  }
  # return the whole matrix
  return(NOEC_comb)
  
}
