##### FUNCTIONS #####

### Generalizable function that works on any value (alpha == 1 and == 2 are limits!)
mux_polyfnx <- function(a.x, x_UL, x_LL) {
  # Validate inputs
  if (length(a.x) != length(x_UL) || length(a.x) != length(x_LL)) {
    stop("a.x, x_UL, and x_LL must have the same length.")
  }

  # Initialize result vector
  mux.poly <- numeric(length(a.x))

  # Loop through each element to handle row-by-row logic
  for (i in seq_along(a.x)) {
    if (is.na(a.x[i]) || is.na(x_UL[i]) || is.na(x_LL[i])) {
      # Handle NA values
      mux.poly[i] <- NA
    } else if (a.x[i] == 1) {
      # Special case: a.x == 1
      if (x_UL[i] > 0 && x_LL[i] > 0) {
        mux.poly[i] <- (x_UL[i] - x_LL[i]) / log(x_UL[i] / x_LL[i])
      } else {
        mux.poly[i] <- NA # Invalid input for log
      }
    } else if (a.x[i] == 2) {
      # Special case: a.x == 2
      epsilon <- 1e-10 # Small value to avoid division by zero
      if (x_UL[i] > 0 && x_LL[i] > 0) {
        mux.poly[i] <- log(x_UL[i] / x_LL[i]) /
          ((x_LL[i] + epsilon)^-1 - (x_UL[i] + epsilon)^-1)
      } else {
        mux.poly[i] <- NA # Invalid input for log
      }
    } else {
      # General case: a.x != 1 and a.x != 2
      if (x_UL[i] > 0 && x_LL[i] > 0) {
        mux.poly[i] <- ((1 - a.x[i]) / (2 - a.x[i])) *
          ((x_UL[i]^(2 - a.x[i]) - x_LL[i]^(2 - a.x[i])) /
            (x_UL[i]^(1 - a.x[i]) - x_LL[i]^(1 - a.x[i])))
      } else {
        mux.poly[i] <- NA # Invalid input for power calculations
      }
    }
  }

  # Return the result
  return(mux.poly)
}


############# VOLUME ############
volumefnx <- function(
  R = NA, # average length-to-width ratio for environment
  H_W_ratio = 0.67, # assumed 0.67 * width per Kooi et al. (2021)
  length, # particle length (always known)
  height = NA, # particle height (if known)
  width = NA # particle width (if known)
) {
  # If width unknown, use L:W ratio
  width <- ifelse(is.na(width), R * length, width)

  # If height unknown, use H:R ratio
  height <- ifelse(is.na(height), H_W_ratio * width, height)

  # Calculate volume
  volume <- (4 / 3) * pi * (length / 2) * (width / 2) * (height / 2)

  return(volume)
}

############### SURFACE AREA ##################
#surface area equation for elongated spheres
SAfnx = function(
  length,
  width = NA,
  height = NA,
  R = NA,
  H_W_ratio = 0.67 # assumed 0.67 * width per Kooi et al. (2021)
) {
  # If width unknown, use L:W ratio
  width <- ifelse(is.na(width), R * length, width)

  # If height unknown, use H:R ratio
  height <- ifelse(is.na(height), H_W_ratio * width, height)

  # a, b, and c are equivalent to 1/2th of the length, width, and height, respectively
  a <- 0.5 * length
  b <- 0.5 * width
  c <- 0.5 * height

  SA = (4 * pi) * ((((a * b)^1.6 + (a * c)^1.6 + (b * c)^1.6) / 3)^(1 / 1.6))
  return(SA)
}

################# MASS ####################
massfnx <- function(v, p) {
  # If either v or p is NA, return NA for those elements
  mass <- ifelse(is.na(v) | is.na(p), NA, p * v * (1 / 1e12) * 1e6) # correction factor (g to µg)
  return(mass)
}

###### SSA #####
SSA.inversefnx = function(
  sa, # average surface area
  m
) {
  #average mass
  SSA.inverse = m / sa
  return(SSA.inverse)
}


#### Ecologically Relevant Metric Functions (used in reactives with user-input params) ####

###function to derive correction factor (CF) from Koelmans et al (equation 2)
CFfnx = function(
  a, #default alpha from Koelmans et al (2020)
  x2D, #set detault values to convert ranges to (1-5,000 um) #5mm is upper defuault
  x1D, #1 um is lower default size
  x2M,
  x1M
) {
  CF = (x2D^(1 - a) - x1D^(1 - a)) / (x2M^(1 - a) - x1M^(1 - a))
  return(CF)
}

#### equations for mu_x_poly (note that there are three depending on certain alphas for limits of equation)
### Generalizable function that works on any value (alpha == 1 and == 2 are limits!)
mux_polyfnx_generalizable <- function(a.x, x_UL, x_LL) {
  # Validate inputs
  if (length(a.x) != length(x_UL) || length(a.x) != length(x_LL)) {
    stop("a.x, x_UL, and x_LL must have the same length.")
  }

  # Initialize result vector
  mux.poly <- numeric(length(a.x))

  # Loop through each element to handle row-by-row logic
  for (i in seq_along(a.x)) {
    if (is.na(a.x[i]) || is.na(x_UL[i]) || is.na(x_LL[i])) {
      # Handle NA values
      mux.poly[i] <- NA
    } else if (a.x[i] == 1) {
      # Special case: a.x == 1
      if (x_UL[i] > 0 && x_LL[i] > 0) {
        mux.poly[i] <- (x_UL[i] - x_LL[i]) / log(x_UL[i] / x_LL[i])
      } else {
        mux.poly[i] <- NA # Invalid input for log
      }
    } else if (a.x[i] == 2) {
      # Special case: a.x == 2
      epsilon <- 1e-10 # Small value to avoid division by zero
      if (x_UL[i] > 0 && x_LL[i] > 0) {
        mux.poly[i] <- log(x_UL[i] / x_LL[i]) /
          ((x_LL[i] + epsilon)^-1 - (x_UL[i] + epsilon)^-1)
      } else {
        mux.poly[i] <- NA # Invalid input for log
      }
    } else {
      # General case: a.x != 1 and a.x != 2
      if (x_UL[i] > 0 && x_LL[i] > 0) {
        mux.poly[i] <- ((1 - a.x[i]) / (2 - a.x[i])) *
          ((x_UL[i]^(2 - a.x[i]) - x_LL[i]^(2 - a.x[i])) /
            (x_UL[i]^(1 - a.x[i]) - x_LL[i]^(1 - a.x[i])))
      } else {
        mux.poly[i] <- NA # Invalid input for power calculations
      }
    }
  }

  # Return the result
  return(mux.poly)
}
#
# #max ingestible specific surface area
SSA.inversefnx = function(
  sa, #surface area, calcaulted elsewhere
  m
) {
  #mass, calculated elsewhere
  SSA.inverse = m / sa
  return(SSA.inverse)
}

#data tidying functions from Ana

############## Levels summary ##################
summarize_and_print <- function(data, column_name) {
  result <- data %>%
    dplyr::group_by({{ column_name }}) %>%
    dplyr::summarise(n_datapoints = n()) %>%
    dplyr::arrange(as.numeric(as.character({{ column_name }}))) %>%
    print(n = 1000)
  return(result)
}


############## Change particle length ##################
update_particle_length <- function(
  data,
  doi,
  length,
  polymer,
  shape,
  new_value
) {
  data$Particle.Length..μm.[
    data$DOI == doi &
      data$Particle.Length..μm. == length &
      data$Polymer == polymer &
      data$Shape == shape
  ] = new_value

  return(data)
}


###### check what is missing #########
generate_structure_checks <- function(data) {
  structure.checks <- data.frame(
    na.counts = sapply(data, function(x) sum(is.na(x))),
    na.percent = round(
      sapply(data, function(x) sum(is.na(x)) / nrow(data) * 100),
      digits = 1
    ),
    n.levels = sapply(data, function(x) length(unique(x)))
  )
  return(structure.checks)
}


#### Define function for presenting thresholds in alternative units (volume, mass, etc.)

#Thresholds are presented in multiple exposure metrics, aligned to ERMs of interest (surface area for translocatable particles, volume for ingestible particles). Data are aligned using the methods described in Kooi et al 2021.
#The environmentally-relevant, ERM-aligned threshold in particles/L is converted to surface area (um2/L), volume (um3/L) and mass (ug/L) using equations 5 and 4 from Kooi et al 2021.
#$EC_{env,ERM, x} = \mu_{x,poly} * EC_{env}$
#Where $\mu_{x,poly}$ is derived using the following equation: $\mu_{x,poly} = \frac{1 - a_{x}}{2 - a_{x}}  \frac{X^{2-a_x}_{UL} - X^{2-a_x}_{LL}}{X^{1-a_x}_{UL} - X^{1-a_x}_{LL}}$
#In this case, the limits LL and UL in $\mu_{x,poly}$ relate to the values of the ERM at the 1 and 5,000 um size limits, respectively (i.e., SA, V, and M), rather than to bioavailability limits. Compartment-specific alpha values are used based on Table S4 of Kooi et al 2021.

convert_units_fxn <- function(thresholds_df, environment, params) {
  if (environment == "Marine") {
    a.sa.input = params$a.sa.marine
    a.m.input = params$a.m.marine
    a.v.input = params$a.v.marine
  }
  if (environment == "Freshwater") {
    a.sa.input = params$a.sa.freshwater
    a.m.input = params$a.m.freshwater
    a.v.input = params$a.v.freshwater
  }

  # report the units they're aligned to (um3/L and um2/L)
  thresholds_df %>%
    # convert tissue translocation thresholds
    dplyr::mutate(
      tissue_translocation_um2.L = `Tissue Translocation (Default)` *
        mux_polyfnx(
          a.x = a.sa.input, #alpha for marine surface water for surface area
          x_UL = x2D_set, #upper limit of default range
          x_LL = x1D_set # lower limit of default range
        )
    ) %>%
    dplyr::mutate(
      tissue_translocation_ug.L = `Tissue Translocation (Default)` *
        mux_polyfnx(
          a.x = a.m.input, #alpha for marine surface water for mass
          x_UL = x2D_set, #upper limit of default range
          x_LL = x1D_set # lower limit of default range
        )
    ) %>%
    ## convert food dilution thresholds
    dplyr::mutate(
      food_dilution_um3.L = `Food Dilution (Default)` *
        mux_polyfnx(
          a.x = a.v.input, #alpha for marine surface water for surface area
          x_UL = x2D_set, #upper limit of default range
          x_LL = x1D_set # lower limit of default range
        )
    ) %>%
    dplyr::mutate(
      food_dilution_ug.L = `Food Dilution (Default)` *
        mux_polyfnx(
          a.x = a.m.input, #alpha for marine surface water for mass
          x_UL = x2D_set, #upper limit of default range
          x_LL = x1D_set # lower limit of default range
        )
    )
}

# example
# convert_units_fxn(thresholds_df = marine_thresholds,
#                   environment = "Marine",
#                   params = params)

#### Function to generate synthetic data for MC alignments
generate_data <- function(
  n_sim = 10,
  ## width to length ratios
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
  alpha.marine.sd = 0.03, #length
  a.sa.marine = 1.50,
  a.sa.marine.sd = 0.009,
  a.v.marine = 1.48,
  a.v.marine.sd = 0.063,
  a.m.marine = 1.32,
  a.m.marine.sd = 0.009,
  a.ssa.marine = 1.98,
  a.ssa.marine.sd = 0.297,

  alpha.freshwater = 2.64,
  alpha.freshwater.sd = 0.01,
  a.sa.freshwater = 2.00,
  a.sa.freshwater.sd = 0.065,
  a.v.freshwater = 1.68,
  a.v.freshwater.sd = 0.081,
  a.m.freshwater = 1.65,
  a.m.freshwater.sd = 0.071,
  a.ssa.freshwater = 2.71,
  a.ssa.freshwater.sd = 0.009,

  # body length estimates
  beta_log10_body_length = 0.9341,
  se_beta_log10_body_length = 0.1376,
  body_length_intercept = 1.1200,
  se_body_length_intercept = 0.3222,
  beta_0 = 1.308344,
  se_beta_0 = 0.3963612,
  beta_1 = -0.01468148,
  se_beta_1 = 0.006657993
) {
  # Generate random samples based on user input or defaults

  ## width to length ratios. If R > 1, then width now becomes length - so it is meaningless. R cannot be negative.
  R.ave.water.marine.samples <- truncnorm::rtruncnorm(
    n_sim,
    mean = R.ave.water.marine,
    a = 0.0001,
    b = 1, #min and max
    sd = R.ave.water.marine.sd
  )

  H_W_ratio.marine_samples <- R.ave.water.marine.samples

  R.ave.water.freshwater.samples <- truncnorm::rtruncnorm(
    n_sim,
    mean = R.ave.water.freshwater,
    a = 0.0001,
    b = 1,
    sd = R.ave.water.freshwater.sd
  )

  H_W_ratio.freshwater_samples <- R.ave.water.freshwater.samples

  R.ave.sediment.marine.samples <- truncnorm::rtruncnorm(
    n_sim,
    mean = R.ave.sediment.marine,
    a = 0.0001,
    b = 1,
    sd = R.ave.sediment.marine.sd
  )
  R.ave.sediment.freshwater.samples <- truncnorm::rtruncnorm(
    n_sim,
    mean = R.ave.sediment.freshwater,
    a = 0.0001,
    b = 1,
    sd = R.ave.sediment.freshwater.sd
  )

  p.ave.marine.samples <- stats::rnorm(
    n_sim,
    mean = p.ave.marine,
    sd = p.ave.marine.sd
  )
  p.ave.freshwater.samples <- stats::rnorm(
    n_sim,
    mean = p.ave.freshwater,
    sd = p.ave.freshwater.sd
  )

  # Marine alignment properties
  alpha.marine.samples <- stats::rnorm(
    n_sim,
    mean = alpha.marine,
    sd = alpha.marine.sd
  )
  a.sa.marine.samples <- stats::rnorm(
    n_sim,
    mean = a.sa.marine,
    sd = a.sa.marine.sd
  )
  a.v.marine.samples <- stats::rnorm(
    n_sim,
    mean = a.v.marine,
    sd = a.v.marine.sd
  )
  a.m.marine.samples <- stats::rnorm(
    n_sim,
    mean = a.m.marine,
    sd = a.m.marine.sd
  )
  a.ssa.marine.samples <- stats::rnorm(
    n_sim,
    mean = a.ssa.marine,
    sd = a.ssa.marine.sd
  )

  # Freshwater alignment properties
  alpha.freshwater.samples <- stats::rnorm(
    n_sim,
    mean = alpha.freshwater,
    sd = alpha.freshwater.sd
  )
  a.sa.freshwater.samples <- stats::rnorm(
    n_sim,
    mean = a.sa.freshwater,
    sd = a.sa.freshwater.sd
  )
  a.v.freshwater.samples <- stats::rnorm(
    n_sim,
    mean = a.v.freshwater,
    sd = a.v.freshwater.sd
  )
  a.m.freshwater.samples <- stats::rnorm(
    n_sim,
    mean = a.m.freshwater,
    sd = a.m.freshwater.sd
  )
  a.ssa.freshwater.samples <- stats::rnorm(
    n_sim,
    mean = a.ssa.freshwater,
    sd = a.ssa.freshwater.sd
  )

  # Simulated body length and intercept
  sim_beta_log10_body_length_samples <- stats::rnorm(
    n_sim,
    mean = beta_log10_body_length,
    sd = se_beta_log10_body_length
  )
  sim_body_length_intercept_samples <- stats::rnorm(
    n_sim,
    mean = body_length_intercept,
    sd = se_body_length_intercept
  )

  # Tissue translocation length calculations
  sim_beta_0 <- stats::rnorm(n_sim * 1.2, mean = beta_0, sd = se_beta_0)
  sim_beta_1 <- stats::rnorm(n_sim * 1.2, mean = beta_1, sd = se_beta_1)
  sim_X50 <- -sim_beta_0 / sim_beta_1
  upper.tissue.trans.size.um.samples <- sim_X50[sim_X50 > 0][1:n_sim]

  # Create parameter values data frame
  param_values <- data.frame(
    alpha.marine = alpha.marine.samples,
    a.sa.marine = a.sa.marine.samples,
    a.v.marine = a.v.marine.samples,
    a.m.marine = a.m.marine.samples,
    a.ssa.marine = a.ssa.marine.samples,
    R.ave.water.marine = R.ave.water.marine.samples,
    H_W_ratio.marine = H_W_ratio.marine_samples,
    alpha.freshwater = alpha.freshwater.samples,
    a.sa.freshwater = a.sa.freshwater.samples,
    a.v.freshwater = a.v.freshwater.samples,
    a.m.freshwater = a.m.freshwater.samples,
    a.ssa.freshwater = a.ssa.freshwater.samples,
    R.ave.water.freshwater = R.ave.water.freshwater.samples,
    H_W_ratio.freshwater = H_W_ratio.freshwater_samples,
    sim_beta_log10_body_length = sim_beta_log10_body_length_samples,
    sim_body_length_intercept = sim_body_length_intercept_samples,
    upper.tissue.trans.size.um = upper.tissue.trans.size.um.samples
  )

  return(param_values)
}


### Function to align data to ERMs
align_data <- function(
  df, #ToMEx database
  x1M_set_input = 1, #um lower size for all alignments
  x1D_set_input = 1, #um lower size for all alignments
  x2D_set_input = 5000, #um
  upper.tissue.trans.size.um_input = 88,
  R.ave.marine_input = 0.77, # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
  H_W_ratio.marine_input = 0.77, # H:W ratio assumed same as width:length ratio (Kooi et al. 2021)
  R.ave.freshwater_input = 0.67, # average length to width ratio of microplastics in freshwater environment (Kooi et al. 2021)
  H_W_ratio.freshwater_input = 0.67, # H:W ratio assumed same as width:length ratio (Kooi et al. 2021)
  R.ave.sediment.marine_input = 0.75, # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
  R.ave.sediment.freshwater_input = 0.70, # average length to width ratio of microplastics in freshwater environment (Kooi et al. 2021)
  p.ave.marine_input = 1.10, #average density in marine surface water
  alpha.marine_input = 2.07, #table s4 for marine surface water. length
  a.sa.marine_input = 1.50, #marine surface area power law
  a.v.marine_input = 1.48, #a_V for marine surface water volume
  a.m.marine_input = 1.32, # upper limit fora_m for mass for marine surface water in table S4
  a.ssa.marine_input = 1.98, # A_SSA for marine surface water
  p.ave.freshwater_input = 1.04, #average density in freshwater surface water
  alpha.freshwater_input = 2.64, #table s4 for freshwater surface water. length
  a.sa.freshwater_input = 2.00, #freshwater surface area power law
  a.v.freshwater_input = 1.68, #a_V for freshwater surface water volume
  a.m.freshwater_input = 1.65, # upper limit fora_m for mass for freshwater surface water in table S4
  a.ssa.freshwater_input = 2.71, # A_SSA for freshwater surface water
  beta_log10_body_length_input = 0.9341, # Jâms, et al 2020 Nature paper
  body_length_intercept_input = 1.1200 # Jâms, et al 2020 Nature paper
) {
  #10 #um #set size for x2M)

  # Check if columns exist and conditionally create them
  if (!"dose.mg.kg.sed.measured" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(dose.mg.kg.sed.measured = measured.dose.mg.kg.sediment)
  }

  if (!"dose.mg.kg.sed.nominal" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(dose.mg.kg.sed.nominal = nominal.dose.mg.kg.sediment)
  }

  if (!"dose.particles.kg.sed.nominal" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        dose.particles.kg.sed.nominal = nominal.dose.particles.kg.sediment
      )
  }

  if (!"environment" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(environment = env_f)
  }

  df <- df %>%
    dplyr::ungroup() %>% #if data grouped - bad things happen
    #assign user  input values
    dplyr::mutate(
      x1M_set = x1M_set_input,
      x1D_set = x1D_set_input,
      x2D_set = x2D_set_input,
      upper.tissue.trans.size.um = upper.tissue.trans.size.um_input,
      beta_log10_body_length = beta_log10_body_length_input,
      body_length_intercept = body_length_intercept_input,
      H_W_ratio.marine = H_W_ratio.marine_input,
      H_W_ratio.freshwater = H_W_ratio.freshwater_input,
      R.ave.marine = R.ave.marine_input,
      R.ave.freshwater = R.ave.freshwater_input,
      R.ave.sediment.marine = R.ave.sediment.marine_input,
      R.ave.sediment.freshwater = R.ave.sediment.freshwater_input,
      p.ave.marine = p.ave.marine_input,
      alpha.marine = alpha.marine_input,
      a.sa.marine = a.sa.marine_input,
      a.v.marine = a.v.marine_input,
      a.m.marine = a.m.marine_input,
      a.ssa.marine = a.ssa.marine_input,
      p.ave.freshwater = p.ave.freshwater_input,
      alpha.freshwater = alpha.freshwater_input,
      a.sa.freshwater = a.sa.freshwater_input,
      a.v.freshwater = a.v.freshwater_input,
      a.m.freshwater = a.m.freshwater_input,
      a.ssa.freshwater = a.ssa.freshwater_input
    )

  # calculate ERM for each species
  aoc_final <- df %>%
    # define environment-specific alpha parameters #
    dplyr::mutate(
      alpha = dplyr::case_when(
        environment == "Marine" ~ alpha.marine,
        environment == "Freshwater" ~ alpha.freshwater
      ),
      a.sa = dplyr::case_when(
        environment == "Marine" ~ a.sa.marine,
        environment == "Freshwater" ~ a.sa.freshwater
      ),
      a.v = dplyr::case_when(
        environment == "Marine" ~ a.v.marine,
        environment == "Freshwater" ~ a.v.freshwater
      ),
      a.m = dplyr::case_when(
        environment == "Marine" ~ a.m.marine,
        environment == "Freshwater" ~ a.m.freshwater
      ),
      a.ssa = dplyr::case_when(
        environment == "Marine" ~ a.ssa.marine,
        environment == "Freshwater" ~ a.ssa.freshwater
      ),
      R.ave = dplyr::case_when(
        environment == "Marine" ~ R.ave.marine,
        environment == "Freshwater" ~ R.ave.freshwater
      ),
      H_W_ratio = dplyr::case_when(
        environment == "Marine" ~ H_W_ratio.marine,
        environment == "Freshwater" ~ H_W_ratio.freshwater
      ),
      p.ave = dplyr::case_when(
        environment == "Marine" ~ p.ave.marine,
        environment == "Freshwater" ~ p.ave.freshwater
      )
    ) %>%

    ###############################################################################
    ###### Determine bioaccesible fractions for polydisperse particle mixtures ####
    ######################################################################################
    # define upper size length for ingestion
    dplyr::mutate(
      x2M_ingest = dplyr::case_when(
        is.na(max.size.ingest.um) ~ x2D_set,
        max.size.ingest.um < x2D_set ~ max.size.ingest.um,
        max.size.ingest.um > x2D_set ~ x2D_set
      )
    ) %>% #set to 5,000 as upper limit or max size ingest, whichever is smaller
    # define upper size length for Translocation
    #set to 88um for upper limit or max size ingest, whichever is smaller
    dplyr::mutate(
      x2M_trans = dplyr::case_when(
        is.na(max.size.ingest.um) ~ upper.tissue.trans.size.um,
        max.size.ingest.um < upper.tissue.trans.size.um ~ max.size.ingest.um,
        max.size.ingest.um >
          upper.tissue.trans.size.um ~ upper.tissue.trans.size.um
      )
    ) %>%
    # tag ingestible/translocatable for monodisperse
    dplyr::mutate(
      ingestible = dplyr::case_when(
        polydispersity == "monodisperse" &
          size.length.um.used.for.conversions <= x2M_ingest ~ "ingestible",
        polydispersity == "monodisperse" &
          size.length.um.used.for.conversions > x2M_ingest ~ "not ingestible"
      ),
      translocatable = dplyr::case_when(
        polydispersity == "monodisperse" &
          size.length.um.used.for.conversions <= x2M_trans ~ "translocatable",
        polydispersity == "monodisperse" &
          size.length.um.used.for.conversions > x2M_trans ~ "not translocatable"
      )
    ) %>%

    ## assign whether polydisperse data are partially, fully, or not bioavailable
    dplyr::mutate(
      ingestible_poly = dplyr::case_when(
        polydispersity == "polydisperse" &
          size.length.max.um.used.for.conversions <= x2M_ingest &
          size.length.min.um.used.for.conversions <=
            x2M_ingest ~ "ingestible (all)",
        polydispersity == "polydisperse" &
          size.length.max.um.used.for.conversions > x2M_ingest &
          size.length.min.um.used.for.conversions <=
            x2M_ingest ~ "ingestible (some)",
        polydispersity == "polydisperse" &
          size.length.max.um.used.for.conversions > x2M_ingest &
          size.length.min.um.used.for.conversions >
            x2M_ingest ~ "not ingestible"
      ),
      translocatable_poly = dplyr::case_when(
        polydispersity == "polydisperse" &
          size.length.max.um.used.for.conversions <= x2M_trans &
          size.length.min.um.used.for.conversions <=
            x2M_trans ~ "translocatable (all)",
        polydispersity == "polydisperse" &
          size.length.max.um.used.for.conversions > x2M_trans &
          size.length.min.um.used.for.conversions <=
            x2M_trans ~ "translocatable (some)",
        polydispersity == "polydisperse" &
          size.length.max.um.used.for.conversions > x2M_trans &
          size.length.min.um.used.for.conversions >
            x2M_trans ~ "not translocatable"
      )
    ) %>%
    # For the partially ingestible/translocatable study, we prepare this data for alignment using a two-step process, in which we first re-calculate #   # the effect concentration (particles/volume) using the Correction Factor equation (Koelmans et al. 2019):
    ####### STEP 1: Re-Calculate Dose  for ingestible/translocatable fractions ####
    # correct for partially translocatable particles
    dplyr::mutate(
      CF_bioavailable_trans = dplyr::case_when(
        translocatable_poly == "translocatable (some)" ~ CFfnx(
          a = alpha,
          x1D = size.length.min.um.used.for.conversions,
          x2D = x2M_trans,
          x1M = size.length.min.um.used.for.conversions,
          x2M = size.length.max.um.used.for.conversions
        ),
        T ~ 1
      )
    ) %>% # all other cases retain original dose
    # now correct the dosage (will be fraction )
    dplyr::mutate(
      dose.particles.mL.trans = dplyr::case_when(
        translocatable_poly == "translocatable (some)" ~ CF_bioavailable_trans *
          dose.particles.mL.master,
        T ~ dose.particles.mL.master
      )
    ) %>%
    # correct for partially ingestible particles
    dplyr::mutate(
      CF_bioavailable_ingest = dplyr::case_when(
        ingestible_poly == "ingestible (some)" ~ CFfnx(
          a = alpha,
          x1D = size.length.min.um.used.for.conversions,
          x2D = x2M_ingest,
          x1M = size.length.min.um.used.for.conversions,
          x2M = size.length.max.um.used.for.conversions
        ),
        T ~ 1
      )
    ) %>%
    dplyr::mutate(
      dose.particles.mL.ingest = dplyr::case_when(
        ingestible_poly == "ingestible (some)" ~ CF_bioavailable_ingest *
          dose.particles.mL.master,
        T ~ dose.particles.mL.master
      )
    ) %>%

    ##### STEP 2: re-assign the min/max sizes of the particle distributions to those that are actually bioavailable within the exposure mixture,             ## labelling them accordingly for use in translocation or food dilution-associated ERM calculations.
    ##### ----- LENGTH ------ ###
    # no need to correct monodisperse. Min for polydispserse remains same #
    ## polydisperse ##
    dplyr::mutate(
      size.length.max.um.trans = dplyr::case_when(
        translocatable_poly == "translocatable (some)" ~ x2M_trans,
        T ~ size.length.max.um.used.for.conversions
      ),
      size.length.max.um.ingest = dplyr::case_when(
        ingestible_poly == "ingestible (some)" ~ x2M_ingest,
        T ~ size.length.max.um.used.for.conversions
      )
    ) %>%
    ##### ----- WIDTH ------ ###
    ## Monodisperse ##
    dplyr::mutate(
      size.width.um.used.for.conversions = dplyr::case_when(
        is.na(size.width.um.used.for.conversions) & shape_f == "Fiber" ~ 15, # assume 15 um width for fibers unless already known (kooi et al. 2021)
        is.na(size.width.um.used.for.conversions) &
          shape_f == "Sphere" ~ size.length.um.used.for.conversions, # W = L for spheres
        is.na(size.width.um.used.for.conversions) &
          shape_f == "Fragment" ~ size.length.um.used.for.conversions * R.ave, #use average width:length ratio for fragments
        is.na(size.width.um.used.for.conversions) &
          shape_f == "Not Reported" ~ size.length.um.used.for.conversions *
          R.ave, #Assume fragment
        T ~ size.width.um.used.for.conversions # if available, use as-is
      )
    ) %>%
    ### Polydisperse ###
    # Min is always same #
    # calculate size parameters using compartment characteristics
    dplyr::mutate(
      size.width.min.um.used.for.conversions = dplyr::case_when(
        shape_f == "sphere" ~ size.length.min.um.used.for.conversions, #all dims same
        shape_f == "fiber" ~ R.ave * size.length.min.um.used.for.conversions, #median holds for all particles (Kooi et al 2021)
        shape_f == "Not Reported" ~ R.ave *
          size.length.min.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
        shape_f == "fragment" ~ R.ave * size.length.min.um.used.for.conversions
      )
    ) %>% # average width to length ratio in the marine environment (kooi et al 2021)
    ### Max depends on ingest/trans limits ###
    # TRANS #
    dplyr::mutate(
      size.width.max.um.trans = dplyr::case_when(
        is.na(size.width.max.um.used.for.conversions) & shape_f == "Fiber" ~ 15, # assume 15 um width for fibers unless already known (kooi et al. 2021)
        shape_f == "Sphere" ~ size.length.max.um.trans, # W = L for spheres
        shape_f == "Fragment" ~ size.length.max.um.trans * R.ave, #use average width:length ratio for fragments
        T ~ size.width.max.um.used.for.conversions # if available, use as-is (fibers only)
      )
    ) %>%
    # INGEST #
    dplyr::mutate(
      size.width.max.um.ingest = dplyr::case_when(
        is.na(size.width.max.um.used.for.conversions) & shape_f == "Fiber" ~ 15, # assume 15 um width for fibers unless already known (kooi et al. 2021)
        shape_f == "Sphere" ~ size.length.max.um.ingest, # W = L for spheres
        shape_f == "Fragment" ~ size.length.max.um.ingest * R.ave, #use average width:length ratio for fragments
        T ~ size.width.max.um.used.for.conversions # if available, use as-is (fibers only)
      )
    ) %>%
    ###### ------ HEIGHT ----- #####
    ## Monodisperse ##
    #estimate height based on shape (data doesn't exist in ToMEx for monodisperse, because never reported)
    dplyr::mutate(
      size.height.um.used.for.conversions = dplyr::case_when(
        shape_f == "Sphere" ~ size.length.um.used.for.conversions, # if spherical, height = length
        shape_f != "Sphere" ~ size.width.um.used.for.conversions * H_W_ratio # if not spherical, height = width * H:W ratio
      )
    ) %>%
    ### Polydisperse ##
    ## Min is always same ##
    dplyr::mutate(
      size.height.min.um.used.for.conversions = dplyr::case_when(
        shape_f == "Sphere" ~ size.length.min.um.used.for.conversions, # if spherical, height = length
        shape_f != "Sphere" ~ size.width.min.um.used.for.conversions * H_W_ratio # if not spherical, height = width * H:W ratio
      )
    ) %>% # environment AND average height to width ratio (kooi et al 2021)
    # trans #
    dplyr::mutate(
      size.height.max.um.trans = dplyr::case_when(
        shape_f == "Sphere" ~ size.width.max.um.trans, # if spherical, height = length
        shape_f != "Sphere" ~ size.width.max.um.trans * H_W_ratio # if not spherical, height = width * H:W ratio
      )
    ) %>%
    # Ingest #
    dplyr::mutate(
      size.height.max.um.ingest = dplyr::case_when(
        shape_f == "Sphere" ~ size.width.max.um.ingest, # if spherical, height = length
        shape_f != "Sphere" ~ size.width.max.um.ingest * H_W_ratio # if not spherical, height = width * H:W ratio
      )
    ) %>%
    ############ ------ Volume ------ ##########
    ###### re-calculate size, surface area, volume, mass based on user-defined R.ave ####
    #### Monodisperse ##
    # calculate volume for monodisperse particles #
    dplyr::mutate(
      particle.volume.um3 = volumefnx(
        R = R.ave,
        length = size.length.um.used.for.conversions,
        width = size.width.um.used.for.conversions,
        height = size.height.um.used.for.conversions
      )
    ) %>%
    #### Polydisperse ##
    dplyr::mutate(
      particle.volume.um3.min = volumefnx(
        R = R.ave,
        length = size.length.min.um.used.for.conversions,
        width = size.width.min.um.used.for.conversions,
        height = size.height.min.um.used.for.conversions
      )
    ) %>%
    ### Trans ##
    # calculate min and max volume when polydisperse particles are used (being sure to use ingestion-restricted sizes)
    # calculate max volume when polydisperse particles are used (translocation-limited)
    dplyr::mutate(
      particle.volume.um3.max.trans = volumefnx(
        R = R.ave,
        length = size.length.max.um.trans,
        width = size.width.max.um.trans,
        height = size.height.max.um.trans
      )
    ) %>%
    ### Ingest  ##
    # calculate max volume when polydisperse particles are used (ingestlocation-limited)
    dplyr::mutate(
      particle.volume.um3.max.ingest = volumefnx(
        R = R.ave,
        length = size.length.max.um.ingest,
        width = size.width.max.um.ingest,
        height = size.height.max.um.ingest
      )
    ) %>%
    ############ ------ Surface Area ------ ##########
    # calculate surface are for monodisperse particles
    dplyr::mutate(
      particle.surface.area.um2 = SAfnx(
        length = size.length.um.used.for.conversions,
        width = size.width.um.used.for.conversions,
        height = size.height.um.used.for.conversions,
        R = R.ave,
        H_W_ratio = H_W_ratio
      )
    ) %>%
    ##### Polydisperse ###
    # calculate min/max SA for polydisperse mixtures (being sure to use translocation/ingestion-restricted polydisperse upper sizes)
    dplyr::mutate(
      particle.surface.area.um2.min = SAfnx(
        length = size.length.min.um.used.for.conversions,
        width = size.width.min.um.used.for.conversions,
        height = size.height.min.um.used.for.conversions,
        R = R.ave,
        H_W_ratio = H_W_ratio
      )
    ) %>%
    ### Trans ##
    dplyr::mutate(
      particle.surface.area.um2.max.trans = SAfnx(
        R = R.ave,
        H_W_ratio = H_W_ratio,
        length = size.length.max.um.trans,
        width = size.width.max.um.trans,
        height = size.height.max.um.trans
      )
    ) %>%
    ### Ingest ###
    dplyr::mutate(
      particle.surface.area.um2.max.ingest = SAfnx(
        R = R.ave,
        H_W_ratio = H_W_ratio,
        length = size.length.max.um.ingest,
        width = size.width.max.um.ingest,
        height = size.height.max.um.ingest
      )
    ) %>%
    #calculate mass for monodisperse particles
    dplyr::mutate(
      mass.per.particle.mg = massfnx(
        v = particle.volume.um3,
        p = density.g.cm3
      ) *
        1e-3
    ) %>% #equation uses g/cm3
    #calculate minimum and maximum mass or polydisperse particles
    dplyr::mutate(
      mass.per.particle.mg.min = massfnx(
        v = particle.volume.um3.min,
        p = density.g.cm3
      ) *
        1e-3
    ) %>% #equation uses g/cm3
    # Trans
    dplyr::mutate(
      mass.per.particle.mg.max.trans = massfnx(
        v = particle.volume.um3.max.trans,
        p = density.g.cm3
      ) *
        1e-3
    ) %>% #equation uses g/cm3
    # Ingest
    dplyr::mutate(
      mass.per.particle.mg.max.ingest = massfnx(
        v = particle.volume.um3.max.ingest,
        p = density.g.cm3
      ) *
        1e-3
    ) %>% #equation uses g/cm3
    ################################################################
    ########################## ALIGNMENTS ##########################
    ###################################################################
    dplyr::mutate(mu.p.mono = 1) %>% #mu_x_mono is always 1 for particles to particles
    #### TISSUE TRANSLOCATION ####
    # calculate effect threshold for particles
    dplyr::mutate(EC_mono_p.particles.mL_trans = dose.particles.mL.trans) %>%
    dplyr::mutate(
      mu.p.poly_trans = mux_polyfnx(
        a.x = alpha,
        x_UL = x2M_trans, #upper translocatable size limit (width of particle)
        x_LL = x1M_set
      )
    ) %>%
    # polydisperse effect threshold for particles
    dplyr::mutate(
      EC_poly_p.particles.mL_trans = (EC_mono_p.particles.mL_trans *
        mu.p.mono) /
        mu.p.poly_trans
    ) %>%
    #calculate CF_bio for all conversions
    dplyr::mutate(
      CF_bio_trans = CFfnx(
        x1M = x1M_set, #lower size bin
        x2M = x2M_trans, #upper translocatable
        x1D = x1D_set, #default
        x2D = x2D_set, #default
        a = alpha
      )
    ) %>%
    ## Calculate environmentally relevant effect threshold for particles
    dplyr::mutate(
      EC_env_p.particles.mL_trans = EC_poly_p.particles.mL_trans * CF_bio_trans
    ) %>% #aligned particle effect concentraiton (1-5000 um)

    #### Surface area ERM ####
    ##--- environmental calculations ---###
    #calculate lower translocatable surface area
    dplyr::mutate(
      x_LL_sa_trans = SAfnx(
        length = x1D_set,
        width = x1D_set,
        height = x1D_set
      ),
      #calculate upper translocatable surface area using spherical assumption
      x_UL_sa_trans = SAfnx(
        length = x2M_trans,
        width = x2M_trans,
        height = x2M_trans
      )
    ) %>%
    #calculate mu_x_poly (env) for surface area
    dplyr::mutate(
      mu.sa.poly_trans = mux_polyfnx(a.sa, x_UL_sa_trans, x_LL_sa_trans)
    ) %>%

    ##--- laboratory calculations ---###
    ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
    #(note that if mixed particles were used, a different equation must be used)
    dplyr::mutate(
      mu.sa.mono.trans = dplyr::case_when(
        polydispersity == "monodisperse" ~ particle.surface.area.um2, # use reported surface area in monodisperse
        polydispersity == "polydisperse" ~ mux_polyfnx(
          a.x = a.sa,
          x_LL = particle.surface.area.um2.min,
          x_UL = particle.surface.area.um2.max.trans
        )
      )
    ) %>%

    #calculate polydisperse effect concentration for surface area (particles/mL)
    dplyr::mutate(
      EC_poly_sa.particles.mL_trans = (EC_mono_p.particles.mL_trans *
        mu.sa.mono.trans) /
        mu.sa.poly_trans
    ) %>%
    #calculate environmentally realistic effect threshold
    dplyr::mutate(
      EC_env_sa.particles.mL_trans = EC_poly_sa.particles.mL_trans *
        CF_bio_trans
    ) %>%

    ##### FOOD DILUTION ####
    # calculate effect threshold for particles
    dplyr::mutate(EC_mono_p.particles.mL_ingest = dose.particles.mL.ingest) %>%
    dplyr::mutate(
      mu.p.poly_ingest = mux_polyfnx(
        a.x = alpha, #alpha for particles
        x_UL = x2M_ingest, #upper ingestible size limit
        x_LL = x1M_set
      )
    ) %>%
    # polydisperse effect threshold for particles
    dplyr::mutate(
      EC_poly_p.particles.mL_ingest = (EC_mono_p.particles.mL_ingest *
        mu.p.mono) /
        mu.p.poly_ingest
    ) %>%
    #calculate CF_bio for all conversions
    dplyr::mutate(
      CF_bio_ingest = CFfnx(
        x1M = x1M_set, #lower size bin
        x2M = x2M_ingest, #upper ingestible length
        x1D = x1D_set, #default
        x2D = x2D_set, #default upper size range
        a = alpha
      )
    ) %>%
    ## Calculate environmentally relevant effect threshold for particles
    dplyr::mutate(
      EC_env_p.particles.mL_ingest = EC_poly_p.particles.mL_ingest *
        CF_bio_ingest
    ) %>% #aligned particle effect concentraiton (1-5000 um)
    #### volume ERM ####
    ##--- environmental calculations ---###
    #calculate lower ingestible volume (assumed spherical)
    dplyr::mutate(
      x_LL_v_ingest = volumefnx(
        length = x1D_set,
        width = x1D_set,
        height = x1D_set
      ),
      # max ingestible volume assumed to be spherical
      x_UL_v_ingest = volumefnx(
        length = x2M_ingest,
        width = x2M_ingest,
        height = x2M_ingest
      )
    ) %>%
    # calculate mu.v.poly
    dplyr::mutate(
      mu.v.poly_ingest = mux_polyfnx(a.v, x_UL_v_ingest, x_LL_v_ingest)
    ) %>%
    ##--- laboratory calculations ---###
    ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
    #(note that if mixed particles were used, a different equation must be used)
    dplyr::mutate(
      mu.v.mono.ingest = dplyr::case_when(
        polydispersity == "monodisperse" ~ particle.volume.um3, # use reported volume in monodisperse
        polydispersity == "polydisperse" ~ mux_polyfnx(
          a.x = a.v,
          x_LL = particle.volume.um3.min,
          x_UL = particle.volume.um3.max.ingest
        )
      )
    ) %>%

    #calculate polydisperse effect concentration for volume (particles/mL)
    dplyr::mutate(
      EC_poly_v.particles.mL_ingest = (EC_mono_p.particles.mL_ingest *
        mu.v.mono.ingest) /
        mu.v.poly_ingest
    ) %>%
    #calculate environmentally realistic effect threshold
    dplyr::mutate(
      EC_env_v.particles.mL_ingest = EC_poly_v.particles.mL_ingest *
        CF_bio_ingest
    ) %>%

    ###### CLEANUP #####
    dplyr::mutate(
      translocatable = ifelse(
        size.length.um.used.for.conversions > x2M_trans,
        "not translocatable",
        "translocatable"
      )
    ) %>%
    dplyr::mutate(
      ingestible = ifelse(
        size.length.um.used.for.conversions > x2M_ingest,
        "not ingestible",
        "ingestible"
      )
    ) %>%
    ## collapse poly/mono bioavailbilities
    dplyr::mutate(
      ingestible = dplyr::case_when(
        !is.na(ingestible_poly) ~ ingestible_poly,
        T ~ ingestible
      ),
      translocatable = dplyr::case_when(
        !is.na(translocatable_poly) ~ translocatable_poly,
        T ~ translocatable
      )
    ) %>%
    #rowwise() %>%
    dplyr::mutate(unique_id = row_number()) %>%
    ##### EASY ID ###
    dplyr::mutate(
      particles.mL.ox.stress = EC_env_sa.particles.mL_trans,
      particles.mL.food.dilution = EC_env_v.particles.mL_ingest
    ) %>%
    dplyr::mutate(
      particles.mL.food.dilution = dplyr::case_when(
        Group == "Algae" ~ NA,
        T ~ particles.mL.food.dilution
      )
    ) %>%
    # ensure particles are within valid range
    dplyr::filter(size.length.um.used.for.conversions >= x1D_set)

  aoc_final
}


##### generate data used in MC-SIM (differs slightly from 'generate data' funciton above as it uses LHS)
### Dynamic Parameters (default from Kooi et al. 2021 and ToMEx 2.0 ecotox risk paper)
param_default_values <- data.frame(
  x1D_set = 1,
  x2D_set = 5000,
  R.ave.water.marine = 0.77,
  R.ave.water.marine.sd = 0.29,
  R.ave.water.freshwater = 0.67,
  R.ave.water.freshwater.sd = 0.28,
  R.ave.sediment.marine = 0.75,
  R.ave.sediment.marine.sd = 0.30,
  R.ave.sediment.freshwater = 0.70,
  R.ave.sediment.freshwater.sd = 0.33,
  p.ave.marine = 1.10,
  p.ave.marine.sd = 0.14,
  alpha.marine = 2.07,
  alpha.marine.sd = 0.03,
  a.sa.marine = 1.50,
  a.sa.marine.sd = 0.009,
  a.v.marine = 1.48,
  a.v.marine.sd = 0.063,
  a.m.marine = 1.32,
  a.m.marine.sd = 0.009,
  a.ssa.marine = 1.98,
  a.ssa.marine.sd = 0.297,
  p.ave.freshwater = 1.04,
  p.ave.freshwater.sd = 0.12,
  alpha.freshwater = 2.64,
  alpha.freshwater.sd = 0.01,
  a.sa.freshwater = 2.00,
  a.sa.freshwater.sd = 0.065,
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

#' Generate Latin hypercube samples for alignment parameters
#'
#' Builds Sobol-based parameter matrices used by the Monte Carlo alignment
#' routine. Defaults mirror the values reported in the ToMEx 2.0 manuscript.
#'
#' @param n_sobol Integer number of Sobol samples to generate.
#' @param params Data frame of parameter means/standard deviations.
#' @param upper.tissue.truncation.limit Maximum tissue translocation size (um).
#' @param x1M_set Lower bound for monodisperse particle length (um).
#' @param x2D_set Upper bound for polydisperse particle length (um).
#'
#' @return A data frame with sampled parameters and a `simulation_id` column.
#' @export
matrix_function <- function(
  n_sobol = 10,
  params = param_default_values,
  upper.tissue.truncation.limit = 500,
  x1M_set = 1,
  x2D_set = 5000
) {
  matrices <- c("A", "B", "AB", "BA")
  first <- total <- "azzini"

  param_names <- c(
    "alpha.marine",
    "a.sa.marine",
    "a.v.marine",
    "a.m.marine",
    "alpha.freshwater",
    "a.sa.freshwater",
    "a.v.freshwater",
    "a.m.freshwater",
    "R.ave.water.marine",
    "R.ave.water.freshwater",
    "H_W_ratio.marine",
    "H_W_ratio.freshwater",
    #  "R.ave.sediment.marine", "R.ave.sediment.freshwater", "a.ssa.marine", "a.ssa.freshwater",
    "sim.beta.log10.body.length",
    "sim.body.length.intercept",
    "upper.tissue.trans.size.um"
  )

  # Generate the Sobol' sequence
  mat <- sensobol::sobol_matrices(N = n_sobol, params = param_names, type = "LHS")

  # Convert to data.table
  # Convert to data.table
  mat <- data.table::as.data.table(mat)

  # Transform each column to the specified probability distribution
  mat[,
    alpha.marine := stats::rnorm(
      .N,
      mean = params$alpha.marine,
      sd = params$alpha.marine.sd
    )
  ]
  mat[,
    a.sa.marine := stats::rnorm(
      .N,
      mean = params$a.sa.marine,
      sd = params$a.sa.marine.sd
    )
  ]
  mat[,
    a.v.marine := stats::rnorm(
      .N,
      mean = params$a.v.marine,
      sd = params$a.v.marine.sd
    )
  ]
  mat[,
    a.m.marine := stats::rnorm(
      .N,
      mean = params$a.m.marine,
      sd = params$a.m.marine.sd
    )
  ]
  #mat[, a.ssa.marine := rnorm(.N,  mean = a_ssa.marine, sd = a.ssa.sd.marine)]
  mat[,
    alpha.freshwater := stats::rnorm(
      .N,
      mean = params$alpha.freshwater,
      sd = params$alpha.freshwater.sd
    )
  ]
  mat[,
    a.sa.freshwater := stats::rnorm(
      .N,
      mean = params$a.sa.freshwater,
      sd = params$a.sa.freshwater.sd
    )
  ]
  mat[,
    a.v.freshwater := stats::rnorm(
      .N,
      mean = params$a.v.freshwater,
      sd = params$a.v.freshwater.sd
    )
  ]
  mat[,
    a.m.freshwater := stats::rnorm(
      .N,
      mean = params$a.m.freshwater,
      sd = params$a.m.freshwater.sd
    )
  ]
  #mat[, a.ssa.freshwater := rnorm(.N,   mean = a_ssa.freshwater, sd = a.ssa.sd.freshwater)]
  mat[,
    R.ave.water.marine := truncnorm::rtruncnorm(
      .N,
      a = 0.0001,
      b = 1,
      mean = params$R.ave.water.marine,
      sd = params$R.ave.water.marine.sd
    )
  ]
  mat[,
    H_W_ratio.marine := truncnorm::rtruncnorm(
      .N,
      a = 0.0001,
      b = 1,
      mean = params$R.ave.water.marine,
      sd = params$R.ave.water.marine.sd
    )
  ]
  mat[,
    R.ave.water.freshwater := truncnorm::rtruncnorm(
      .N,
      a = 0.0001,
      b = 1,
      mean = params$R.ave.water.freshwater,
      sd = params$R.ave.water.freshwater.sd
    )
  ]
  mat[,
    H_W_ratio.freshwater := truncnorm::rtruncnorm(
      .N,
      a = 0.0001,
      b = 1,
      mean = params$R.ave.water.freshwater,
      sd = params$R.ave.water.freshwater.sd
    )
  ]
  mat[,
    sim.beta.log10.body.length := stats::rnorm(
      .N,
      mean = params$beta_log10_body_length,
      sd = params$se_beta_log10_body_length
    )
  ]
  mat[,
    sim.body.length.intercept := stats::rnorm(
      .N,
      mean = params$body_length_intercept,
      sd = params$se_beta_log10_body_length
    )
  ]
  #mat[, "R.ave.sediment.marine" := rtruncnorm(.N, a = 0.0001, b= 0.9999, mean = 0.75, sd = 0.30)]
  #mat[, "R.ave.sediment.freshwater" := rtruncnorm(.N, a = 0.0001, b= 0.9999, mean = 0.70, sd = 0.33)]

  # Parameters from tissue translocaiton logistic regression model derived in Coffin et al. (2026)
  # beta_0 <- 1.308344
  # beta_1 <- -0.01468148
  # se_beta_0 <- 0.3963612
  # se_beta_1 <- 0.006657993

  # Simulate beta_0 and beta_1
  sim_beta_0 <- stats::rnorm(
    nrow(mat) * 1.4,
    mean = params$beta_0,
    sd = params$se_beta_0
  )
  sim_beta_1 <- stats::rnorm(
    nrow(mat) * 1.4,
    mean = params$beta_1,
    sd = params$se_beta_1
  )

  # Calculate X50 for each simulation
  sim_X50 <- -sim_beta_0 / sim_beta_1

  # Truncate distribution to not fall below 0
  upper.tissue.trans.size.um_samples <- sim_X50 %>%
    data.frame() %>%
    # filter invalid values
    dplyr::filter(. > x1M_set & . < x2D_set) %>% # truncate distribution between default sizes (1 and 5,000 um)
    dplyr::filter(. < upper.tissue.truncation.limit) %>% #truncate distribution below biologically plausible limit
    dplyr::slice(1:nrow(mat))

  upper.tissue.trans.size.um_samples <- as.numeric(
    upper.tissue.trans.size.um_samples$.
  )

  mat[, upper.tissue.trans.size.um := upper.tissue.trans.size.um_samples]

  # Convert the data.table to a data.frame
  mat <- as.data.frame(mat) %>%
    dplyr::mutate(simulation_id = paste0("sim", row_number()))

  return(mat)
}


#### MC-SIM Align Data and generate thresholds (preps data for PSSD++ AND runs MC-SIM SSD method)
model_wrapper_sobol_parallel <- function(params, simulation_id) {
  # perform alignments using parameters for simulation
  aoc_MC_iter <- align_data(
    aoc_aligned,
    upper.tissue.trans.size.um_input = params$upper.tissue.trans.size.um[1],
    body_length_intercept_input = params$sim_body_length_intercept[1],
    beta_log10_body_length_input = params$sim_beta_log10_body_length[1],
    R.ave.marine_input = params$R.ave.water.marine[1],
    H_W_ratio.marine_input = params$H_W_ratio.marine[1],
    alpha.marine_input = params$alpha.marine[1],
    a.sa.marine_input = params$a.sa.marine[1],
    a.v.marine_input = params$a.v.marine[1],
    a.m.marine_input = params$a.m.marine[1],
    R.ave.freshwater_input = params$R.ave.water.freshwater[1],
    H_W_ratio.freshwater_input = params$H_W_ratio.freshwater[1],
    alpha.freshwater_input = params$alpha.freshwater[1],
    a.sa.freshwater_input = params$a.sa.freshwater[1],
    a.v.freshwater_input = params$a.v.freshwater[1],
    a.m.freshwater_input = params$a.m.freshwater[1]
  )

  # Filter out risk criteria and calculate thresholds for different environments
  aoc_aligned <- aoc_MC_iter %>%
    tidyr::drop_na(effect.metric) %>%
    dplyr::filter(tier_zero_tech_f == "Red Criteria Passed")

  # Calculate thresholds for marine and freshwater
  marine_thresholds <- process_environment_data(
    aoc_aligned,
    "Marine",
    upper.tissue.trans.size.um = sim.upper.tissue.trans.size.um,
    x1D_set = x1D_set,
    x2D_set = x2D_set
  ) %>%
    ## convert particles/L thresholds to volume, mass, surface area (ERM-dependent)
    convert_units_fxn(environment = "Marine", params = params[1, ])

  freshwater_thresholds <- process_environment_data(
    aoc_aligned,
    "Freshwater",
    upper.tissue.trans.size.um = sim.upper.tissue.trans.size.um,
    x1D_set = x1D_set,
    x2D_set = x2D_set
  ) %>%
    ## convert particles/L thresholds to volume, mass, surface area (ERM-dependent)
    convert_units_fxn(environment = "Freshwater", params = params[1, ])

  # Store results in list
  sobol_result <- list(
    tox_vals = list(
      simulation_id = simulation_id,
      particles_mL_ox_stress = aoc_aligned$particles.mL.ox.stress,
      particles_mL_food_dilution = aoc_aligned$particles.mL.food.dilution,
      ingestible = aoc_aligned$ingestible,
      ingestible_poly = aoc_aligned$ingestible_poly,
      translocatable = aoc_aligned$translocatable,
      translocatable_poly = aoc_aligned$translocatable_poly,
      species = aoc_aligned$Species,
      group = aoc_aligned$Group,
      environment = aoc_aligned$environment,
      size.length.um.used.for.conversions = aoc_aligned$size.length.um.used.for.conversions,
      polydispersity = aoc_aligned$polydispersity,
      shape_f = aoc_aligned$shape_f,
      poly_f = aoc_aligned$poly_f,
      shape_f = aoc_aligned$shape_f,
      doi = aoc_aligned$doi,
      effect.metric = aoc_aligned$effect.metric,
      af.noec = aoc_aligned$af.noec,
      af.time = aoc_aligned$af.time,
      bio_f = aoc_aligned$bio_f,
      risk.13 = aoc_aligned$risk.13,
      dose.particles.mL.master = aoc_aligned$dose.particles.mL.master,
      unique_id = aoc_aligned$unique_id
    ), #row id that allows matching back to database
    base_thresholds = list(
      simulation_id = simulation_id,
      marine = marine_thresholds,
      freshwater = freshwater_thresholds
    )
  )

  return(sobol_result)
}


####### parallel process alignments (NO SSDs - just used to prep data for PSSD++)#####
## wrapper function for Monte Carlo simulations
MC_sim_align_wrapper <- function(tox_data, params, simulation_id) {
  # perform alignments using parameters for simulation
  aligned_tox_data <- align_data(
    tox_data,
    upper.tissue.trans.size.um_input = params$upper.tissue.trans.size.um[1],
    body_length_intercept_input = params$sim_body_length_intercept[1],
    beta_log10_body_length_input = params$sim_beta_log10_body_length[1],
    R.ave.marine_input = params$R.ave.water.marine[1],
    H_W_ratio.marine_input = params$H_W_ratio.marine[1],
    alpha.marine_input = params$alpha.marine[1],
    a.sa.marine_input = params$a.sa.marine[1],
    a.v.marine_input = params$a.v.marine[1],
    a.m.marine_input = params$a.m.marine[1],
    R.ave.freshwater_input = params$R.ave.water.freshwater[1],
    H_W_ratio.freshwater_input = params$H_W_ratio.freshwater[1],
    alpha.freshwater_input = params$alpha.freshwater[1],
    a.sa.freshwater_input = params$a.sa.freshwater[1],
    a.v.freshwater_input = params$a.v.freshwater[1],
    a.m.freshwater_input = params$a.m.freshwater[1]
  )

  # Store results in list
  result <- list(
    tox_vals = list(
      simulation_id = simulation_id,
      particles_mL_ox_stress = aligned_tox_data$particles.mL.ox.stress,
      particles_mL_food_dilution = aligned_tox_data$particles.mL.food.dilution,
      ingestible = aligned_tox_data$ingestible,
      ingestible_poly = aligned_tox_data$ingestible_poly,
      translocatable = aligned_tox_data$translocatable,
      translocatable_poly = aligned_tox_data$translocatable_poly,
      species = aligned_tox_data$Species,
      group = aligned_tox_data$Group,
      environment = aligned_tox_data$environment,
      size.length.um.used.for.conversions = aligned_tox_data$size.length.um.used.for.conversions,
      polydispersity = aligned_tox_data$polydispersity,
      shape_f = aligned_tox_data$shape_f,
      poly_f = aligned_tox_data$poly_f,
      shape_f = aligned_tox_data$shape_f,
      doi = aligned_tox_data$doi,
      effect.metric = aligned_tox_data$effect.metric,
      af.noec = aligned_tox_data$af.noec,
      af.time = aligned_tox_data$af.time,
      bio_f = aligned_tox_data$bio_f,
      risk.13 = aligned_tox_data$risk.13,
      dose.particles.mL.master = aligned_tox_data$dose.particles.mL.master,
      unique_id = aligned_tox_data$unique_id #row id that allows matching back to database
    )
  )
  return(result)
}

# Parallel process function
#' Monte Carlo alignments in parallel
#'
#' Runs the alignment workflow across Sobol-sampled parameter sets, returning
#' simulation-ready microplastic ecotoxicity data.
#'
#' @param tox_data Data frame of toxicity records (defaults to the bundled dataset).
#' @param param_matrix Parameter matrix created by \link[=matrix_function]{matrix_function()}.
#' @param x1D_set Lower particle-size bound (um) for alignment.
#' @param x2D_set Upper particle-size bound (um) for alignment.
#' @param n_sim Number of alignment simulations to run.
#' @param num_cores Number of workers (`\"auto\"` uses all cores minus one).
#' @param suppress_warnings Logical; silence alignment warnings.
#'
#' @return A data frame with aligned doses for each simulation.
#' @export
MC_sim_align_parallel <- function(
  tox_data = aoc_risk_paper, # data to be processed (raw ToMEx 2.0 as default)
  param_matrix = mat, # matrix params for simulation
  x1D_set = 1, # minimum particle size
  x2D_set = 5000, # maximum particle size
  n_sim = 10,
  num_cores = "auto",
  suppress_warnings = TRUE
) {
  if (suppress_warnings) {
    old_warn <- getOption("warn")
    options(warn = -1)
    on.exit(options(warn = old_warn), add = TRUE)
  }

  # assign cores
  if (num_cores == "auto") {
    num_cores <- max(1, parallel::detectCores() - 1)
  }
  # make cluster
  cl <- parallel::makeCluster(num_cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  doSNOW::registerDoSNOW(cl)
  # list number of cores used
  message(crayon::yellow(paste0(
    "Using ",
    num_cores,
    " cores for parallel processing"
  )))
  # Set up a progress bar
  pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  on.exit(try(close(pb), silent = TRUE), add = TRUE)
  # Initialize vector to store sobol results
  results_list <- vector("list", n_sim)
  # Initialize vector to store sobol results
  results_parallel <- vector("list", n_sim)
  tictoc::tic("Parallel processing alignments")

  x1D_set = x1D_set
  x2D_set = x2D_set
  tox_data = tox_data

  cat(crayon::yellow(paste0(
    "Running alignments for ",
    n_sim,
    " simulation(s)"
  )))
  # Run alignments in parallel
  results_list <- foreach::foreach(
    i = 1:n_sim,
    .packages = c(
      "magrittr",
      "dplyr",
      "tidyr",
      "purrr",
      "foreach",
      "doParallel",
      "truncnorm",
      "sensobol"
    ),
    .options.snow = opts,
    .export = c(
      "MC_sim_align_wrapper",
      "tox_data",
      "align_data",
      "generate_structure_checks",
      "mux_polyfnx",
      "CFfnx",
      "volumefnx",
      "SAfnx",
      "massfnx",
      "SSA.inversefnx",
      "n_sim",
      "x1D_set",
      "x2D_set",
      "param_matrix"
    ),
    .combine = rbind
  ) %dopar%
    {
      param_matrix <- param_matrix[i, ]
      # Ensure all extracted parameters are correctly coerced to numeric
      # Coerce to numeric directly while accessing the first element of potential list
      alpha.marine <- as.numeric(param_matrix$alpha.marine[1])
      a.sa.marine <- as.numeric(param_matrix$a.sa.marine[1])
      a.v.marine <- as.numeric(param_matrix$a.v.marine[1])
      a.m.marine <- as.numeric(param_matrix$a.m.marine[1])
      # a.ssa.marine <- as.numeric(param_matrix$a.ssa.marine[1])
      alpha.freshwater <- as.numeric(param_matrix$alpha.freshwater[1])
      a.sa.freshwater <- as.numeric(param_matrix$a.sa.freshwater[1])
      a.v.freshwater <- as.numeric(param_matrix$a.v.freshwater[1])
      a.m.freshwater <- as.numeric(param_matrix$a.m.freshwater[1])
      #  a.ssa.freshwater <- as.numeric(param_matrix$a.ssa.freshwater[1])
      R_ave_water_marine <- as.numeric(param_matrix$R.ave.water.marine[1])
      R_ave_water_freshwater <- as.numeric(param_matrix$R.ave.water.freshwater[
        1
      ])
      #  R_ave_sediment_marine <- as.numeric(param_matrix$R.ave.sediment.marine[1])
      #  R_ave_sediment_freshwater <- as.numeric(param_matrix$R.ave.sediment.freshwater[1])
      sim_beta_log10_body_length <- as.numeric(param_matrix$sim.beta.log10.body.length[
        1
      ])
      sim_body_length_intercept <- as.numeric(param_matrix$sim.body.length.intercept[
        1
      ])
      sim.upper.tissue.trans.size.um <- as.numeric(param_matrix$upper.tissue.trans.size.um[
        1
      ])
      simulation_id <- as.character(param_matrix$simulation_id[1])

      # Perform the model wrapping
      results_parallel[[i]] <- MC_sim_align_wrapper(
        tox_data = tox_data,
        param_matrix,
        simulation_id = simulation_id
      )
    }

  elapsed <- tictoc::toc()
  message(crayon::green(paste0(
    "Probabilistic alignments completed in ",
    round(elapsed$toc - elapsed$tic, 2),
    " sec"
  )))
  # Close progress bar and stop the cluster
  close(pb)
  parallel::stopCluster(cl)

  # Separate individual toxicity values from Sims in MC_results #
  results_df <- do.call(
    rbind,
    lapply(results_list, function(x) {
      data.frame(
        simulation_id = x$simulation_id,
        unique_id = factor(x$unique_id), #toxicity data point (row)
        particles_L_ox_stress = x$particles_mL_ox_stress * 1000,
        particles_L_food_dilution = x$particles_mL_food_dilution * 1000,
        species = x$species,
        group = x$group,
        shape_f = x$shape_f,
        poly_f = x$poly_f,
        environment = x$environment,
        effect.metric = x$effect.metric,
        bio_f = x$bio_f,
        risk.13 = x$risk.13,
        af.time = x$af.time,
        af.noec = x$af.noec,
        ingestible = factor(x$ingestible),
        ingestible_poly = x$ingestible_poly,
        translocatable = factor(x$translocatable),
        translocatable_poly = factor(x$translocatable_poly)
      )
    })
  ) %>%
    dplyr::mutate(
      ingestible = dplyr::case_when(
        !is.na(ingestible_poly) ~ ingestible_poly,
        T ~ ingestible
      ),
      translocatable = dplyr::case_when(
        !is.na(translocatable_poly) ~ translocatable_poly,
        T ~ translocatable
      )
    ) |>
    dplyr::rename(Group = group, Species = species) %>%
    dplyr::mutate(
      sim_id_numeric = as.numeric(str_replace(simulation_id, "sim", ""))
    )
  return(results_df)
} # close MC_sim_align_parallel function


## PSSD++ wrapper function ###
generate_plots_and_summary <- function(
  tier,
  environment,
  erm,
  results_df,
  quantile_type = 7,
  color_palette = global_color_palette,
  sim = 300, # Default number of simulations
  num_iterations = 300, # Default number of iterations
  cv_dp = NULL,
  cv_uf = 0.5, # Default value for coefficient of variation for uncertainty factors
  rmore_method = "step", # step or lognormal
  species_data_source,
  output_path,
  presentation_path,
  debug = FALSE,
  debug_dir = "../output/pssd_debug",
  silent = FALSE, # progress reporting for single-threading
  progress = NULL # report-out progress for parallel processing
) {
  combo_id <- paste0("Tier", tier, "_", environment, "_", erm)

  if (debug) {
    if (!dir.exists(debug_dir)) {
      dir.create(debug_dir, recursive = TRUE)
    }
    cat("\n🐞 DEBUG MODE ON for:", combo_id, "\n")
  }

  # ---------- STEP 0: Save raw input ----------
  if (debug) {
    saveRDS(
      list(
        tier = tier,
        environment = environment,
        erm = erm,
        results_df = results_df,
        species_data_source = species_data_source,
        quantile_type = quantile_type,
        color_palette = color_palette,
        sim = sim,
        num_iterations = num_iterations,
        cv_dp = cv_dp,
        cv_uf = cv_uf,
        rmore_method = rmore_method
      ),
      file = file.path(debug_dir, paste0("step0_inputs_", combo_id, ".rds"))
    )
  }

  # ---------- STEP 1: Prepare data ----------
  MCdf <- tryCatch(
    {
      if (debug) {
        cat("🔧 STEP 1: prep_data()...\n")
      }
      x <- prep_data(results_df)
      if (debug) {
        cat("  MCdf names:", paste(names(x), collapse = ", "), "\n")
        saveRDS(
          x,
          file = file.path(debug_dir, paste0("step1_MCdf_", combo_id, ".rds"))
        )
      }
      x
    },
    error = function(e) {
      if (debug) {
        saveRDS(
          list(
            error = e,
            results_df = results_df
          ),
          file = file.path(
            debug_dir,
            paste0("ERROR_step1_prep_data_", combo_id, ".rds")
          )
        )
        cat("❌ ERROR @ STEP 1 (prep_data):", e$message, "\n")
      }
      stop(e)
    }
  )

  # ---------- STEP 2: Run PSSD analysis ----------
  pssd_results <- tryCatch(
    {
      if (debug) {
        cat("🔧 STEP 2: run_pSSD_analysis()...\n")
      }
      x <- run_pSSD_analysis(
        data_matrices = MCdf,
        num_iterations = num_iterations,
        sim = sim,
        cv_dp = cv_dp,
        cv_uf = cv_uf,
        rmore_method = rmore_method,
        data_name = paste0("tier", tier, "_", erm),
        silent = silent, # silencing for parallel processing
        progress = progress # report out for parallel processing
      )
      if (debug) {
        cat("  pSSD dims:", dim(x$pSSD), "\n")
        saveRDS(
          x,
          file = file.path(
            debug_dir,
            paste0("step2_pssd_results_", combo_id, ".rds")
          )
        )
      }
      x
    },
    error = function(e) {
      if (debug) {
        saveRDS(
          list(
            error = e,
            MCdf = MCdf
          ),
          file = file.path(
            debug_dir,
            paste0("ERROR_step2_pssd_", combo_id, ".rds")
          )
        )
        cat("❌ ERROR @ STEP 2 (run_pSSD_analysis):", e$message, "\n")
      }
      stop(e)
    }
  )

  # ---------- STEP 3: Prepare data for ggplot ----------
  plot_prep <- tryCatch(
    {
      if (debug) {
        cat("🔧 STEP 3: prepare_plot_data()...\n")
      }
      x <- prepare_plot_data(
        pSSD = pssd_results$pSSD,
        original_data = results_df,
        data_matrices = MCdf,
        quantile_type = quantile_type,
        species_data_source = species_data_source,
        debug = debug
      )
      if (debug) {
        cat("  quantiles_df rows:", nrow(x$quantiles_df), "\n")
        cat("  species_data rows:", nrow(x$species_data), "\n")
        cat("  all_NOEC rows:", nrow(x$all_NOEC), "\n")
        saveRDS(
          x,
          file = file.path(
            debug_dir,
            paste0("step3_plot_prep_", combo_id, ".rds")
          )
        )
      }
      x
    },
    error = function(e) {
      if (debug) {
        saveRDS(
          list(
            error = e,
            pSSD = pssd_results$pSSD,
            MCdf = MCdf,
            results_df = results_df
          ),
          file = file.path(
            debug_dir,
            paste0("ERROR_step3_prepare_plot_data_", combo_id, ".rds")
          )
        )
        cat("❌ ERROR @ STEP 3 (prepare_plot_data):", e$message, "\n")
      }
      stop(e)
    }
  )

  # ---------- STEP 4a: Generate PSSD plot ----------
  pSSD_plot <- tryCatch(
    {
      if (debug) {
        cat("🔧 STEP 4a: pSSD_plot_fnx()...\n")
      }
      x <- pSSD_plot_fnx(
        quantiles_df = plot_prep$quantiles_df,
        species_data = plot_prep$species_data,
        all_NOEC = plot_prep$all_NOEC,
        tier_name = paste("Tier", tier),
        environment_name = environment,
        ERM_name = erm,
        color_palette = color_palette,
        debug = debug,
        debug_dir = debug_dir,
        combo_id = combo_id
      )
      if (debug) {
        saveRDS(
          x,
          file = file.path(
            debug_dir,
            paste0("step4a_pSSD_plot_", combo_id, ".rds")
          )
        )
      }
      x
    },
    error = function(e) {
      if (debug) {
        saveRDS(
          list(
            error = e,
            quantiles_df = plot_prep$quantiles_df,
            species_data = plot_prep$species_data,
            all_NOEC = plot_prep$all_NOEC
          ),
          file = file.path(
            debug_dir,
            paste0("ERROR_step4a_pSSD_plot_", combo_id, ".rds")
          )
        )
        cat("❌ ERROR @ STEP 4a (pSSD_plot_fnx):", e$message, "\n")
      }
      stop(e)
    }
  )

  # ---------- STEP 4b: PNEC plots ----------
  PNEC_05 <- apply(pssd_results$pSSD, 2, function(x) {
    quantile(x, probs = 0.05, type = quantile_type)
  })
  PNEC_05_df <- data.frame(PNEC = PNEC_05)
  if (debug) {
    saveRDS(
      PNEC_05_df,
      file = file.path(debug_dir, paste0("step4b_PNEC05_df_", combo_id, ".rds"))
    )
  }

  PNEC_plot_05 <- tryCatch(
    {
      if (debug) {
        cat("🔧 STEP 4b: make_PNEC_plot() HC5...\n")
      }
      p <- make_PNEC_plot(
        pssd_results = pssd_results,
        hcx = 0.05,
        quantile_type = quantile_type,
        data_name = paste("Tier", tier, "(", erm, ") - HC5")
      ) +
        ggplot2::labs(
          title = "Tier 3 (HC5)",
          x = "PNEC (Particles/L; 1 to 5,000 μm)"
        ) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank())
      if (debug) {
        saveRDS(
          p,
          file = file.path(
            debug_dir,
            paste0("step4b_PNEC_plot_05_", combo_id, ".rds")
          )
        )
      }
      p
    },
    error = function(e) {
      if (debug) {
        saveRDS(
          list(error = e, pssd_results = pssd_results),
          file = file.path(
            debug_dir,
            paste0("ERROR_step4b_PNEC_plot_05_", combo_id, ".rds")
          )
        )
        cat("❌ ERROR @ STEP 4b (make_PNEC_plot HC5):", e$message, "\n")
      }
      stop(e)
    }
  )

  PNEC_plot_10 <- tryCatch(
    {
      if (debug) {
        cat("🔧 STEP 4b: make_PNEC_plot() HC10...\n")
      }
      p <- make_PNEC_plot(
        pssd_results = pssd_results,
        hcx = 0.10,
        quantile_type = quantile_type,
        data_name = paste("Tier", tier, "(", erm, ") - HC10")
      ) +
        ggplot2::labs(
          title = "Tier 4 (HC10)",
          x = "PNEC (Particles/L; 1 to 5,000 μm)"
        ) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank()
        )
      if (debug) {
        saveRDS(
          p,
          file = file.path(
            debug_dir,
            paste0("step4b_PNEC_plot_10_", combo_id, ".rds")
          )
        )
      }
      p
    },
    error = function(e) {
      if (debug) {
        saveRDS(
          list(error = e, pssd_results = pssd_results),
          file = file.path(
            debug_dir,
            paste0("ERROR_step4b_PNEC_plot_10_", combo_id, ".rds")
          )
        )
        cat("❌ ERROR @ STEP 4b (make_PNEC_plot HC10):", e$message, "\n")
      }
      stop(e)
    }
  )

  # ---------- STEP 5: PNEC summaries ----------
  PNEC_summary_05 <- PNEC_data_summary(
    pssd_results = pssd_results,
    hcx = 0.05,
    quantile_type = quantile_type,
    data_name = paste("Tier", tier, "(", erm, ") - HC5")
  )
  PNEC_summary_10 <- PNEC_data_summary(
    pssd_results = pssd_results,
    hcx = 0.10,
    quantile_type = quantile_type,
    data_name = paste("Tier", tier, "(", erm, ") - HC10")
  )
  if (debug) {
    saveRDS(
      list(
        PNEC_summary_05 = PNEC_summary_05,
        PNEC_summary_10 = PNEC_summary_10
      ),
      file = file.path(
        debug_dir,
        paste0("step5_PNEC_summaries_", combo_id, ".rds")
      )
    )
  }

  # ---------- STEP 6: Arrange plots ----------
  arranged_plot <- tryCatch(
    {
      if (debug) {
        cat("🔧 STEP 6: ggarrange()...\n")
      }
      p <- ggpubr::ggarrange(
        pSSD_plot,
        ggpubr::ggarrange(
          PNEC_plot_05,
          PNEC_plot_10,
          ncol = 2,
          widths = c(0.5, 0.5)
        ),
        nrow = 2,
        heights = c(0.6, 0.4)
      )
      if (debug) {
        saveRDS(
          p,
          file = file.path(
            debug_dir,
            paste0("step6_arranged_plot_", combo_id, ".rds")
          )
        )
      }
      p
    },
    error = function(e) {
      if (debug) {
        saveRDS(
          list(
            error = e,
            pSSD_plot = pSSD_plot,
            PNEC_plot_05 = PNEC_plot_05,
            PNEC_plot_10 = PNEC_plot_10
          ),
          file = file.path(
            debug_dir,
            paste0("ERROR_step6_arranged_plot_", combo_id, ".rds")
          )
        )
        cat("❌ ERROR @ STEP 6 (ggarrange):", e$message, "\n")
      }
      stop(e)
    }
  )

  # ---------- STEP 7: Save figures ----------
  tryCatch(
    {
      if (debug) {
        cat("💾 STEP 7: ggsave() plot...\n")
      }
      ggplot2::ggsave(
        filename = paste0(erm, "_", environment, "_arranged_plot.jpg"),
        dpi = 300,
        path = output_path,
        plot = arranged_plot,
        width = 10,
        height = 8,
        units = "in"
      )
    },
    error = function(e) {
      if (debug) {
        saveRDS(
          list(error = e, arranged_plot = arranged_plot),
          file = file.path(
            debug_dir,
            paste0("ERROR_step7_ggsave_plot_", combo_id, ".rds")
          )
        )
        cat("❌ ERROR @ STEP 7 (ggsave plot):", e$message, "\n")
      }
      stop(e)
    }
  )

  # ---------- RETURN ----------
  if (debug) {
    cat("✅ generate_plots_and_summary() complete for", combo_id, "\n")
  }

  return(list(
    pSSD_plot = pSSD_plot,
    PNEC_plot_05 = PNEC_plot_05,
    PNEC_plot_10 = PNEC_plot_10,
    arranged_plot = arranged_plot,
    summary_05 = PNEC_summary_05,
    summary_10 = PNEC_summary_10,
    quantiles = plot_prep$quantiles_df,
    PNEC_df = PNEC_05_df
  ))
}


# Dynamically generate a color palette based on the levels of Group
generate_color_palette <- function(data) {
  unique_groups <- unique(data$Group) # Extract unique levels of Group
  #palette <- cols4all::c4a("seaborn.dark", n = length(unique_groups)) # Generate colors
  palette <- cols4all::c4a("cols4all.area9d", n = length(unique_groups)) # Generate colors
  names(palette) <- unique_groups # Assign names to the colors
  return(palette)
}


pSSD_plot_fnx <- function(
  quantiles_df,
  all_NOEC,
  species_data,
  data_name = NULL,
  environment_name,
  tier_name,
  ERM_name,
  color_palette = color_palette,
  debug = FALSE,
  debug_dir = "../output/pssd_debug",
  combo_id = "UNKNOWN"
) {
  # ========== DEBUG SETUP =====================================================
  if (debug) {
    if (!dir.exists(debug_dir)) {
      dir.create(debug_dir, recursive = TRUE)
    }
    cat("\n📊 DEBUG: pSSD_plot_fnx() for", combo_id, "\n")
    cat("  Species_data rows:", nrow(species_data), "\n")
    cat("  all_NOEC rows:", nrow(all_NOEC), "\n")
    cat("  N NA NOEC values:", sum(is.na(all_NOEC$NOEC_median)), "\n")
    cat("  N Inf NOEC values:", sum(is.infinite(all_NOEC$NOEC_median)), "\n")
    cat(
      "  Range median NOEC:",
      paste(range(all_NOEC$NOEC_median, na.rm = TRUE), collapse = " to "),
      "\n"
    )

    saveRDS(
      list(
        quantiles_df = quantiles_df,
        all_NOEC = all_NOEC,
        species_data = species_data
      ),
      file = file.path(debug_dir, paste0("pSSD_plot_inputs_", combo_id, ".rds"))
    )
  }

  # ========== X LIMITS =========================================================
  max_species_value <- max(species_data$NOEC_log10, na.rm = TRUE)
  x_limits <- c(
    min(c(quantiles_df$iv, species_data$NOEC_log10), na.rm = TRUE) - 1,
    max_species_value + 6
  )

  if (debug) {
    cat("  x_limits:", paste(x_limits, collapse = " to "), "\n")
    below <- species_data %>% dplyr::filter(NOEC_log10 < x_limits[1])
    above <- species_data %>% dplyr::filter(NOEC_log10 > x_limits[2])
    cat("  N labels < x_min:", nrow(below), "\n")
    cat("  N labels > x_max:", nrow(above), "\n")
    saveRDS(
      list(below = below, above = above),
      file = file.path(
        debug_dir,
        paste0("pSSD_label_x_checks_", combo_id, ".rds")
      )
    )
  }

  # ========== LEGEND LABELS + GROUPS ==========================================
  legend_labels <- setNames(
    paste0(
      "<span style='color:",
      color_palette,
      "'>",
      names(color_palette),
      "</span>"
    ),
    names(color_palette)
  )

  # only groups that actually appear in the data AND exist in the palette
  groups_present <- sort(unique(na.omit(all_NOEC$Group)))
  groups_present <- intersect(groups_present, names(color_palette))

  legend_labels_present <- legend_labels[groups_present]

  if (debug) {
    cat("  Groups present:", paste(groups_present, collapse = ", "), "\n")
  }

  # ========== SAFE STAGE WRAPPER ==============================================
  safe_stage <- function(plot_obj, stage_label) {
    if (!debug) {
      return(plot_obj)
    }
    cat("    🧱 Building stage:", stage_label, "...\n")
    tryCatch(
      {
        ggplot2::ggplot_build(plot_obj)
        saveRDS(
          plot_obj,
          file = file.path(
            debug_dir,
            paste0("pSSD_plot_stage_", stage_label, "_", combo_id, ".rds")
          )
        )
        cat("    ✅ Stage OK:", stage_label, "\n")
        plot_obj
      },
      error = function(e) {
        saveRDS(
          list(error = e, plot = plot_obj),
          file = file.path(
            debug_dir,
            paste0("ERROR_pSSD_plot_stage_", stage_label, "_", combo_id, ".rds")
          )
        )
        cat("    ❌ ERROR in stage", stage_label, ":", e$message, "\n")
        stop(e)
      }
    )
  }

  # ========== BASE =============================================================
  p <- ggplot2::ggplot()
  p <- safe_stage(p, "base")

  # ========== POINTS ===========================================================
  p <- p +
    ggplot2::geom_point(
      data = all_NOEC,
      ggplot2::aes(x = NOEC_001th, y = Prop, color = Group, fill = Group),
      size = 1,
      alpha = 0.03
    ) +
    ggplot2::geom_point(
      data = all_NOEC,
      ggplot2::aes(x = NOEC_01th, y = Prop, color = Group, fill = Group),
      size = 1.5,
      alpha = 0.07
    ) +
    ggplot2::geom_point(
      data = all_NOEC,
      ggplot2::aes(x = NOEC_1th, y = Prop, color = Group, fill = Group),
      size = 2,
      alpha = 0.1
    ) +
    ggplot2::geom_point(
      data = all_NOEC,
      ggplot2::aes(x = NOEC_10th, y = Prop, color = Group, fill = Group),
      size = 2.5,
      alpha = 0.2
    ) +
    ggplot2::geom_point(
      data = all_NOEC,
      ggplot2::aes(x = NOEC_median, y = Prop, color = Group, fill = Group),
      size = 3,
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = all_NOEC,
      ggplot2::aes(x = NOEC_90th, y = Prop, color = Group, fill = Group),
      size = 2.5,
      alpha = 0.2
    ) +
    ggplot2::geom_point(
      data = all_NOEC,
      ggplot2::aes(x = NOEC_99th, y = Prop, color = Group, fill = Group),
      size = 2,
      alpha = 0.1
    ) +
    ggplot2::geom_point(
      data = all_NOEC,
      ggplot2::aes(x = NOEC_999th, y = Prop, color = Group, fill = Group),
      size = 1.5,
      alpha = 0.07
    ) +
    ggplot2::geom_point(
      data = all_NOEC,
      ggplot2::aes(x = NOEC_9999th, y = Prop, color = Group, fill = Group),
      size = 1,
      alpha = 0.03
    )

  p <- safe_stage(p, "points")

  # ========== LABELS LOWER =====================================================
  p <- p +
    ggrepel::geom_label_repel(
      data = species_data %>% dplyr::filter(Prop < 0.6),
      ggplot2::aes(x = NOEC_log10, y = Prop, label = Species, color = Group),
      direction = "y",
      hjust = 0,
      nudge_x = 7,
      fontface = "italic",
      size = 3,
      fill = "white",
      alpha = 0.85,
      segment.alpha = 0.4,
      max.overlaps = Inf,
      max.time = 5,
      max.iter = 1e7,
      na.rm = TRUE
    )

  p <- safe_stage(p, "labels_lower")

  # ========== LABELS UPPER =====================================================
  p <- p +
    ggrepel::geom_label_repel(
      data = species_data %>% dplyr::filter(Prop >= 0.6),
      ggplot2::aes(x = NOEC_log10, y = Prop, label = Species, color = Group),
      direction = "y",
      hjust = 1,
      nudge_x = -7,
      fontface = "italic",
      size = 3,
      fill = "white",
      alpha = 0.85,
      segment.alpha = 0.4,
      max.overlaps = Inf,
      max.time = 5,
      max.iter = 1e7,
      na.rm = TRUE
    )

  p <- safe_stage(p, "labels_upper")

  # ========== MEAN ECDF LINE ===================================================
  p <- p +
    ggplot2::geom_line(
      data = quantiles_df,
      ggplot2::aes(x = iv, y = mean),
      color = "black",
      linewidth = 1
    )

  p <- safe_stage(p, "meanline")

  # ========== SCALES + THEME (NO VECTOR override.aes) ==========================
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = pretty(x_limits, n = 6),
      labels = scales::math_format(10^.x),
      limits = x_limits
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::ylab("Cumulative probability") +
    ggplot2::xlab("NOEC (particles/L; 1 to 5,000 μm)") +

    # Global palette, but legend only includes groups_present
    ggplot2::scale_color_manual(
      values = color_palette,
      breaks = groups_present,
      labels = legend_labels_present,
      drop = FALSE,
      name = "Organism Group"
    ) +
    ggplot2::scale_fill_manual(
      values = color_palette,
      breaks = groups_present,
      labels = legend_labels_present,
      drop = FALSE,
      name = "Organism Group",
      guide = "none"
    ) +

    # IMPORTANT: override.aes does NOT contain per-key vectors anymore
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          shape = 21,
          size = 3,
          alpha = 1
        )
      )
    ) +

    ggplot2::labs(
      title = "Microplastics PSSD++",
      subtitle = paste0(ERM_name, " (", environment_name, ")")
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      legend.text = ggtext::element_markdown(size = 12)
    )

  p <- safe_stage(p, "scales_theme")

  if (debug) {
    cat("✅ pSSD_plot_fnx() finished for", combo_id, "\n")
  }

  return(p)
}

# Calculate PNEC1 if not already calculated
make_PNEC_plot <- function(
  pssd_results,
  hcx,
  data_name,
  quantile_type = 7 #Type 1 is no interpolation (default for pSSD), # type  7 inovlves interpolation (note that type 7 results in higher values for mean/median, but the mode is around the same as the median would be with type 1. It's necessary to use type 7 to get distinct values between HC5 and HC10 for freshwater bc there's so few species)
) {
  PNEC <- apply(pssd_results$pSSD, 2, function(x) {
    quantile(
      x,
      probs = hcx, #HC value
      type = quantile_type
    )
  })
  # Convert PNEC1 to a dataframe for ggplot
  PNEC_df <- data.frame(PNEC = PNEC)

  # Calculate specific percentiles
  median_val <- quantile(PNEC, probs = 0.5, type = quantile_type, na.rm = TRUE)
  p001_val <- quantile(PNEC, probs = 0.01, type = quantile_type, na.rm = TRUE)
  p5_val <- quantile(PNEC, probs = 0.05, type = quantile_type, na.rm = TRUE)
  p95_val <- quantile(PNEC, probs = 0.95, type = quantile_type, na.rm = TRUE)
  p999_val <- quantile(
    PNEC,
    probs = 0.999999999999,
    type = quantile_type,
    na.rm = TRUE
  ) # 99.9th percentile for x-axis limits

  # Define the x-axis range dynamically
  x_range <- p999_val - p001_val # Replace max_x and min_x with your actual x-axis limits

  # Create a data frame for the percentile labels
  label_data <- data.frame(
    x = c(p5_val, median_val, (p95_val - 0.001 * x_range)),
    #  y = c(0.00001, 0.00001, 0.00001),  # Adjust y position for labels
    label = c(
      paste0("**5th %** <br>", signif(p5_val, 2)),
      paste0("**Median** <br>", signif(median_val, 2)),
      paste0("**95th %** <br>", signif(p95_val, 2))
    )
  )

  # Plotting
  # Plotting
  ggplot2::ggplot(PNEC_df, ggplot2::aes(x = PNEC)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      fill = grDevices::rgb(255, 127, 80, max = 255, alpha = 50),
      color = "coral",
      bins = 350,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      x = "PNEC (Particles/L)",
      y = "Probability density",
      title = "Probability density of the PNEC",
      subtitle = data_name
    ) +
    ggplot2::scale_x_log10(limits = c(p001_val, p999_val)) + # Adjust x-axis limits to include 99.9th percentile

    ggplot2::geom_density(col = "red", linewidth = 1, na.rm = TRUE) + # Optionally add a density line
    # Add dotted vertical lines for the percentiles
    ggplot2::geom_vline(
      xintercept = median_val,
      linetype = "dotted",
      color = "black",
      linewidth = 1.5
    ) +
    ggplot2::geom_vline(
      xintercept = p5_val,
      linetype = "dotted",
      color = "black",
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = p95_val,
      linetype = "dotted",
      color = "black",
      linewidth = 1
    ) +
    # Add text labels for the percentiles above the distribution
    # Add text labels for the percentiles
    ggtext::geom_richtext(
      data = label_data,
      aes(x = x, label = label),
      y = 0.01,
      size = 5,
      color = "black",
      vjust = -1,
      hjust = -0.01,
      fill = NA,
      label.color = NA
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      #axis.text.y = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_blank()
    )
}

prepare_plot_data <- function(
  original_data = results_df_food_marine,
  pSSD = tier1_2_food_pssd$pSSD,
  data_matrices = tier1_2_food_matrices,
  species_data_source = tier1_2_food,
  quantile_type = 7,
  debug = FALSE
) {
  # Ensure that the required matrices are present and correctly formatted
  if (
    is.null(data_matrices[['DP']]) ||
      is.null(data_matrices[['UFdd']]) ||
      is.null(data_matrices[['UFt']])
  ) {
    stop("Data matrices for DP, UFdd, or UFt are not available or incorrect.")
  }

  # Optional Debugging
  if (debug) {
    cat("\n🔧 DEBUG: prepare_plot_data() input checks...\n")
    cat("  DP dimensions:", dim(data_matrices[['DP']]), "\n")
    cat("  UFdd dimensions:", dim(data_matrices[['UFdd']]), "\n")
    cat("  UFt dimensions:", dim(data_matrices[['UFt']]), "\n")
    cat("  pSSD matrix dimensions:", dim(pSSD), "\n")
    cat("  original_data rows:", nrow(original_data), "\n")

    bad1 <- sum(is.na(original_data$dose_new_particles_L))
    bad2 <- sum(original_data$af.noec == 0)
    bad3 <- sum(original_data$af.time == 0)

    cat("  Missing dose values:", bad1, "\n")
    cat("  Zero af.noec:", bad2, "\n")
    cat("  Zero af.time:", bad3, "\n")
  }

  # Define interpolation values
  iv <- seq(-5, 10, 0.001)
  ECDF.data <- matrix(NA, ncol(pSSD), length(iv))

  # Compute ECDF for each column in pSSD
  # for (i in 1:10000) {
  for (i in 1:ncol(pSSD)) {
    # uncomment if not beta testing!
    the.ecdf.f <- ecdf(log(pSSD[, i], base = 10))
    ECDF.data[i, ] <- the.ecdf.f(iv)
  }

  # Preparing dataframe for ggplot
  iv_vec <- seq(-5, 10, length.out = ncol(ECDF.data))
  ecdf_df <- data.frame(
    iv = rep(iv_vec, times = nrow(ECDF.data)),
    ecdf_value = as.vector(t(ECDF.data))
  )

  # Calculating quantiles for shading
  quantiles_df <- ecdf_df %>%
    dplyr::group_by(iv) %>%
    dplyr::summarise(
      min = min(ecdf_value, na.rm = TRUE),
      max = max(ecdf_value, na.rm = TRUE),
      q05 = quantile(ecdf_value, 0.05, na.rm = TRUE, type = quantile_type),
      q95 = quantile(ecdf_value, 0.95, na.rm = TRUE, type = quantile_type),
      q25 = quantile(ecdf_value, 0.25, na.rm = TRUE, type = quantile_type),
      q75 = quantile(ecdf_value, 0.75, na.rm = TRUE, type = quantile_type),
      mean = mean(ecdf_value, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  # Calculate deterministic values of NOEC
  NOEC.det <- data_matrices[['DP']] /
    (data_matrices[['UFdd']] * data_matrices[['UFt']])

  # Calculate the geometric mean of NOEC
  if (ncol(NOEC.det) > 0 && nrow(NOEC.det) > 0) {
    NOEC.gmean <- apply(NOEC.det, 2, function(x) {
      exp(mean(log(x), na.rm = TRUE))
    })
  } else {
    stop("NOEC.det matrix is empty or not properly formatted.")
  }

  ###### GEOMEAN APPROACH - Create Single point for each species ####
  # Creating a dataframe for ggplot with species data
  species_data <- data.frame(
    Species = colnames(NOEC.det),
    NOEC_log10 = log10(NOEC.gmean),
    NOEC = NOEC.gmean,
    Prop = rank(NOEC.gmean) / (length(NOEC.gmean) + 1)
  )

  # Join metadata to species data
  species_data <- species_data %>%
    dplyr::left_join(
      species_data_source %>% dplyr::distinct(Group, Species),
      by = "Species"
    ) %>%
    dplyr::left_join(
      species_data_source %>% dplyr::distinct(Species, environment),
      by = "Species"
    )

  #### Include all Points for each species ####
  # Create a dataframe with all individual unique NOEC values (collapsing alignment df)

  # Optional Debugging
  if (debug) {
    ### debugging
    cat("\n🧪 DEBUG: Calculating all_NOEC...\n")

    tmp <- original_data %>%
      dplyr::mutate(
        NOEC_raw = dose_new_particles_L / (af.noec * af.time),
        NOEC_log10 = log10(NOEC_raw)
      )

    cat("  N NOEC_raw <= 0:", sum(tmp$NOEC_raw <= 0, na.rm = TRUE), "\n")
    cat("  N NaN NOEC_log10:", sum(is.nan(tmp$NOEC_log10)), "\n")
    cat("  N Inf NOEC_log10:", sum(is.infinite(tmp$NOEC_log10)), "\n")

    bad_rows <- tmp %>%
      dplyr::filter(
        is.nan(NOEC_log10) | is.infinite(NOEC_log10) | NOEC_raw <= 0
      )

    if (nrow(bad_rows) > 0) {
      cat("  ❗ BAD ROW EXAMPLE:\n")
      print(head(bad_rows, 3))
    }
  }

  ### perform data summarization
  all_NOEC <- original_data %>%
    dplyr::mutate(
      NOEC_log10 = log10(dose_new_particles_L / (af.noec * af.time))
    ) %>%
    dplyr::group_by(
      Species,
      Group,
      poly_f,
      shape_f,
      effect.metric,
      af.noec,
      af.time,
      bio_f
    ) %>%
    dplyr::filter(is.finite(NOEC_log10)) %>% # Ensure no NAs
    # for each unique POD, get the 5th, 50th, and 95th percentile alignment value for plotting
    dplyr::summarize(
      NOEC_001th = quantile(NOEC_log10, 0.0001, type = 7),
      NOEC_01th = quantile(NOEC_log10, 0.001, type = 7),
      NOEC_1th = quantile(NOEC_log10, 0.01, type = 7),
      NOEC_10th = quantile(NOEC_log10, 0.10, type = 7),
      NOEC_mean = mean(NOEC_log10, na.rm = TRUE),
      NOEC_median = median(NOEC_log10),
      NOEC_90th = quantile(NOEC_log10, 0.90, type = 7),
      NOEC_99th = quantile(NOEC_log10, 0.99, type = 7),
      NOEC_999th = quantile(NOEC_log10, 0.999, type = 7),
      NOEC_9999th = quantile(NOEC_log10, 0.9999, type = 7),
      .groups = "drop"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      species_data %>% dplyr::distinct(Species, Prop, NOEC_log10),
      by = "Species"
    )

  if (debug) {
    n_bad <- sum(is.na(all_NOEC$NOEC_001th))
    if (n_bad > 0) {
      cat(
        "⚠ DEBUG:",
        n_bad,
        "POD groups produced NA quantiles (likely sparse Tier 3/4 data)\n"
      )
    }
  }

  if (debug) {
    n_empty <- sum(is.na(all_NOEC$NOEC_mean))
    if (n_empty > 0) {
      cat(
        "⚠ DEBUG:",
        n_empty,
        "groups had zero valid NOECs after Tier 3 stratification\n"
      )
    }
  }

  list(
    quantiles_df = quantiles_df,
    species_data = species_data,
    all_NOEC = all_NOEC
  )
}

prep_data <- function(data) {
  # Base part of the names used for returned list elements
  names <- c(
    "DP",
    "DP.SD",
    "UFt",
    "UFdd",
    "doseDescriptor",
    "polymer",
    "shape",
    "environment",
    "group"
  )

  # Data processing steps
  matrices <- list(
    DP = data %>%
      dplyr::mutate(id = row_number()) %>%
      dplyr::mutate(particles_L = dose_new_particles_L) %>% # ERM-specific - converted appropriately upstream
      dplyr::select(id, Species, particles_L) %>%
      reshape2::acast(id ~ Species, value.var = "particles_L") %>%
      Matrix::Matrix(sparse = TRUE),

    # Conditionally create DP.SD only if dose_new_particles_L_sd exists
    DP.SD = if ("dose_new_particles_L_sd" %in% colnames(data)) {
      data %>%
        dplyr::mutate(id = row_number()) %>%
        dplyr::mutate(particles_L_sd = dose_new_particles_L_sd) %>%
        dplyr::select(id, Species, particles_L_sd) %>% # ERM-specific - converted appropriately upstream
        reshape2::acast(id ~ Species, value.var = "particles_L_sd") %>%
        Matrix::Matrix(sparse = TRUE)
    } else {
      NULL
    },

    UFt = data %>%
      dplyr::mutate(id = row_number()) %>%
      dplyr::select(id, Species, af.time) %>%
      reshape2::acast(id ~ Species, value.var = "af.time") %>%
      Matrix::Matrix(sparse = T),
    UFdd = data %>%
      dplyr::mutate(id = row_number()) %>%
      dplyr::select(id, Species, af.noec) %>%
      reshape2::acast(id ~ Species, value.var = "af.noec") %>%
      Matrix::Matrix(sparse = T),
    doseDescriptor <- data %>%
      dplyr::mutate(id = row_number()) %>%
      dplyr::select(id, Species, effect.metric) %>%
      reshape2::acast(id ~ Species, value.var = "effect.metric") %>%
      Matrix::Matrix(sparse = T),
    polymer_matrix <- data %>%
      dplyr::mutate(id = row_number()) %>%
      dplyr::select(id, Species, poly_f) %>%
      reshape2::acast(id ~ Species, value.var = "poly_f") %>%
      Matrix::Matrix(sparse = T),
    shape_matrix <- data %>%
      dplyr::mutate(id = row_number()) %>%
      select(id, Species, shape_f) %>%
      reshape2::acast(id ~ Species, value.var = "shape_f") %>%
      Matrix::Matrix(sparse = T),
    env_matrix <- data %>%
      dplyr::mutate(id = row_number()) %>%
      select(id, Species, environment) %>%
      reshape2::acast(id ~ Species, value.var = "environment") %>%
      Matrix::Matrix(sparse = T),
    group_matrix <- data %>%
      dplyr::mutate(id = row_number()) %>%
      select(id, Species, Group) %>%
      reshape2::acast(id ~ Species, value.var = "Group") %>%
      Matrix::Matrix(sparse = T)
  )

  # Rename the list elements based on constructed names
  names(matrices) <- names

  # Return the named list of matrices
  matrices
}


# minimum PSSD function
run_pSSD_analysis <- function(
  data_matrices,
  num_iterations,
  sim,
  cv_dp,
  cv_uf,
  data_name,
  rmore_method = "step",
  silent = FALSE,
  progress = NULL
) {

  # Initialize an empty list to store the results
  pSSD_list <- vector("list", num_iterations)

  # Define checkpoints for reporting
  checkpoints <- floor(c(0.25, 0.5, 0.75, 1.0) * num_iterations)

  if (!silent) {
    cat(crayon::blue(sprintf(
      "🚀 Starting PSSD analysis with %d sims each...\n",
      num_iterations
    )))
  }

  # Loop through each iteration
  for (i in 1:num_iterations) {
    # Run the modified PSSD simulation
    pSSD_list[[i]] <- do.pSSD_mod(
      DP = data_matrices[["DP"]],
      DP.SD = data_matrices[["DP.SD"]],
      UFt = data_matrices[["UFt"]],
      UFdd = data_matrices[["UFdd"]],
      SIM = sim,
      CV.DP = cv_dp,
      CV.UF = cv_uf,
      rmore_method = rmore_method #lognormal or step
    )

    if (!silent) {
      # Progress reporting at checkpoints
      if (i %in% checkpoints) {
        pct_done <- round(100 * i / num_iterations)
        cat(crayon::yellow(sprintf(
          "⏳ %.0f%% complete (%d of %d simulations done)\n",
          pct_done,
          i,
          num_iterations
        )))
      }
    }
    # report out progress when running in parallel with wrapper
    if (!is.null(progress)) {
      progress()
    }
  }

  # Combine all the results into a single object
  pSSD <- do.call(cbind, pSSD_list)

  # Calculate corrected endpoints
  corr_endpoints <- data_matrices[["DP"]] /
    (data_matrices[["UFt"]] * data_matrices[["UFdd"]])

  cat(crayon::green(sprintf(
    "✅ PSSD analysis for %s complete and saved.\n",
    data_name
  )))

  return(list(pSSD = pSSD, corr_endpoints = corr_endpoints))
}

# Define dataset mapping for PSSD++ (commented out as new version below)
# get_results_df <- function(tier, environment, erm) {
#   if (environment == "Marine" && erm == "Food Dilution") {
#     if (tier < 3) {
#       return(results_df_food_marine)
#     } else {
#       return(results_df_food_t3_t4_marine)
#     }
#   } else if (environment == "Marine" && erm == "Tissue Translocation") {
#     if (tier < 3) {
#       return(results_df_tissue_marine)
#     } else {
#       return(results_df_tissue_t3_t4_marine)
#     }
#   } else if (environment == "Freshwater" && erm == "Food Dilution") {
#     if (tier < 3) {
#       return(results_df_food_freshwater)
#     } else {
#       return(results_df_food_t3_t4_freshwater)
#     }
#   } else if (environment == "Freshwater" && erm == "Tissue Translocation") {
#     if (tier < 3) {
#       return(results_df_tissue_freshwater)
#     } else {
#       return(results_df_tissue_t3_t4_freshwater)
#     }
#   } else {
#     stop("Invalid combination of tier, environment, and ERM")
#   }
# }

# ------------------------------------------------------------------
# Resolve tier → dataset key
# ------------------------------------------------------------------
resolve_tier_key <- function(tier) {
  if (tier %in% c(1, 2)) {
    return("base")
  }
  if (tier %in% c(3, 4)) {
    return("t3_t4")
  }
  stop("Unsupported tier: ", tier)
}

# ------------------------------------------------------------------
# Main data accessor used by PSSD++
# ------------------------------------------------------------------

get_results_df <- function(tier, environment, erm, erm_registry) {
  if (!erm %in% names(erm_registry)) {
    stop("Unknown ERM: ", erm)
  }

  tier_key <- resolve_tier_key(tier)

  df <- erm_registry[[erm]][[tier_key]]

  if (is.null(df)) {
    stop("No data found for ERM = ", erm, ", tier = ", tier)
  }

  df_env <- df %>%
    dplyr::filter(environment == !!environment)

  if (nrow(df_env) == 0) {
    stop(
      "Empty dataset after filtering: ",
      "ERM=",
      erm,
      ", tier=",
      tier,
      ", environment=",
      environment
    )
  }

  return(df_env)
}
# ------------------------------------------------------------------
# Sanity check all combinations before running PSSD++
# ------------------------------------------------------------------

validate_pssd_inputs <- function(tiers, environments, erms, erm_registry) {
  for (tier in tiers) {
    for (environment in environments) {
      for (erm in erms) {
        tryCatch(
          {
            df <- get_results_df(tier, environment, erm, erm_registry)
            message(
              sprintf(
                "✔ OK: Tier %d | %s | %s (n = %d)",
                tier,
                environment,
                erm,
                nrow(df)
              )
            )
          },
          error = function(e) {
            stop(
              sprintf(
                "❌ Validation failed: Tier %d | %s | %s\n%s",
                tier,
                environment,
                erm,
                e$message
              )
            )
          }
        )
      }
    }
  }
}

# Mode calculation using density
Mode_Y <- function(x) {
  dens <- density(x)
  ind <- which(dens$y == max(dens$y))
  return(dens$x[ind])
}

#PNEC Data summary function
PNEC_data_summary <- function(
  pssd_results,
  hcx,
  data_name,
  quantile_type = 7 #Type 1 is no interpolation (default for pSSD), # type  7 inovlves interpolation (note that type 7 results in higher values for mean/median, but the mode is around the same as the median would be with type 1. It's necessary to use type 7 to get distinct values between HC5 and HC10 for freshwater bc there's so few species)
) {
  PNEC <- apply(pssd_results$pSSD, 2, function(x) {
    quantile(
      x,
      probs = hcx, #HC value
      type = quantile_type
    )
  })

  #PNEC <- apply(pssd_results$pSSD, 2, function(x) quantile(x, probs = hcx, type = 1))
  PNEC_df <- data.frame(PNEC = PNEC)

  # Calculate statistics for PNEC1
  Stat_PNEC <- data.frame(
    Min = min(PNEC),
    Q5 = quantile(PNEC, 0.05, type = quantile_type),
    Q25 = quantile(PNEC, 0.25, type = quantile_type),
    Mean = mean(PNEC),
    #Median = median(PNEC),
    Median = quantile(PNEC, 0.5, type = quantile_type),
    Mode = Mode_Y(PNEC),
    Q75 = quantile(PNEC, 0.75, type = quantile_type),
    Q95 = quantile(PNEC, 0.95, type = quantile_type),
    Max = max(PNEC)
  ) %>%
    dplyr::mutate_if(is.numeric, ~ signif(., 3))

  # Transposing for similar structure to your matrix
  Stat_PNEC_t <- t(Stat_PNEC)
  colnames(Stat_PNEC_t) <- paste0(data_name, " - HC", hcx * 100)

  return(list("stats" = as.data.frame(Stat_PNEC_t), "df" = PNEC_df))
}

# function to generate PSSDs and plots with caching and optional debugging
#' Run the PSSD++ workflow across tiers, environments, and ERMs
#'
#' Generates probabilistic SSDs, summary statistics, and figures with caching
#' for each combination of tier, environment, and ERM.
#'
#' @param MC_sim_df Monte Carlo alignment output (data frame).
#' @param tiers Numeric vector of threshold tiers to evaluate.
#' @param environments Character vector of environments to include.
#' @param erms Character vector of ERM names.
#' @param erm_registry Named list of ERM-specific data frames created from
#'   the aligned data.
#' @param sim Number of PSSD simulations per combination.
#' @param cv_uf Coefficient of variation for uncertainty factors.
#' @param rmore_method Distribution method for `rmore` ("step" or "lognormal").
#' @param quantile_type Quantile type used when summarizing outputs.
#' @param debug_option Logical; saves intermediate objects for debugging.
#' @param parallel Logical; run combinations in parallel.
#' @param workers Number of workers when `parallel = TRUE`.
#' @param base_cache_dir Directory for cached PSSD objects.
#' @param base_output_path Directory for figures.
#' @param overwrite_cache Logical; recompute even if cache exists.
#' @param progress Logical; show progress bars. If NULL, disables in RStudio or
#'   knitr.
#'
#' @return A named list of PSSD results (one per combination).
#' @export
make_all_pSSDs <- function(
  MC_sim_df = NULL, # MC-sim dataframe
  tiers = c(3), #mehinto et al tiers
  environments = c("Freshwater", "Marine"), #environments to loop through
  erms = c("Food Dilution", "Tissue Translocation"), #ERMs to loop through
  erm_registry = NULL,
  sim = 10,
  cv_uf = 0.5,
  rmore_method = "step",
  quantile_type = 8,
  debug_option = FALSE,
  parallel = TRUE,
  workers = max(1, parallel::detectCores() - 1),
  base_cache_dir = file.path(tempdir(), "pssd_cache"),
  base_output_path = file.path(tempdir(), "pssd_figures"),
  overwrite_cache = FALSE,
  progress = NULL
) {
  # ---- Argument validation ----
  if (is.null(MC_sim_df)) {
    stop("MC_sim_df must be provided explicitly.", call. = FALSE)
  }
  if (is.null(erm_registry)) {
    stop("erm_registry must be provided explicitly.", call. = FALSE)
  }

  progress_enabled <- progress
  if (is.null(progress_enabled)) {
    progress_enabled <- !isTRUE(getOption("knitr.in.progress")) &&
      !identical(Sys.getenv("RSTUDIO"), "1")
  }
  progress_enabled <- isTRUE(progress_enabled) &&
    isTRUE(getOption("progressr.enable", TRUE))
  if (!progress_enabled) {
    old_progressr_opts <- options(progressr.enable = FALSE)
    on.exit(options(old_progressr_opts), add = TRUE)
  }
  # Create output directories if they don't exist
  output_path <- file.path(
    base_output_path,
    paste0(rmore_method, "_", sim, "sims")
  )
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  cache_dir <- file.path(base_cache_dir, paste0(rmore_method, "_", sim, "sims"))
  # Create cache directory if it doesn't exist
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = T)
  }
  # generate color palette for unique levels of species/groups
  all_species <- MC_sim_df |>
    dplyr::filter(
      environment %in% environments,
      !Group %in% c("Insect", "Annelida") # not used in this analysis
    ) |>
    dplyr::distinct(Species, Group, environment) |>
    base::droplevels() |>
    dplyr::arrange(environment, Group, Species)

  # Generate a consistent color palette
  global_color_palette <- generate_color_palette(all_species)

  # Initialize tracking
  total_iterations <- length(tiers) * length(environments) * length(erms)
  iteration_times <- numeric(total_iterations)
  iteration_count <- 0
  tictoc::tic("PSSD++ For Loop Begins...")
  # Progress bar handlers with progressr package
  if (progress_enabled) {
    tryCatch(
      {
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar") # console-safe
      },
      error = function(e) {
        progress_enabled <<- FALSE
      }
    )
  }

  # Start loop
  results <- list()
  # list combinations to loop through
  combo_tbl <- expand.grid(
    tier = tiers,
    environment = environments,
    erm = erms,
    stringsAsFactors = FALSE
  )
  # calculate length of progress bar
  n_combos <- nrow(combo_tbl)
  total_steps <- n_combos * sim

  # Silence repeated package startup messages inside futures
  old_future_opts <- options(future.startup.messages = FALSE)
  on.exit(options(old_future_opts), add = TRUE)
  # set up parllel plan
  if (parallel) {
    silent = TRUE #removes redundant progress update
    cat(sprintf(
      "⚙️ Running %d PSSD combinations in parallel (%d workers)\n",
      nrow(combo_tbl),
      workers
    ))
    # ensure required packages are available on workers (including magrittr for legacy pipes)
    old_future_pkgs <- getOption("future.packages")
    options(
      future.packages = unique(c(
        "PSSDplusplus",
        "magrittr",
        "dplyr",
        "tidyr",
        "ggplot2",
        "reshape2",
        "Matrix",
        "stats",
        "utils",
        "purrr"
      ))
    )
    on.exit(options(future.packages = old_future_pkgs), add = TRUE)
    future::plan(future::multisession, workers = workers)
  } else {
    silent = FALSE # sequential progress update
    cat("⚙️ Running PSSD combinations sequentially\n")
    future::plan(future::sequential)
  }

  on.exit(future::plan(future::sequential), add = TRUE)
  progressr::with_progress({
    p <- progressr::progressor(steps = total_steps)

    results_list <- future.apply::future_lapply(
      seq_len(nrow(combo_tbl)),
      function(i) {
        tier <- combo_tbl$tier[i]
        environment <- combo_tbl$environment[i]
        erm <- combo_tbl$erm[i]

        combo_id <- paste0("Tier", tier, "_", environment, "_", erm)
        cache_file <- file.path(cache_dir, paste0(combo_id, ".rds"))

        # ---- Cache logic ----
        if (file.exists(cache_file) && !overwrite_cache) {
          cat(crayon::blue(sprintf("🔁 Using cached: %s\n", combo_id)))
          return(readRDS(cache_file))
        }

        if (file.exists(cache_file) && overwrite_cache) {
          cat(crayon::yellow(sprintf("♻️ Overwriting cached: %s\n", combo_id)))
        }

        # ---- Run safely ----
        result <- tryCatch(
          {
            results_df <- get_results_df(tier, environment, erm, erm_registry)

            res <- generate_plots_and_summary(
              tier = tier,
              environment = environment,
              erm = erm,
              color_palette = global_color_palette,
              results_df = results_df,
              sim = sim,
              num_iterations = sim,
              quantile_type = quantile_type,
              cv_dp = NULL,
              cv_uf = cv_uf,
              rmore_method = rmore_method,
              species_data_source = MC_sim_df,
              output_path = output_path,
              presentation_path = NULL,
              debug = debug_option,
              silent = silent,
              progress = p # pass progress callback
            )

            # Atomic cache write
            tmp <- paste0(cache_file, ".tmp")
            saveRDS(res, tmp)
            file.rename(tmp, cache_file)

            cat(crayon::yellow(sprintf("✅ Completed: %s\n", combo_id)))
            res
          },
          error = function(e) {
            cat(crayon::red(sprintf(
              "❌ ERROR in %s: %s\n",
              combo_id,
              e$message
            )))
            NULL
          }
        )

        result
      },
      future.seed = TRUE,
      future.packages = c("PSSDplusplus")
    )
  })

  names(results_list) <- paste0(
    "Tier",
    combo_tbl$tier,
    "_",
    combo_tbl$environment,
    "_",
    combo_tbl$erm
  )

  return(results_list)

  # Wrap up
  pSSD_time <- tictoc::toc()
  total_runtime_sec <- pSSD_time$toc - pSSD_time$tic
  cat(crayon::blue(sprintf(
    "\n🎉 PSSD++ complete for all combinations! Total time: %.2f hours\n",
    total_runtime_sec / 3600
  )))
}

#### PNEC Summarization ####
#' Summarize PSSD-derived PNECs
#'
#' Extracts HC5/HC10 summaries from the PSSD result list produced by
#' \link[=make_all_pSSDs]{make_all_pSSDs()} and returns a tidy table.
#'
#' @param pSSDs Named list returned by \link[=make_all_pSSDs]{make_all_pSSDs()}.
#'
#' @return Data frame of summary statistics with combination, environment, ERM,
#'   and HCX labels.
#' @export
summarize_PNECs <- function(pSSDs = NULL) {
  # Initialize an empty list to store the standardized summaries
  PNEC_summaries <- list()

  # Loop through the results list to extract and standardize PNEC summaries
  for (combination in names(pSSDs)) {
    # Extract the PNEC summary statistics for HC5 (hcx = 0.05)
    PNEC_stats_05 <- pSSDs[[combination]]$summary_05$stats
    PNEC_stats_05_long <- data.frame(
      Statistic = rownames(PNEC_stats_05),
      Value = as.numeric(PNEC_stats_05[, 1]), # Extract the first column (values)
      HCX = "HC5", # Add a column to indicate HC5
      Tier = "Tier3" # Assign Tier3 for HC5
    )

    # Extract the PNEC summary statistics for HC10 (hcx = 0.10)
    PNEC_stats_10 <- pSSDs[[combination]]$summary_10$stats
    PNEC_stats_10_long <- data.frame(
      Statistic = rownames(PNEC_stats_10),
      Value = as.numeric(PNEC_stats_10[, 1]), # Extract the first column (values)
      HCX = "HC10", # Add a column to indicate HC10
      Tier = "Tier4" # Assign Tier4 for HC10
    )

    # Combine HC5 and HC10 into a single data frame
    PNEC_stats_long <- rbind(PNEC_stats_05_long, PNEC_stats_10_long)

    # Add the combination name and split it into environment and ERM
    PNEC_stats_long$Combination <- combination
    PNEC_stats_long$Environment <- sub(".*_(.*)_(.*)", "\\1", combination) # Extract environment
    PNEC_stats_long$ERM <- sub(".*_(.*)", "\\1", combination) # Extract ERM

    # Append to the list
    PNEC_summaries[[combination]] <- PNEC_stats_long
  }

  # Combine all standardized summaries into a single data frame
  PNEC_summary_table <- do.call(rbind, PNEC_summaries)

  # Reorder columns for clarity
  PNEC_summary_table <- PNEC_summary_table[, c(
    "Tier",
    "Environment",
    "ERM",
    "HCX",
    "Statistic",
    "Value"
  )]

  # Sort the table by Environment, ERM, Tier, and HCX
  PNEC_summary_table <- PNEC_summary_table[
    order(
      PNEC_summary_table$Environment,
      PNEC_summary_table$ERM,
      PNEC_summary_table$Tier,
      PNEC_summary_table$HCX
    ),
  ]

  # Remove row names
  rownames(PNEC_summary_table) <- NULL

  # Pivot the table wider
  PNEC_summary_table_wide <- PNEC_summary_table %>%
    tidyr::pivot_wider(
      names_from = c(Statistic), # Use both Statistic and HCX for new column names
      values_from = Value # Use the Value column for the values
    ) %>%
    dplyr::select(-HCX)

  return(PNEC_summary_table_wide)
}
