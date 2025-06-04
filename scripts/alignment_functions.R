###########################################
###########################################
############ ALIGNMENT FUNCTION ###########
###########################################
###########################################
# This script provides all necessary equations and functions to perform the food
# dilution and tissue translocation alignments ###

######## Functions ######
#function to derive correction factor (CF) from Koelmans et al (equation 2)
CFfnx = function(a, #default alpha from Koelmans et al (2020)
                 x2D, #set detault values to convert ranges to (1-5,000 um) #5mm is upper defuault 
                 x1D, #1 um is lower default size
                 x2M, x1M){
  
  CF = (x2D^(1-a)-x1D^(1-a))/(x2M^(1-a)-x1M^(1-a))
  
  return(CF)
}

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
        mux.poly[i] <- NA  # Invalid input for log
      }
    } else if (a.x[i] == 2) {
      # Special case: a.x == 2
      epsilon <- 1e-10  # Small value to avoid division by zero
      if (x_UL[i] > 0 && x_LL[i] > 0) {
        mux.poly[i] <- log(x_UL[i] / x_LL[i]) /
          ((x_LL[i] + epsilon)^-1 - (x_UL[i] + epsilon)^-1)
      } else {
        mux.poly[i] <- NA  # Invalid input for log
      }
    } else {
      # General case: a.x != 1 and a.x != 2
      if (x_UL[i] > 0 && x_LL[i] > 0) {
        mux.poly[i] <- ((1 - a.x[i]) / (2 - a.x[i])) *
          ((x_UL[i]^(2 - a.x[i]) - x_LL[i]^(2 - a.x[i])) /
             (x_UL[i]^(1 - a.x[i]) - x_LL[i]^(1 - a.x[i])))
      } else {
        mux.poly[i] <- NA  # Invalid input for power calculations
      }
    }
  }
  
  # Return the result
  return(mux.poly)
}


############# VOLUME ############
volumefnx <- function(R = NA, # average length-to-width ratio for environment
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
SAfnx = function(length,
                 width = NA, 
                 height = NA,
                 R = NA,
                 H_W_ratio = 0.67# assumed 0.67 * width per Kooi et al. (2021)
) {
  # If width unknown, use L:W ratio
  width <- ifelse(is.na(width), R * length, width)
  
  # If height unknown, use H:R ratio
  height <- ifelse(is.na(height), H_W_ratio * width, height)
  
  # a, b, and c are equivalent to 1/2th of the length, width, and height, respectively
  a <- 0.5 * length
  b <- 0.5 * width
  c <- 0.5 * height
  
  SA = (4 * pi) * ((((a*b)^1.6 + (a*c)^1.6 + (b*c)^1.6) / 3) ^ (1/1.6))
  return(SA)}

################# MASS ####################
massfnx <- function(v, p) {
  # If either v or p is NA, return NA for those elements
  mass <- ifelse(is.na(v) | is.na(p), NA, p * v * (1 / 1e12) * 1e6) # correction factor (g to µg)
  return(mass)
}

###### SSA #####
SSA.inversefnx = function(sa, # average surface area
                          m){ #average mass
  SSA.inverse = m/sa
  return(SSA.inverse)}


######################################################
############### MASTER ALIGNMENT FUNCTION ###############
########################################################

preparation_fxn <- function(df, #dataframe input
                          beta_log10_body_length =  0.9341, # Jâms, et al 2020 Nature paper
                          body_length_intercept = 1.1200, # Jâms, et al 2020 Nature paper
                          H_W_ratio = 0.67, # H:W ratio across all environments/particle types
                          R.ave.water.marine = 0.77, # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
                          R.ave.water.freshwater = 0.67, # average length to width ratio of microplastics in freshwater environment (Kooi et al. 2021)
                          R.ave.sediment.marine = 0.75, # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
                          R.ave.sediment.freshwater = 0.70, # average length to width ratio of microplastics in freshwater environment (Kooi et al. 2021)
                          p.ave.marine = 1.10, #average density in marine surface water
                          alpha.marine = 2.07, #table s4 for marine surface water. length
                          a.sa.marine = 1.50, #marine surface area power law
                          a.v.marine = 1.48, #a_V for marine surface water volume
                          a.m.marine = 1.32, # upper limit fora_m for mass for marine surface water in table S4 
                          a.ssa.marine = 1.98, # A_SSA for marine surface water
                          p.ave.freshwater = 1.04, #average density in freshwater surface water
                          alpha.freshwater = 2.64, #table s4 for freshwater surface water. length
                          a.sa.freshwater = 2.00, #freshwater surface area power law
                          a.v.freshwater = 1.68, #a_V for freshwater surface water volume
                          a.m.freshwater = 1.65, # upper limit fora_m for mass for freshwater surface water in table S4 
                          a.ssa.freshwater = 2.71 # A_SSA for freshwater surface water
                          ){
  
  # Check if columns exist and conditionally create them
    if (!"dose.mg.kg.sed.measured" %in% names(df)) {
      df <- df %>%
        mutate(dose.mg.kg.sed.measured = measured.dose.mg.kg.sediment)
    }
  
  if (!"dose.mg.kg.sed.nominal" %in% names(df)) {
    df <- df %>%
      mutate(dose.mg.kg.sed.nominal = nominal.dose.mg.kg.sediment)
  }
  
  if (!"dose.particles.kg.sed.nominal" %in% names(df)) {
    df <- df %>%
      mutate(dose.particles.kg.sed.nominal = nominal.dose.particles.kg.sediment)
  }
  
  if (!"environment" %in% names(df)){
    df <- df %>% 
      mutate(environment = env_f)
    }
  
  df_prepared <- df %>% 
    mutate(beta_log10_body_length = !!beta_log10_body_length,
           body_length_intercept = !!body_length_intercept,
           H_W_ratio = !!H_W_ratio,
           R.ave.water.marine = !!R.ave.water.marine,                                 
           R.ave.water.freshwater = !!R.ave.water.freshwater, 
           R.ave.sediment.marine = !!R.ave.sediment.marine,
           R.ave.sediment.freshwater = !!R.ave.sediment.freshwater, 
           p.ave.marine = !!p.ave.marine,
           alpha.marine = !!alpha.marine,
           a.sa.marine = !!a.sa.marine, 
           a.v.marine = !!a.v.marine, 
           a.m.marine = !!a.m.marine, 
           a.ssa.marine = !!a.ssa.marine, 
           p.ave.freshwater = !!p.ave.freshwater, 
           alpha.freshwater = !!alpha.freshwater,
           a.sa.freshwater = !!a.sa.freshwater, 
           a.v.freshwater = !!a.v.freshwater, 
           a.m.freshwater = !!a.m.freshwater,                     
           a.ssa.freshwater = !!a.ssa.freshwater) %>% 
    
    ##### define environment-specific alpha parameters ############### 
    mutate(alpha = case_when(environment == "Marine"  & exposure.route == "water" ~ alpha.marine,
                             environment == "Freshwater" & exposure.route == "water" ~ alpha.freshwater),
           a.sa = case_when(environment == "Marine" & exposure.route == "water" ~ a.sa.marine,
                            environment == "Freshwater" & exposure.route == "water" ~ a.sa.freshwater),
           a.v = case_when(environment == "Marine" & exposure.route == "water" ~ a.v.marine,
                           environment == "Freshwater" & exposure.route == "water" ~ a.v.freshwater),
           a.m = case_when(environment == "Marine" & exposure.route == "water" ~ a.m.marine,
                           environment == "Freshwater" & exposure.route == "water" ~ a.m.freshwater),
           a.ssa = case_when(environment == "Marine" & exposure.route == "water" ~ a.ssa.marine,
                             environment == "Freshwater" & exposure.route == "water" ~ a.ssa.freshwater),
           R.ave = case_when(environment == "Marine" & exposure.route == "water" ~ R.ave.water.marine,
                             environment == "Freshwater" & exposure.route == "water" ~ R.ave.water.freshwater,
                             environment == "Marine" & exposure.route == "sediment" ~ R.ave.sediment.marine,
                             environment == "Freshwater" & exposure.route == "sediment" ~ R.ave.sediment.freshwater),
           p.ave = case_when(environment == "Marine" & exposure.route == "water" ~ p.ave.marine,
                             environment == "Freshwater" & exposure.route == "water" ~ p.ave.freshwater,
                             environment == "Marine" & exposure.route == "sediment" ~ p.ave.marine,
                             environment == "Freshwater" & exposure.route == "sediment" ~ p.ave.freshwater),
           H_W_ratio = H_W_ratio # currently just using one value for all environments
    ) %>% 
    mutate(shape = shape_f) %>% 
    ## calculate size parameters using compartment characteristics
    mutate(size.width.min.um.used.for.conversions = case_when(
      shape == "sphere" ~ size.length.min.um.used.for.conversions, #all dims same
      shape == "fiber" ~ R.ave * size.length.min.um.used.for.conversions, #median holds for all particles (Kooi et al 2021)
      shape == "Not Reported" ~ R.ave * size.length.min.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
      shape == "fragment" ~ R.ave * size.length.min.um.used.for.conversions)) %>% # average width to length ratio in the marine environment (kooi et al 2021)
    mutate(size.height.min.um.used.for.conversions = case_when(
      shape == "sphere" ~ size.length.min.um.used.for.conversions, #all dims same
      shape == "Not Reported" ~ R.ave * H_W_ratio * size.length.min.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
      shape == "fiber" ~  R.ave * size.length.min.um.used.for.conversions, #height same as width for fibers
      shape == "fragment" ~ R.ave * H_W_ratio * size.length.min.um.used.for.conversions)) %>% # average width to length ratio in the marine environment AND average height to width ratio (kooi et al 2021)
    # maxima
    mutate(size.length.max.um.used.for.conversions = case_when(
      is.na(size.length.max.mm.measured) ~ size.length.max.mm.nominal * 1000,
      !is.na(size.length.max.mm.measured) ~ size.length.max.mm.measured * 1000)) %>% 
    mutate(size.width.max.um.used.for.conversions = case_when(
      shape == "sphere" ~ size.length.max.um.used.for.conversions, #all dims same
      shape == "fiber" ~ R.ave * size.length.max.um.used.for.conversions, #median holds for all particles (Kooi et al 2021) #there are no fibers
      shape == "Not Reported" ~ R.ave * size.length.max.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
      shape == "fragment" ~ R.ave * size.length.max.um.used.for.conversions)) %>% # average width to length ratio in the marine environment (kooi et al 2021)
    mutate(size.height.max.um.used.for.conversions = case_when(
      shape == "sphere" ~ size.length.max.um.used.for.conversions, #all dims same
      shape == "Not Reported" ~ R.ave * H_W_ratio * size.length.max.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
      shape == "fiber" ~ R.ave * size.length.max.um.used.for.conversions, #hieght same as width
      shape == "fragment" ~ R.ave * H_W_ratio * size.length.max.um.used.for.conversions)) %>%  # average width to length ratio in the marine environment AND average height to width ratio (kooi et al 2021)
    # first ensure that width and height are filled out (sometimes there are reported, but if not, some default assumptions are made to estimate them)
    mutate(size.width.um.used.for.conversions = case_when(
      is.na(size.width.um.used.for.conversions) & shape_f == "Fiber" ~ 15, # assume 15 um width for fibers unless already known (kooi et al. 2021)
      is.na(size.width.um.used.for.conversions) & shape_f == "Sphere" ~ size.length.um.used.for.conversions, # W = L for spheres
      is.na(size.width.um.used.for.conversions) & shape_f == "Fragment" ~ size.length.um.used.for.conversions * R.ave, #use average width:length ratio for fragments
      T ~ size.width.um.used.for.conversions # if available, use as-is
    )) %>% 
    #estimate height based on shape (data doesn't exist in ToMEx for monodisperse, because never reported)
    mutate(size.height.um.used.for.conversions = case_when(
      shape_f == "Sphere" ~ size.length.um.used.for.conversions, # if spherical, height = length
      shape_f != "Sphere" ~ size.width.um.used.for.conversions * H_W_ratio # if not spherical, height = width * H:W ratio
    )) %>% 
    #### use assessment factors to get chronic NOECs ##
    #mutate(dose.particles.mL.AF.corrected = dose.particles.mL.master / (af.time * af.noec)) %>%
    # if NA, then there's missing AFs - data are considered invalid #
    #drop_na(dose.particles.mL.AF.corrected) %>%
    # rename
    #mutate(dose.particles.mL.master = dose.particles.mL.AF.corrected) %>% 
    # calculate volume for monodisperse particles #
    mutate(particle.volume.um3 = volumefnx(R = R.ave,
                                           length = size.length.um.used.for.conversions, 
                                           width = size.width.um.used.for.conversions,
                                           height = size.height.um.used.for.conversions
    )) %>% 
    # calculate min and max volume when polydisperse particles are used (being sure to use ingestion-restricted sizes)
    mutate(particle.volume.um3.min = volumefnx(R = R.ave, 
                                               length = size.length.min.um.used.for.conversions,
                                               width = size.width.min.um.used.for.conversions, 
                                               height = size.height.min.um.used.for.conversions),
           particle.volume.um3.max = volumefnx(R = R.ave,
                                               length = size.length.max.um.used.for.conversions,
                                               width = size.width.max.um.used.for.conversions, 
                                               height = size.height.max.um.used.for.conversions)) %>% 
    # calculate surface are for monodisperse particles
    mutate(particle.surface.area.um2 = SAfnx(length = size.length.um.used.for.conversions,
                                             width = size.width.um.used.for.conversions,
                                             height = size.height.um.used.for.conversions,
                                             R = R.ave,
                                             H_W_ratio = H_W_ratio)) %>% 
    # calculate min/max SA for polydisperse mixtures (being sure to use translocation-restricted polydisperse upper sizes)
    mutate(particle.surface.area.um2.min = SAfnx(length = size.length.min.um.used.for.conversions,
                                                 width = size.width.min.um.used.for.conversions,
                                                 height = size.height.min.um.used.for.conversions,
                                                 R = R.ave,
                                                 H_W_ratio = H_W_ratio),
           particle.surface.area.um2.max = SAfnx(length = size.length.max.um.used.for.conversions,
                                                 width = size.width.max.um.used.for.conversions,
                                                 height = size.height.max.um.used.for.conversions,
                                                 R = R.ave,
                                                 H_W_ratio = H_W_ratio)) %>% 
    #calculate minimum and maximum mass for polydisperse particles
    mutate(mass.per.particle.mg.min = massfnx(v = particle.volume.um3.min, p = density.g.cm3) * 1e-3) %>% #equation uses g/cm3
    mutate(mass.per.particle.mg.max = massfnx(v = particle.volume.um3.max, p = density.g.cm3) * 1e-3) %>%   #equation uses g/cm3
    mutate(mass.per.particle.mg = massfnx(v = particle.volume.um3, p = density.g.cm3) * 1e-3) %>%   #equation uses g/cm3
    ########## DOSE METRICS ################

    #calcualte dose metrics accordingly
    mutate(dose.surface.area.um2.mL.master = particle.surface.area.um2 * dose.particles.mL.master) %>% 
    mutate(particle.surface.area.um2.mg = particle.surface.area.um2 / mass.per.particle.mg) %>% 
    
    #Sediment-based concentration metrics
    mutate(dose.mg.kg.sediment.master = if_else(!is.na(dose.mg.kg.sed.measured), dose.mg.kg.sed.measured, dose.mg.kg.sed.nominal)) %>% #Create master column with measured concentrations preferred
    mutate(dose.particles.kg.sediment.master = dose.particles.kg.sed.nominal) %>% #Create master column with measured concentrations preferred (only nominal concentrations available)
    
    #Create reported vs. converted columns for sediment-based metrics
    mutate(dose.mg.kg.sediment.master.converted.reported = if_else(!is.na(dose.mg.kg.sediment.master), "reported", NA_character_)) %>% 
    mutate(dose.particles.kg.sediment.master.converted.reported = if_else(!is.na(dose.particles.kg.sediment.master), "reported", NA_character_)) %>%  
    
    #Sediment Mass (converted)
    mutate(dose.mg.kg.sediment.master = ifelse(is.na(dose.mg.kg.sediment.master), (dose.particles.kg.sediment.master)*mass.per.particle.mg, dose.mg.kg.sediment.master)) %>% 
    mutate(dose.mg.kg.sediment.master.converted.reported = factor(ifelse((!is.na(dose.mg.kg.sediment.master)&is.na(dose.mg.kg.sediment.master.converted.reported)), "converted", dose.mg.kg.sediment.master.converted.reported))) %>% 
    
    #Sediment Count (converted)
    mutate(dose.particles.kg.sediment.master = ifelse(is.na(dose.particles.kg.sediment.master), (dose.mg.kg.sediment.master)/mass.per.particle.mg, dose.particles.kg.sediment.master)) %>% 
    mutate(dose.particles.kg.sediment.master.converted.reported = factor(ifelse((!is.na(dose.particles.kg.sediment.master)&is.na(dose.particles.kg.sediment.master.converted.reported)), "converted", dose.particles.kg.sediment.master.converted.reported))) %>%  
    
    #Volume
    mutate(dose.um3.mL.master = particle.volume.um3 * dose.particles.mL.master) %>%  #calculate volume/mL
    mutate(dose.um3.kg.sediment.master = particle.volume.um3 * dose.particles.kg.sediment.master) %>% #calculate volume/kg sediment
    
    #Surface Area
    mutate(dose.um2.mL.master = as.numeric(particle.surface.area.um2) * dose.particles.mL.master) %>% 
    mutate(dose.um2.kg.sediment.master = as.numeric(particle.surface.area.um2) * dose.particles.kg.sediment.master) %>% 
    
    #Specific Surface Area
    mutate(dose.um2.ug.mL.master = dose.um2.mL.master / (mass.per.particle.mg / 1000)) %>% #correct mg to ug
    mutate(dose.um2.ug.kg.sediment.master = dose.um2.kg.sediment.master/(mass.per.particle.mg / 1000)) %>% 
    ################## BIOAVAILABILITY ##############
    #### Estimate ingestible plastic size
    mutate(max.size.ingest.um = 1000 * (10^(beta_log10_body_length * log10(body.length.cm * 10) - body_length_intercept)))#(Jâms, et al 2020 Nature paper)correction for cm to mm
    
  
  return(df_prepared)
}

######################### ALIGNMENT FUNCTION #############
alignment_fxn <- function(df_prepared, #dataframe input
                            x1D_set = 1, # lower default distribution (typically 1 um)
                            x2D_set = 5000, # upper default distribution (typically 5000 um)
                            x1M_set = 1, # lower bioavailable size limit (typically 1 um)
                          
                            upper.tissue.trans.size.um = 88, # Coffin et al. (2025)
                            beta_log10_body_length =  0.9341, # Jâms, et al 2020 Nature paper 
                            body_length_intercept = 1.1200, # Jâms, et al 2020 Nature paper
                          
                            H_W_ratio = 0.67, # H:W ratio across all environments/particle types
                          
                            R.ave.water.marine = 0.77, # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
                            R.ave.water.freshwater = 0.67, # average length to width ratio of microplastics in freshwater environment (Kooi et al. 2021)
                            R.ave.sediment.marine = 0.75, # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
                            R.ave.sediment.freshwater = 0.70, # average length to width ratio of microplastics in freshwater environment (Kooi et al. 2021)
                          
                            p.ave.marine = 1.10, #average density in marine surface water
                            alpha.marine = 2.07, #table s4 for marine surface water. length
                            a.sa.marine = 1.50, #marine surface area power law
                            a.v.marine = 1.48, #a_V for marine surface water volume
                            a.m.marine = 1.32, # upper limit fora_m for mass for marine surface water in table S4 
                            a.ssa.marine = 1.98, # A_SSA for marine surface water
                          
                            p.ave.freshwater = 1.04, #average density in freshwater surface water
                            alpha.freshwater = 2.64, #table s4 for freshwater surface water. length
                            a.sa.freshwater = 2.00, #freshwater surface area power law
                            a.v.freshwater = 1.68, #a_V for freshwater surface water volume
                            a.m.freshwater = 1.65, # upper limit fora_m for mass for freshwater surface water in table S4 
                            a.ssa.freshwater = 2.71 # A_SSA for freshwater surface water
){
  
  ### Create EC_mono_p.particles.mL only if it doesn't already exist (ensures compatibility with Shiny app) ##
  if (!"EC_mono_p.particles.mL" %in% names(df_prepared)) {
    df_prepared <- df_prepared %>% 
      mutate(EC_mono_p.particles.mL = dose.particles.mL.master)
  }
  
  # Check if columns exist and conditionally create them
  if (!"dose.mg.kg.sed.measured" %in% names(df_prepared)) {
    df_prepared <- df_prepared %>%
      mutate(dose.mg.kg.sed.measured = measured.dose.mg.kg.sediment)
  }
  
  if (!"dose.mg.kg.sed.nominal" %in% names(df_prepared)) {
    df_prepared <- df_prepared %>%
      mutate(dose.mg.kg.sed.nominal = nominal.dose.mg.kg.sediment)
  }
  
  if (!"dose.particles.kg.sed.nominal" %in% names(df_prepared)) {
    df_prepared <- df_prepared %>%
      mutate(dose.particles.kg.sed.nominal = nominal.dose.particles.kg.sediment)
  }
  
  if (!"environment" %in% names(df)){
    df_prepared <- df_prepared %>% 
      mutate(environment = env_f)
  }
  
  df_aligned <- df_prepared %>% 
    mutate(x1D_set = !!x1D_set,
           x2D_set = !!x2D_set,
           x1M_set = !!x1M_set,
           upper.tissue.trans.size.um = !!upper.tissue.trans.size.um,
           beta_log10_body_length = !!beta_log10_body_length,
           body_length_intercept = !!body_length_intercept,
           H_W_ratio = !!H_W_ratio,
           R.ave.water.marine = !!R.ave.water.marine,                                 
           R.ave.water.freshwater = !!R.ave.water.freshwater, 
           R.ave.sediment.marine = !!R.ave.sediment.marine,
           R.ave.sediment.freshwater = !!R.ave.sediment.freshwater, 
           p.ave.marine = !!p.ave.marine,
           alpha.marine = !!alpha.marine,
           a.sa.marine = !!a.sa.marine, 
           a.v.marine = !!a.v.marine, 
           a.m.marine = !!a.m.marine, 
           a.ssa.marine = !!a.ssa.marine, 
           p.ave.freshwater = !!p.ave.freshwater, 
           alpha.freshwater = !!alpha.freshwater,
           a.sa.freshwater = !!a.sa.freshwater, 
           a.v.freshwater = !!a.v.freshwater, 
           a.m.freshwater = !!a.m.freshwater,                     
           a.ssa.freshwater = !!a.ssa.freshwater) %>% 
    ## assign alpha values
    mutate(environment = env_f) %>% 
    mutate(alpha = case_when(environment == "Marine"  & exposure.route == "water" ~ alpha.marine,
                             environment == "Freshwater" & exposure.route == "water" ~ alpha.freshwater),
           a.sa = case_when(environment == "Marine" & exposure.route == "water" ~ a.sa.marine,
                            environment == "Freshwater" & exposure.route == "water" ~ a.sa.freshwater),
           a.v = case_when(environment == "Marine" & exposure.route == "water" ~ a.v.marine,
                           environment == "Freshwater" & exposure.route == "water" ~ a.v.freshwater),
           a.m = case_when(environment == "Marine" & exposure.route == "water" ~ a.m.marine,
                           environment == "Freshwater" & exposure.route == "water" ~ a.m.freshwater),
           a.ssa = case_when(environment == "Marine" & exposure.route == "water" ~ a.ssa.marine,
                             environment == "Freshwater" & exposure.route == "water" ~ a.ssa.freshwater),
           R.ave = case_when(environment == "Marine" & exposure.route == "water" ~ R.ave.water.marine,
                             environment == "Freshwater" & exposure.route == "water" ~ R.ave.water.freshwater,
                             environment == "Marine" & exposure.route == "sediment" ~ R.ave.sediment.marine,
                             environment == "Freshwater" & exposure.route == "sediment" ~ R.ave.sediment.freshwater),
           p.ave = case_when(environment == "Marine" & exposure.route == "water" ~ p.ave.marine,
                             environment == "Freshwater" & exposure.route == "water" ~ p.ave.freshwater,
                             environment == "Marine" & exposure.route == "sediment" ~ p.ave.marine,
                             environment == "Freshwater" & exposure.route == "sediment" ~ p.ave.freshwater),
           H_W_ratio = H_W_ratio # currently just using one value for all environments
    ) %>% 
      ##### filter out undesired data ####
  # ensure no nanoparticle (or simply particles outside desired range) studies are used
  filter(size.length.um.used.for.conversions >= x1D_set,
         size.length.um.used.for.conversions <= x2D_set) %>% 
    mutate(nanoparticle_polydisperse = case_when(
      polydispersity == "polydisperse" & size.length.min.um.used.for.conversions < x1D_set ~ "nanoparticle_filter out",
      T ~ "don't filter out")) %>% 
    filter(nanoparticle_polydisperse == "don't filter out") %>% 
    ############################
  ######### BIOAVAILABILITY##########
  #################################
  mutate(max.size.ingest.um = 1000 * (10^(beta_log10_body_length * log10(body.length.cm * 10) - body_length_intercept))) %>% #(Jâms, et al 2020 Nature paper)correction for cm to mm
    # max ingestible particle size for bioavailability restriction
    mutate(x2M_ingest = case_when(max.size.ingest.um < x2D_set ~ max.size.ingest.um, 
                                  max.size.ingest.um >= x2D_set ~ x2D_set), #if max ingestible size is bigger than default distribution size, then just use that
           # max translocatable particle size for bioavailability restriction
           x2M_trans = case_when(max.size.ingest.um < upper.tissue.trans.size.um ~  max.size.ingest.um,
                                 max.size.ingest.um >= upper.tissue.trans.size.um ~ upper.tissue.trans.size.um)) %>%  
    # tag whether monodisperse particles are ingestible/translocatable for each species #
    mutate(ingestible = case_when(
      polydispersity == "monodisperse" & size.length.um.used.for.conversions <= x2M_ingest ~ "ingestible",
      polydispersity == "monodisperse" & size.length.um.used.for.conversions > x2M_ingest ~ "not ingestible"),
      translocatable = case_when(
        polydispersity == "monodisperse" & size.length.um.used.for.conversions <= x2M_trans ~ "translocatable",
        polydispersity == "monodisperse" & size.length.um.used.for.conversions > x2M_trans ~ "not translocatable"),
      #### tag whether polydipserse particles are ingestible/translocatable
      ingestible_poly = case_when(
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions <= x2M_ingest & size.length.min.um.used.for.conversions <= x2M_ingest ~ "ingestible (all)",
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions > x2M_ingest & size.length.min.um.used.for.conversions <= x2M_ingest ~ "ingestible (some)",
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions > x2M_ingest & size.length.min.um.used.for.conversions > x2M_ingest ~ "not ingestible"),
      translocatable_poly = case_when(
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions <= x2M_trans & size.length.min.um.used.for.conversions <= x2M_trans ~ "translocatable (all)",
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions > x2M_trans & size.length.min.um.used.for.conversions <= x2M_trans ~ "translocatable (some)",
        polydispersity == "polydisperse" & size.length.max.um.used.for.conversions > x2M_trans & size.length.min.um.used.for.conversions > x2M_trans ~ "not translocatable")
    ) %>%  
    # filter out studies that are not translocatable (as they will not be used in any alignments)
    filter(!grepl("not", translocatable_poly)) %>% 
    filter(!grepl("not", translocatable)) %>% 
    # correct for partially translocatable particles  by calculating translocatable fraction
    mutate(CF_bioavailable_trans = case_when(
      translocatable_poly == "translocatable (some)" ~ CFfnx(a = alpha,
                                                             x1D = size.length.min.um.used.for.conversions,
                                                             x2D = x2M_trans,
                                                             x1M = size.length.min.um.used.for.conversions,
                                                             x2M = size.length.max.um.used.for.conversions),
      T ~ 1)) %>% 
    # calculate translocatable dose using correction factor (fraction)
    mutate(EC_mono_p.particles.mL_trans = case_when(
      translocatable_poly == "translocatable (some)" ~ CF_bioavailable_trans * EC_mono_p.particles.mL,
      T ~ EC_mono_p.particles.mL)
      ) %>% 
    # determine fraction of partially ingestible particles using correction factor
    mutate(CF_bioavailable_ingest = case_when(
      ingestible_poly == "ingestible (some)" ~ CFfnx(a = alpha,
                                                     x1D = size.length.min.um.used.for.conversions,
                                                     x2D = x2M_ingest,
                                                     x1M = size.length.min.um.used.for.conversions,
                                                     x2M = size.length.max.um.used.for.conversions),
      T ~ 1 )) %>% 
    # calculate ingestible particle effect concentration using correction factoor (fraction)
    mutate(EC_mono_p.particles.mL_ingest = case_when(
      ingestible_poly == "ingestible (some)" ~ CF_bioavailable_ingest * EC_mono_p.particles.mL,
      T ~ dose.particles.mL.master)) %>% 
    ## restrict size range of partially bioavailable polydisperse mixtures to bioavailable fractions (just the max sizes in each dimension)
    mutate(size.length.max.um.trans = case_when(
      translocatable_poly == "translocatable (some)" ~ x2M_trans,
      T ~ size.length.max.um.used.for.conversions),
      # adjust width based on existing L:W ratio, and fraction of length that's used
      size.width.max.um.trans = case_when(
        translocatable_poly == "translocatable (some)" ~ size.width.max.um.used.for.conversions * (x2M_trans / size.length.max.um.trans),
        T ~ size.width.max.um.used.for.conversions),
      # adjust height based on existing H:L ratio, and fraction of length that's used
      size.height.max.um.trans = case_when(
        translocatable_poly == "translocatable (some)" ~ size.height.max.um.used.for.conversions * (x2M_trans / size.length.max.um.trans),
        T ~ size.height.max.um.used.for.conversions),
      #### perform same operations for partially ingestible particles
      size.length.max.um.ingest = case_when(
        ingestible_poly == "ingestible (some)" ~ x2M_ingest,
        T ~ size.length.max.um.used.for.conversions),
      # adjust width based on existing L:W ratio, and fraction of length that's used
      size.width.max.um.ingest = case_when(
        ingestible_poly == "ingestible (some)" ~ size.width.max.um.used.for.conversions * (x2M_ingest / size.length.max.um.ingest),
        T ~ size.width.max.um.used.for.conversions),
      # adjust height based on existing H:L ratio, and fraction of length that's used
      size.height.max.um.ingest = case_when(
        ingestible_poly == "ingestible (some)" ~ size.height.max.um.used.for.conversions * (x2M_ingest / size.length.max.um.ingest),
        T ~ size.height.max.um.used.for.conversions),
    ) %>% 
    ###################################################
  ###################### Re-calculate surface area/volume with bioavailable polydisperse fractions
  ############################################################
  # calculate min and max volume when polydisperse particles are used (being sure to use ingestion-restricted sizes)
  mutate(particle.volume.um3.min = volumefnx(R = R.ave, 
                                             length = size.length.min.um.used.for.conversions,
                                             width = size.width.min.um.used.for.conversions, 
                                             height = size.height.min.um.used.for.conversions),
         particle.volume.um3.max = volumefnx(R = R.ave,
                                             length = size.length.max.um.ingest,
                                             width = size.width.max.um.ingest, 
                                             height = size.height.max.um.ingest)) %>% 
    # calculate surface are for monodisperse particles
    mutate(particle.surface.area.um2 = SAfnx(length = size.length.um.used.for.conversions,
                                             width = size.width.um.used.for.conversions,
                                             height = size.height.um.used.for.conversions,
                                             R = R.ave,
                                             H_W_ratio = H_W_ratio)) %>% 
    # calculate min/max SA for polydisperse mixtures (being sure to use translocation-restricted polydisperse upper sizes)
    mutate(particle.surface.area.um2.min = SAfnx(length = size.length.min.um.used.for.conversions,
                                                 width = size.width.min.um.used.for.conversions,
                                                 height = size.height.min.um.used.for.conversions,
                                                 R = R.ave,
                                                 H_W_ratio = H_W_ratio),
           particle.surface.area.um2.max = SAfnx(length = size.length.max.um.trans,
                                                 width = size.width.max.um.trans,
                                                 height = size.height.max.um.trans,
                                                 R = R.ave,
                                                 H_W_ratio = H_W_ratio)) %>% 
  
    #####################################################
    ##################### ALIGNMENTS ####################
  #########################################################
    ##### Determine CF_bio for ERM of interest ###
    # calculate CF_bio for translocation
    mutate(CF_bio_trans = CFfnx(x1M = x1M_set,#lower size bin
                                x2M = x2M_trans, #upper translocatable
                                x1D = x1D_set, #default
                                x2D = x2D_set,  #default
                                a = alpha),
           CF_bio_ingest = CFfnx(x1M = x1M_set,#lower size bin
                                 x2M = x2M_ingest, #upper ingestible length
                                 x1D = x1D_set, #default
                                 x2D = x2D_set,  #default upper size range
                                 a = alpha)) %>% 
    
    ################################################
  ############## Particle ERM #####################
  ################################################
  ### Particle ERM ###
  # calculate effect threshold for particles
    mutate(mu.p.mono = 1) %>% #mu_x_mono is always 1 for particles to particles
    #### translocation-limited ####
    mutate(mu.p.poly_trans = mux_polyfnx(a.x = alpha, x_UL = x2M_trans, x_LL = x1M_set)) %>% 
    # polydisperse effect threshold for particles
    mutate(EC_poly_p.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.p.mono)/mu.p.poly_trans) %>% 
    ## Calculate environmentally relevant effect threshold for particles
    mutate(EC_env_p.particles.mL_trans = EC_poly_p.particles.mL_trans * CF_bio_trans) %>%  #aligned particle effect concentration (1-5000 um)
    
    ##### Ingestion-limited ###
    mutate(mu.p.poly_ingest = mux_polyfnx(a.x = alpha, x_UL = x2M_ingest, x_LL = x1M_set)) %>% 
    # polydisperse effect threshold for particles
    mutate(EC_poly_p.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.p.mono)/mu.p.poly_ingest) %>% 
    ## Calculate environmentally relevant effect threshold for particles
    mutate(EC_env_p.particles.mL_ingest = EC_poly_p.particles.mL_ingest * CF_bio_ingest) %>%  #aligned particle effect concentration (1-5000 um)
    
  ################################################
################### FOOD DILUTION ERM ##########
###############################################
     # calculate volume for monodisperse particles #
    mutate(particle.volume.um3 = volumefnx(R = R.ave,
                                           length = size.length.um.used.for.conversions, 
                                           width = size.width.um.used.for.conversions,
                                           height = size.height.um.used.for.conversions
    )) %>% 
    # calculate min and max volume when polydisperse particles are used (being sure to use ingestion-restricted sizes)
    mutate(particle.volume.um3.min = volumefnx(R = R.ave, 
                                               length = size.length.min.um.used.for.conversions,
                                               width = size.width.min.um.used.for.conversions, 
                                               height = size.height.min.um.used.for.conversions),
           particle.volume.um3.max = volumefnx(R = R.ave,
                                               length = size.length.max.um.ingest,
                                               width = size.width.max.um.ingest, 
                                               height = size.height.max.um.ingest)) %>% 
    # now determine mu.v.mono for monodisperse and polydisperse lab exposure particles
    mutate(mu.v.mono = case_when(
      polydispersity == "monodisperse" ~ particle.volume.um3, # use reported volume in monodisperse
      polydispersity == "polydisperse" ~ mux_polyfnx(a.x = a.v, 
                                                     x_LL = particle.volume.um3.min,
                                                     x_UL = particle.volume.um3.max))) %>% 
    #### INGESTION-LIMITED ####
    #ingestion-limited lower/upper bioavailability limits
    mutate(x_LL_v_ingest = volumefnx(length = x1D_set,
                                     width = x1D_set,
                                     height = x1D_set),
           x_UL_v_ingest = volumefnx(length = x2M_ingest, 
                                     width = x2M_ingest,
                                     height = x2M_ingest)) %>%
    # translate to environmental
    mutate(mu.v.poly_ingest = mux_polyfnx(a.v, x_UL_v_ingest, x_LL_v_ingest),
           EC_poly_v.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.v.mono)/mu.v.poly_ingest,
           EC_env_v.particles.mL_ingest = EC_poly_v.particles.mL_ingest * CF_bio_ingest#,
           # final ingestion-based ERM excludes algae
        #   EC_env_v.particles.mL_ingest =  case_when(
         #    Group == "Algae" ~ NA,
  #           T ~ EC_env_v.particles.mL_ingest)
  ) %>% 
    mutate(particles.mL.food.dilution = EC_env_v.particles.mL_ingest) %>% 
    
    #### TRANSLOCATION_LIMITED###
    mutate(x_LL_v_trans = volumefnx(length = x1D_set,
                                     width = x1D_set,
                                     height = x1D_set),
           x_UL_v_trans = volumefnx(length = x2M_trans, 
                                     width = x2M_trans,
                                     height = x2M_trans)) %>%
    # translate to environmental
    mutate(mu.v.poly_trans = mux_polyfnx(a.v, x_UL_v_trans, x_LL_v_trans),
           EC_poly_v.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.v.mono)/mu.v.poly_trans,
           EC_env_v.particles.mL_trans = EC_poly_v.particles.mL_trans * CF_bio_trans) %>% 
  ###############################################
  ################ TISSUE TRANSLOCATION ERM #######
  ##################################################
 # step 1: ensure particle surface area is estimated for exposure particles
mutate(particle.surface.area.um2 = SAfnx(length = size.length.um.used.for.conversions,
                                         width = size.width.um.used.for.conversions,
                                         height = size.height.um.used.for.conversions,
                                         R = R.ave,
                                         H_W_ratio = H_W_ratio)) %>% 
  # calculate min/max SA for polydisperse mixtures (being sure to use translocation-restricted polydisperse upper sizes)
  mutate(particle.surface.area.um2.min = SAfnx(length = size.length.min.um.used.for.conversions,
                                               width = size.width.min.um.used.for.conversions,
                                               height = size.height.min.um.used.for.conversions,
                                               R = R.ave,
                                               H_W_ratio = H_W_ratio),
         particle.surface.area.um2.max = SAfnx(length = size.length.max.um.trans,
                                               width = size.width.max.um.trans,
                                               height = size.height.max.um.trans,
                                               R = R.ave,
                                               H_W_ratio = H_W_ratio)) %>% 
  # calculate mu.sa.mono for mono and polydisperse particles
  # mu.x_poly equation must be used in case of polydisperse exposure concentrations 
  mutate(mu.sa.mono = case_when(
    polydispersity == "monodisperse" ~ particle.surface.area.um2, # use reported surface area in monodisperse
    polydispersity == "polydisperse" ~  mux_polyfnx(a.x = a.sa, 
                                                    x_LL = particle.surface.area.um2.min,
                                                    x_UL = particle.surface.area.um2.max))) %>% 
  #TRANSLOCATION-LIMTED #
    #calculate lower translocatable surface area using spherical assumption
  mutate(x_LL_sa_trans = SAfnx(length = x1D_set, 
                               width = x1D_set, 
                               height = x1D_set),
  #calculate upper translocatable surface area using spherical assumption
         x_UL_sa_trans = SAfnx(length = x2M_trans, 
                               width = x2M_trans, 
                               height = x2M_trans)) %>%  
  #calculate SA mu_poly and EC_poly_SA (translocatation-limited)
  mutate(mu.sa.poly_trans = mux_polyfnx(a.sa, x_UL_sa_trans, x_LL_sa_trans),
         EC_poly_sa.particles.mL_trans = (EC_mono_p.particles.mL_trans  * mu.sa.mono)/mu.sa.poly_trans,
         # calculate EC_env_sa (translocation limited). This is the Tissue Translocation ERM Value!
         EC_env_sa.particles.mL_trans = EC_poly_sa.particles.mL_trans * CF_bio_trans,
         particles.mL.tissue.translocation = EC_env_sa.particles.mL_trans) %>% 
    #INGESTION-LIMITED
    #calculate lower ingestible surface area using spherical assumption
    mutate(x_LL_sa_ingest = SAfnx(length = x1D_set, 
                                 width = x1D_set, 
                                 height = x1D_set),
           #calculate upper ingestible surface area using spherical assumption
           x_UL_sa_ingest = SAfnx(length = x2M_ingest, 
                                 width = x2M_ingest, 
                                 height = x2M_ingest)) %>%  
    #calculate SA mu_poly and EC_poly_SA (ingestion-limited)
    mutate(mu.sa.poly_ingest = mux_polyfnx(a.sa, x_UL_sa_ingest, x_LL_sa_ingest),
           EC_poly_sa.particles.mL_ingest = (EC_mono_p.particles.mL_ingest  * mu.sa.mono)/mu.sa.poly_ingest,
           # calculate EC_env_sa (ingestion limited). This is the Tissue ingestlocation ERM Value!
           EC_env_sa.particles.mL_ingest = EC_poly_sa.particles.mL_ingest * CF_bio_ingest) %>% 
    
    
    #################################################
    ################ MASS ERM ######################
  ####################################################
  #calculate lower ingestible mass (translocation-limited)
  mutate(x_LL_m_trans = massfnx(v = x_LL_v_trans, p = p.ave),
         x_UL_m_trans = massfnx(v = x_UL_v_trans, p = p.ave)) %>% #average density
    #calculate lower ingestible mass (ingestion-limited)
    mutate(x_LL_m_ingest = massfnx(v = x_LL_v_ingest, p = p.ave),
           x_UL_m_ingest = massfnx(v = x_UL_v_ingest, p = p.ave)) %>% #average density
    # calculate mu.m.poly (trans / ingest)
    mutate(mu.m.poly_trans = mux_polyfnx(a.m, x_UL_m_trans, x_LL_m_trans),
           mu.m.poly_ingest = mux_polyfnx(a.m, x_UL_m_ingest, x_LL_m_ingest)
           ) %>% 
    ##--- laboratory calculations ---###
    ## define mu_x_mono OR mu_x_poly (lab) for alignment to ERM  #
    #(note that if mixed particles were used, a different equation must be used)
    mutate(mu.m.mono = case_when(
      polydispersity == "monodisperse" ~  mass.per.particle.mg * 1000, # use reported volume in monodisperse
      polydispersity == "polydisperse" ~ mux_polyfnx(a.x = a.m, 
                                                     x_UL = mass.per.particle.mg.max * 1000,
                                                     x_LL = mass.per.particle.mg.min * 1000))) %>% 
    #calculate polydisperse effect concentration for volume (particles/mL)
    mutate(EC_poly_m.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.m.mono)/mu.m.poly_ingest,
           EC_poly_m.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.m.mono)/mu.m.poly_trans) %>%
    #calculate environmentally realistic effect threshold
    mutate(EC_env_m.particles.mL_trans = EC_poly_m.particles.mL_trans * CF_bio_trans,
           EC_env_m.particles.mL_ingest = EC_poly_m.particles.mL_ingest * CF_bio_ingest) %>% 
  
    #################################################
    ############# SPECIFIC SURFACE AREA ERM #########
    ####################################################
  ##### specific surface area ERM ####
  mutate(mu.ssa.mono = mu.sa.mono/mu.m.mono) %>% #define mu_x_mono for alignment to ERM (um^2/ug)
    ### translocation-limited ####
    #calculate lower translocatable 1/SSA
    mutate(x_LL_ssa_trans = SSA.inversefnx(sa = x_LL_sa_trans, #surface area
                                           m = x_LL_m_trans), #mass
           x_UL_ssa_trans = SSA.inversefnx(sa = x_UL_sa_trans, #surface area
                                           m = x_UL_m_trans) #mass
    ) %>% 
    #calculate mu_x_poly for specific surface area
    #note that mu were calcaulted for polydisperse particles before, so not special case needed here
    mutate(mu.ssa.inverse.poly_trans = mux_polyfnx(a.ssa, x_UL_ssa_trans, x_LL_ssa_trans)) %>% 
    #calculate polydisperse effect concentration for specific surface area (particles/mL)
    mutate(mu.ssa.poly_trans = 1 / mu.ssa.inverse.poly_trans) %>%  #calculate mu_SSA from inverse
    mutate(EC_poly_ssa.particles.mL_trans = (EC_mono_p.particles.mL_trans * mu.ssa.mono)/mu.ssa.poly_trans) %>% 
    #calculate environmentally realistic effect threshold
    mutate(EC_env_ssa.particles.mL_trans = EC_poly_ssa.particles.mL_trans * CF_bio_trans) %>% 
    
    ### Ingestion-limited ####
  #calculate lower ingestlocatable 1/SSA
  mutate(x_LL_ssa_ingest = SSA.inversefnx(sa = x_LL_sa_ingest, #surface area
                                         m = x_LL_m_ingest), #mass
         x_UL_ssa_ingest = SSA.inversefnx(sa = x_UL_sa_ingest, #surface area
                                         m = x_UL_m_ingest) #mass
  ) %>% 
    #calculate mu_x_poly for specific surface area
    #note that mu were calcaulted for polydisperse particles before, so not special case needed here
    mutate(mu.ssa.inverse.poly_ingest = mux_polyfnx(a.ssa, x_UL_ssa_ingest, x_LL_ssa_ingest)) %>% 
    #calculate polydisperse effect concentration for specific surface area (particles/mL)
    mutate(mu.ssa.poly_ingest = 1 / mu.ssa.inverse.poly_ingest) %>%  #calculate mu_SSA from inverse
    mutate(EC_poly_ssa.particles.mL_ingest = (EC_mono_p.particles.mL_ingest * mu.ssa.mono)/mu.ssa.poly_ingest) %>% 
    #calculate environmentally realistic effect threshold
    mutate(EC_env_ssa.particles.mL_ingest = EC_poly_ssa.particles.mL_ingest * CF_bio_ingest) %>% 
    
    ##################################################
    ### Convert to Metrics other than particles/mL ###
  #######################################################
    ## convert all environmentally realistic thresholds to surface area ##
    # particle count to surface area #
    mutate(EC_env_p.um2.mL_ingest =  EC_env_p.particles.mL_ingest * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # surface area to surface area #
    mutate(EC_env_sa.um2.mL_ingest =  EC_env_sa.particles.mL_ingest * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # volume to surface area #
    mutate(EC_env_v.um2.mL_ingest =  EC_env_v.particles.mL_ingest * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>%
    # mass to surface area #
    mutate(EC_env_m.um2.mL_ingest =  EC_env_m.particles.mL_ingest * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # specific surface area to surface area #
    mutate(EC_env_ssa.um2.mL_ingest =  EC_env_ssa.particles.mL_ingest * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    
    ## convert all environmentally realistic thresholds to volume ##
    # particle count to volume #
    mutate(EC_env_p.um3.mL_ingest =  EC_env_p.particles.mL_ingest * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # surface area to volume #
    mutate(EC_env_sa.um3.mL_ingest =  EC_env_sa.particles.mL_ingest * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # volume to volume #
    mutate(EC_env_v.um3.mL_ingest =  EC_env_v.particles.mL_ingest * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>%
    # mass to volume #
    mutate(EC_env_m.um3.mL_ingest =  EC_env_m.particles.mL_ingest * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # specific surface area to volume #
    mutate(EC_env_ssa.um3.mL_ingest =  EC_env_ssa.particles.mL_ingest * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    
    ## convert all environmentally realistic thresholds to mass ##
    # particle count to mass #
    mutate(EC_env_p.ug.mL_ingest =  EC_env_p.particles.mL_ingest * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # surface area to mass #
    mutate(EC_env_sa.ug.mL_ingest =  EC_env_sa.particles.mL_ingest * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # volume to mass #
    mutate(EC_env_v.ug.mL_ingest =  EC_env_v.particles.mL_ingest * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>%
    # mass to mass #
    mutate(EC_env_m.ug.mL_ingest =  EC_env_m.particles.mL_ingest * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # specific surface area to mass #
    mutate(EC_env_ssa.ug.mL_ingest =  EC_env_ssa.particles.mL_ingest * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    
    ## convert all environmentally realistic thresholds to specific surface area ##
    # particle count to specific surface area #
    mutate(EC_env_p.um2.ug.mL_ingest =  EC_env_p.particles.mL_ingest * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # surface area to specific surface area #
    mutate(EC_env_sa.um2.ug.mL_ingest =  EC_env_sa.particles.mL_ingest * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # volume to specific surface area #
    mutate(EC_env_v.um2.ug.mL_ingest =  EC_env_v.particles.mL_ingest * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>%
    # mass to specific surface area #
    mutate(EC_env_m.um2.ug.mL_ingest =  EC_env_m.particles.mL_ingest * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # specific surface area to specific surface area #
    mutate(EC_env_ssa.um2.ug.mL_ingest =  EC_env_ssa.particles.mL_ingest * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>%
    
    #### TRANSLOCATION ####
  # particle count to surface area #
  mutate(EC_env_p.um2.mL_trans =  EC_env_p.particles.mL_trans * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # surface area to surface area #
    mutate(EC_env_sa.um2.mL_trans =  EC_env_sa.particles.mL_trans * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # volume to surface area #
    mutate(EC_env_v.um2.mL_trans =  EC_env_v.particles.mL_trans * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>%
    # mass to surface area #
    mutate(EC_env_m.um2.mL_trans =  EC_env_m.particles.mL_trans * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # specific surface area to surface area #
    mutate(EC_env_ssa.um2.mL_trans =  EC_env_ssa.particles.mL_trans * mux_polyfnx(a.x = a.sa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    
    ## convert all environmentally realistic thresholds to volume ##
    # particle count to volume #
    mutate(EC_env_p.um3.mL_trans =  EC_env_p.particles.mL_trans * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # surface area to volume #
    mutate(EC_env_sa.um3.mL_trans =  EC_env_sa.particles.mL_trans * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # volume to volume #
    mutate(EC_env_v.um3.mL_trans =  EC_env_v.particles.mL_trans * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>%
    # mass to volume #
    mutate(EC_env_m.um3.mL_trans =  EC_env_m.particles.mL_trans * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # specific surface area to volume #
    mutate(EC_env_ssa.um3.mL_trans =  EC_env_ssa.particles.mL_trans * mux_polyfnx(a.x = a.v, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    
    ## convert all environmentally realistic thresholds to mass ##
    # particle count to mass #
    mutate(EC_env_p.ug.mL_trans =  EC_env_p.particles.mL_trans * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # surface area to mass #
    mutate(EC_env_sa.ug.mL_trans =  EC_env_sa.particles.mL_trans * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # volume to mass #
    mutate(EC_env_v.ug.mL_trans =  EC_env_v.particles.mL_trans * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>%
    # mass to mass #
    mutate(EC_env_m.ug.mL_trans =  EC_env_m.particles.mL_trans * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # specific surface area to mass #
    mutate(EC_env_ssa.ug.mL_trans =  EC_env_ssa.particles.mL_trans * mux_polyfnx(a.x = a.m, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    
    ## convert all environmentally realistic thresholds to specific surface area ##
    # particle count to specific surface area #
    mutate(EC_env_p.um2.ug.mL_trans =  EC_env_p.particles.mL_trans * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # surface area to specific surface area #
    mutate(EC_env_sa.um2.ug.mL_trans =  EC_env_sa.particles.mL_trans * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # volume to specific surface area #
    mutate(EC_env_v.um2.ug.mL_trans =  EC_env_v.particles.mL_trans * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>%
    # mass to specific surface area #
    mutate(EC_env_m.um2.ug.mL_trans =  EC_env_m.particles.mL_trans * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    # specific surface area to specific surface area #
    mutate(EC_env_ssa.um2.ug.mL_trans =  EC_env_ssa.particles.mL_trans * mux_polyfnx(a.x = a.ssa, x_UL = x2D_set, x_LL = x1D_set)) %>% 
    
    #annotate aligned ERM of interest for user interpretability
    mutate("Surface-Area Aligned, Ingestion-Limited, Exposure Concentration (particles/mL)" = EC_env_sa.particles.mL_ingest,
           "Volume Aligned, Ingestion-Limited, Exposure Concentration (particles/mL)" = EC_env_v.particles.mL_ingest,
           "Mass Aligned, Ingestion-Limited, Exposure Concentration (particles/mL)" = EC_env_m.particles.mL_ingest,
           "Specific Surface Area Aligned, Ingestion-Limited, Exposure Concentration (particles/mL)" = EC_env_ssa.particles.mL_ingest,
           
           "Surface-Area Aligned, Translocation-Limited, Exposure Concentration (particles/mL)" = EC_env_sa.particles.mL_trans,
           "Volume Aligned, Translocation-Limited, Exposure Concentration (particles/mL)" = EC_env_v.particles.mL_trans,
           "Mass Aligned, Translocation-Limited, Exposure Concentration (particles/mL)" = EC_env_m.particles.mL_trans,
           "Specific Surface Area Aligned, Translocation-Limited, Exposure Concentration (particles/mL)" = EC_env_ssa.particles.mL_trans
           ) %>% 
    #nudge to front
    dplyr::relocate("Surface-Area Aligned, Ingestion-Limited, Exposure Concentration (particles/mL)",
                    "Volume Aligned, Ingestion-Limited, Exposure Concentration (particles/mL)",
                    "Mass Aligned, Ingestion-Limited, Exposure Concentration (particles/mL)",
                    "Specific Surface Area Aligned, Ingestion-Limited, Exposure Concentration (particles/mL)",
                    
                    "Surface-Area Aligned, Translocation-Limited, Exposure Concentration (particles/mL)",
                    "Volume Aligned, Translocation-Limited, Exposure Concentration (particles/mL)",
                    "Mass Aligned, Translocation-Limited, Exposure Concentration (particles/mL)",
                    "Specific Surface Area Aligned, Translocation-Limited, Exposure Concentration (particles/mL)"
                    ) 
  
  return(df_aligned)
}

cat("Alignment function and dependencies loaded!")

#### EXAMPLE APPLICATION ####
# aoc_z <- readRDS("data/input/aoc_z_tomex2.RDS")
# 
# # step 1: prepare data for alignment
# aoc_prepared <- preparation_fxn(aoc_z)
# # step 2: align data
# aoc_aligned <- alignment_fxn(aoc_prepared)
# 
# #inspect
# aoc_aligned %>%
#   select(Group, Species, size.length.um.used.for.conversions, shape_f, dose.particles.mL.master, `Volume Aligned, Ingestion-Limited, Exposure Concentration (particles/mL)`, 
#          `Surface-Area Aligned, Translocation-Limited, Exposure Concentration (particles/mL)`) %>%
#   head() %>%
#   t()
