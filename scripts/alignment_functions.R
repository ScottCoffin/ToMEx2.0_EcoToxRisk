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
massfnx = function(v, p){
  mass = p * #density (g/cm^3)
    v *
    1/1e12 * 1e6 #correction factor
  return(mass)}

######################################################
############### MASTER ALIGNMENT FUNCTION ###############
########################################################

alignment_fxn <- function(df, #dataframe input
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
  df_aligned <- df %>% 
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
    mutate(dose.particles.mL.AF.corrected = dose.particles.mL.master / (af.time * af.noec)) %>%
    # if NA, then there's missing AFs - data are considered invalid #
    drop_na(dose.particles.mL.AF.corrected) %>%
    # correct from paritcles/mL to particles/L
    mutate(dose.particles.L = dose.particles.mL.AF.corrected * 1000) %>% 
    # ensure no nanoparticle (or simply particles outside desired range) studies are used
    filter(size.length.um.used.for.conversions >= x1D_set,
           size.length.um.used.for.conversions <= x2D_set) %>% 
    mutate(nanoparticle_polydisperse = case_when(
      polydispersity == "polydisperse" & size.length.min.um.used.for.conversions < x1D_set ~ "nanoparticle_filter out",
      T ~ "don't filter out")) %>% 
    filter(nanoparticle_polydisperse == "don't filter out") %>% 
    ##### define environment-specific alpha parameters ############### 
    mutate(environment = env_f) %>% 
    mutate(alpha = case_when(environment == "Marine" ~ alpha.marine,
                             environment == "Freshwater" ~ alpha.freshwater),
           a.sa = case_when(environment == "Marine" ~ a.sa.marine,
                            environment == "Freshwater" ~ a.sa.freshwater),
           a.v = case_when(environment == "Marine" ~ a.v.marine,
                           environment == "Freshwater" ~ a.v.freshwater),
           a.m = case_when(environment == "Marine" ~ a.m.marine,
                           environment == "Freshwater" ~ a.m.freshwater),
           a.ssa = case_when(environment == "Marine" ~ a.ssa.marine,
                             environment == "Freshwater" ~ a.ssa.freshwater),
           R.ave = case_when(environment == "Marine" ~ R.ave.water.marine,
                             environment == "Freshwater" ~ R.ave.water.freshwater),
           p.ave = case_when(environment == "Marine" ~ p.ave.marine,
                             environment == "Freshwater" ~ p.ave.freshwater),
           H_W_ratio = H_W_ratio # currently just using one value for all environments
           ) %>% 
    ################## BIOAVAILABILITY ##############
    #### Estimate ingestible plastic size
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
  mutate(dose.particles.L_trans = case_when(
    translocatable_poly == "translocatable (some)" ~ CF_bioavailable_trans * dose.particles.L,
    T ~ dose.particles.L)) %>% 
  # determine fraciton of partially ingestible particles using correction factor
  mutate(CF_bioavailable_ingest = case_when(
    ingestible_poly == "ingestible (some)" ~ CFfnx(a = alpha,
                                                   x1D = size.length.min.um.used.for.conversions,
                                                   x2D = x2M_ingest,
                                                   x1M = size.length.min.um.used.for.conversions,
                                                   x2M = size.length.max.um.used.for.conversions),
    T ~ 1 )) %>% 
  # calculate ingestible particle effect concentration using correction factoor (fraction)
  mutate(dose.particles.L_ingest = case_when(
    ingestible_poly == "ingestible (some)" ~ CF_bioavailable_ingest * dose.particles.L,
    T ~ dose.particles.L)) %>% 
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
    #ingestion-limited lower/upper bioavailability limits
    mutate(x_LL_v_ingest = volumefnx(length = x1D_set,
                                     width = x1D_set,
                                     height = x1D_set),
           x_UL_v_ingest = volumefnx(length = x2M_ingest, 
                                     width = x2M_ingest,
                                     height = x2M_ingest)) %>%
    # translate to environmental
    mutate(mu.v.poly_ingest = mux_polyfnx(a.v, x_UL_v_ingest, x_LL_v_ingest),
           EC_poly_v.particles.L_ingest = (dose.particles.L_ingest * mu.v.mono)/mu.v.poly_ingest,
           EC_env_v.particles.L_ingest = EC_poly_v.particles.L_ingest * CF_bio_ingest,
           # final Food dilution ERM excludes algae
           particles.L.food.dilution = case_when(
             Group == "Algae" ~ NA,
             T ~ EC_env_v.particles.L_ingest)) %>% 
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
         EC_poly_sa.particles.L_trans = (dose.particles.L_trans  * mu.sa.mono)/mu.sa.poly_trans,
         # calculate EC_env_sa (translocation limited). This is the Tissue Translocation ERM Value!
         EC_env_sa.particles.L_trans = EC_poly_sa.particles.L_trans * CF_bio_trans,
         particles.L.tissue.translocation = EC_env_sa.particles.L_trans)
  
  return(df_aligned)
}


#### EXAMPLE APPLICATION ####
# aoc_z <- readRDS("data/input/aoc_z_tomex2.RDS") 
# 
# aoc_aligned <- alignment_fxn(aoc_z)
# 
# #inspect
# aoc_aligned %>% 
#   select(Group, Species, size.length.um.used.for.conversions, shape_f, dose.particles.L, particles.L.food.dilution, particles.L.tissue.translocation) %>% 
#   head() %>% 
#   t()
