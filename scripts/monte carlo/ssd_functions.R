##### SSD FUNCTIONS #####

library(ssdtools)

#### define function for SSD generation for tier 1 ####
SSD_function_t1 <- function(filtered.data, hcxlcl){
  set.seed(99)
  #data collapse
  collapsed <- filtered.data %>% 
    group_by(Species, Group) %>% 
    summarize(Conc = quantile(dose_new, 0.25),
              .groups = "drop")
  # metadata
  metadata <- filtered.data %>% 
    summarize(n_species = n_distinct(Species),
              n_groups = n_distinct(Group),
              n_datapoints = n(),
              .groups = "drop")
  # If there are less than 6 rows, create and return a list with NA and metadata
  if (nrow(collapsed) < 6) {
    result <- list(
      hc5lcl = NA,
      n_species = metadata$n_species,
      n_groups = metadata$n_groups,
      n_datapoints = metadata$n_datapoints
    )
  } else {
  #fit distributions
  dists <- ssd_fit_dists(collapsed, left = "Conc", dists = c("llogis", "lnorm", #"gamma",
                                                             "lgumbel"), computable = FALSE, silent = FALSE) 
  #use average distribution with weighthing based on AICC
  preds <- predict(dists, average = TRUE, 
                   #ic = "aicc",
                   nboot = nboot, ci= TRUE) 
  #report HC metrics of interest
  hc5lcl <- c(preds$lcl[hcxlcl]) #CI05
  #values to extract
  list(hc5lcl = hc5lcl,
       n_species = metadata$n_species,
       n_groups = metadata$n_groups,
       n_datapoints = metadata$n_datapoints)
  } #ifelse
}# function

# Leave-one-out by study function
sensitivity_t1 <- Vectorize(function(x, y) {
  train <- filtered.data[!filtered.data$doi %in% x,]
  SSD_function_t1(train, hcxlcl = y)
})

#### define function for SSD generation for tier 2 ####
SSD_function_t2 <- function(filtered.data, hcx){
  set.seed(99)
  #data collapse
  collapsed <- filtered.data %>% 
    group_by(Species, Group) %>% 
    summarize(Conc = quantile(dose_new, 0.25),
              .groups = "drop")
  # metadata
  metadata <- filtered.data %>% 
    summarize(n_species = n_distinct(Species),
              n_groups = n_distinct(Group),
              n_datapoints = n(),
              .groups = "drop")
  
  # If there are less than 6 rows, create and return a list with NA and metadata
  if (nrow(collapsed) < 6) {
    result <- list(
      hcx_est = NA,
      hcx05cl = NA,
      hcx95cl = NA,
      n_species = metadata$n_species,
      n_groups = metadata$n_groups,
      n_datapoints = metadata$n_datapoints
    )
  } else {
    
  #fit distributions
  dists <- ssd_fit_dists(collapsed, left = "Conc", 
                         dists = c( "llogis", "lnorm", #"gamma", 
                                    "lgumbel"), computable = FALSE, silent = FALSE) 
  #use average distribution with weighthing based on AICC
  preds <- predict(dists, average = TRUE,# ic = "aicc",
                   nboot = nboot, ci= TRUE) 
  #report HC metrics of interest
  hcx_est <- c(preds$est[hcx]) #HC5
  hcx05cl <- c(preds$lcl[hcx]) #CI05
  hcx95cl <- c(preds$ucl[hcx]) #CI95
  #values to extract
  list(hcx_est = hcx_est,
       hcx05cl = hcx05cl,
       hcx95cl = hcx95cl,
       n_species = metadata$n_species,
       n_groups = metadata$n_groups,
       n_datapoints = metadata$n_datapoints)
  } #close ifelse
} #close function

#### define function for SSD generation for tiers 3 and 4 ####
SSD_function_t3_4 <- function(filtered.data, hcx){
  set.seed(99)
  #data collapse
  collapsed <- filtered.data %>% 
    #filter specific things for tiers 3 and 4
    filter(risk.13 != 1,
           bio_f %in% c("Organism", "Population")) %>% 
    group_by(Species, Group) %>% 
    summarize(Conc = quantile(dose_new, 0.50),
              .groups = "drop")
  # metadata
  metadata <- filtered.data %>% 
    #filter specific things for tiers 3 and 4
    filter(risk.13 != 1,
           bio_f %in% c("Organism", "Population")) %>% 
    summarize(n_species = n_distinct(Species),
              n_groups = n_distinct(Group),
              n_datapoints = n(),
              .groups = "drop")
  
  # If there are less than 6 rows, create and return a list with NA and metadata
  if (nrow(collapsed) < 6) {
    result <- list(
      hcx_est = NA,
      hcx05cl = NA,
      hcx95cl = NA,
      n_species = metadata$n_species,
      n_groups = metadata$n_groups,
      n_datapoints = metadata$n_datapoints
    )
  } else {
    
  #fit distributions
  dists <- ssd_fit_dists(collapsed, left = "Conc", dists = c( "llogis", "lnorm",# "gamma", 
                                                              "lgumbel"), computable = FALSE, silent = FALSE) 
  #use average distribution with weighthing based on AICC
  preds <- predict(dists, average = TRUE, 
                   #ic = "aicc", 
                   nboot = nboot, ci= TRUE) 
  #report HC metrics of interest
  hcx_est <- c(preds$est[hcx]) #HC5
  hcx05cl <- c(preds$lcl[hcx]) #CI05
  hcx95cl <- c(preds$ucl[hcx]) #CI95
  #values to extract
  list(hcx_est = hcx_est,
       hcx05cl = hcx05cl,
       hcx95cl = hcx95cl,
       n_species = metadata$n_species,
       n_groups = metadata$n_groups,
       n_datapoints = metadata$n_datapoints)
  } #ifelse
} #function


## function to generate return filtered data ready for SSDs
filter_environment_data <- function(data, env_filter, upper.tissue.trans.size.um, x1D_set, x2D_set) {
  # Filter and process data for Tissue Translocation
  filtered_data_small_default_t1_2 <- data %>%
    ungroup() %>% 
    mutate(dose_new = particles.mL.ox.stress / (af.time * af.noec)) %>%
    drop_na(dose_new) %>%
    mutate(dose_new = dose_new * 1000) %>% # Convert particles/mL to particles/L
    filter(between(size.length.um.used.for.conversions, x1D_set, upper.tissue.trans.size.um),
           shape_f != "Not Reported",
           poly_f != "Not Reported",
           environment %in% env_filter,
           Group != "Bacterium",
           Group != "Plant",
           effect.metric != "HONEC",
           #only consider studies in which particles are below x2M_trans (translocation only!)
           translocatable == "translocatable"
    ) 
  
  filtered_data_small_default_t3_4 <- filtered_data_small_default_t1_2 %>%
    filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
  
   # Filter and process data for Food Dilution
  filtered_data_large_default_t1_2 <- data %>%
    filter(Group != "Algae") %>%
    mutate(dose_new = particles.mL.food.dilution / (af.time * af.noec)) %>%
    drop_na(dose_new) %>%
    mutate(dose_new = dose_new * 1000) %>% # Convert particles/mL to particles/L
    filter(between(size.length.um.used.for.conversions, x1D_set, x2D_set),
           poly_f != "Not Reported",
           environment %in% env_filter,
           Group != "Bacterium",
           Group != "Plant",
           effect.metric != "HONEC",
           ingestible == "ingestible") #only consider studies in which particles are below x2M_ingest )
  
  filtered_data_large_default_t3_4 <- filtered_data_large_default_t1_2 %>%
    filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
  
  
  return(list(filtered_data_small_default_t1_2 =  filtered_data_small_default_t1_2,
              filtered_data_small_default_t3_4 = filtered_data_small_default_t3_4,
              filtered_data_large_default_t1_2 = filtered_data_large_default_t1_2,
              filtered_data_large_default_t3_4 = filtered_data_large_default_t3_4
          ))
}



## function to generate threshold based on different environment filterings (Marine or Freshwater,etc.)
process_environment_data <- function(data, 
                                     env_filter = "Marine",
                                     upper.tissue.trans.size.um = 88, 
                                     x1D_set = 1,
                                     x2D_set = 5000
                                     ) {
  
  
  ###  only if it doesn't already exist (ensures compatibility with Shiny app) ##
  if (!"particles.mL.ox.stress" %in% names(data)) {
    data <- data %>% 
      mutate(particles.mL.ox.stress = EC_env_sa.particles.mL_trans)
  }
  
  if (!"particles.mL.food.dilution" %in% names(data)) {
    data <- data %>% 
      mutate(particles.mL.food.dilution = EC_env_v.particles.mL_ingest)
  }

  
  # Filter and process data for Tissue Translocation
  filtered_data_small_default_t1_2 <- data %>%
    ungroup() %>% 
    mutate(dose_new = particles.mL.ox.stress / (af.time * af.noec)) %>%
    drop_na(dose_new) %>%
    mutate(dose_new = dose_new * 1000) %>% # Convert particles/mL to particles/L
    filter(
          between(size.length.um.used.for.conversions, x1D_set, upper.tissue.trans.size.um),
           shape_f != "Not Reported",
           poly_f != "Not Reported",
           environment %in% env_filter,
           Group != "Bacterium",
           Group != "Plant",
           effect.metric != "HONEC",
          ingestible == "ingestible",
           #only consider studies in which particles are below x2M_trans (translocation only!)
           translocatable == "translocatable",
          dose_new > 0
    ) 
  
  filtered_data_small_default_t3_4 <- filtered_data_small_default_t1_2 %>%
    filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
  
  # Compute thresholds for tissue translocation
  small_default_t1 <- SSD_function_t1(filtered.data = filtered_data_small_default_t1_2, hcxlcl = 5)
  small_default_t2 <- SSD_function_t2(filtered.data = filtered_data_small_default_t1_2, hcx = 5)
  small_default_t3 <- SSD_function_t3_4(filtered.data = filtered_data_small_default_t3_4, hcx = 5)
  small_default_t4 <- SSD_function_t3_4(filtered.data = filtered_data_small_default_t3_4, hcx = 10)
  
  # Filter and process data for Food Dilution
  filtered_data_large_default_t1_2 <- data %>%
    filter(Group != "Algae") %>%
    mutate(dose_new = particles.mL.food.dilution / (af.time * af.noec)) %>%
    drop_na(dose_new) %>%
    mutate(dose_new = dose_new * 1000) %>% # Convert particles/mL to particles/L
    filter(between(size.length.um.used.for.conversions, x1D_set, x2D_set),
           poly_f != "Not Reported",
           environment %in% env_filter,
           Group != "Bacterium",
           Group != "Plant",
           effect.metric != "HONEC",
           ingestible == "ingestible", #only consider studies in which particles are below x2M_ingest )
           dose_new > 0) 
  
  filtered_data_large_default_t3_4 <- filtered_data_large_default_t1_2 %>%
    filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
  
  # Compute thresholds for food dilution
  large_default_t1 <- SSD_function_t1(filtered.data = filtered_data_large_default_t1_2, hcxlcl = 5)
  large_default_t2 <- SSD_function_t2(filtered.data = filtered_data_large_default_t1_2, hcx = 5)
  large_default_t3 <- SSD_function_t3_4(filtered.data = filtered_data_large_default_t3_4, hcx = 5)
  large_default_t4 <- SSD_function_t3_4(filtered.data = filtered_data_large_default_t3_4, hcx = 10)
  
  # Compile thresholds into a tibble
  base_thresholds <- tibble(
    Tier = c('Tier1', 'Tier2', 'Tier3', 'Tier4'),
    "Tissue Translocation (Default)" = c(small_default_t1$hc5lcl, small_default_t2$hcx_est, small_default_t3$hcx_est, small_default_t4$hcx_est),
    "Food Dilution (Default)" = c(large_default_t1$hc5lcl, large_default_t2$hcx_est, large_default_t3$hcx_est, large_default_t4$hcx_est),
    "Tissue Translocation (5th %)" = c(NA, small_default_t2$hcx05cl, small_default_t3$hcx05cl, small_default_t4$hcx05cl),
    "Tissue Translocation (95th %)" = c(NA, small_default_t2$hcx95cl, small_default_t3$hcx95cl, small_default_t4$hcx95cl),
    "Food Dilution (5th %)" = c(NA, large_default_t2$hcx05cl, large_default_t3$hcx05cl, large_default_t4$hcx05cl),
    "Food Dilution (95th %)" = c(NA, large_default_t2$hcx95cl, large_default_t3$hcx95cl, large_default_t4$hcx95cl),
    "N Species (tissue trans)" = c(small_default_t1$n_species, small_default_t2$n_species, small_default_t3$n_species, small_default_t4$n_species),
    "N Species (food dilution)" = c(large_default_t1$n_species, large_default_t2$n_species, large_default_t3$n_species, large_default_t4$n_species),
    "N Groups (tissue trans)" = c(small_default_t1$n_groups, small_default_t2$n_groups, small_default_t3$n_groups, small_default_t4$n_groups),
    "N Groups (food dilution)" = c(large_default_t1$n_groups, large_default_t2$n_groups, large_default_t3$n_groups, large_default_t4$n_groups),
    "N Datapoints (tissue trans)" = c(small_default_t1$n_datapoints, small_default_t2$n_datapoints, small_default_t3$n_datapoints, small_default_t4$n_datapoints),
    "N Datapoints (food dilution)" = c(large_default_t1$n_datapoints, large_default_t2$n_datapoints, large_default_t3$n_datapoints, large_default_t4$n_datapoints)
  )
  
  
  return(base_thresholds)
}

# # # Assuming aoc_risk_paper is your initial dataset
 # marine_thresholds <- process_environment_data(aoc_risk_paper, "Marine", 83, 1, 5000)
 # freshwater_thresholds <- process_environment_data(aoc_risk_paper, "Freshwater", 83, 1, 5000)



# ##### Base Thresholds
# #filter out risk criteria (not done above)#
# aoc_risk_paper <- aoc_final %>%
#   drop_na(effect.metric) %>%
#   filter(tier_zero_tech_f == ("Red Criteria Passed"))
# 
# ####---- TISSUE TRANSLOCATION ------#####
# filtered.data.small.default_t1.2 <- aoc_risk_paper %>%
#   mutate(dose_new = particles.mL.ox.stress / (af.time * af.noec)) %>%
#   drop_na(dose_new) %>%
#   mutate(dose_new = dose_new * 1000) %>% #convert particles/mL to particles/L
#   filter(between(size.length.um.used.for.conversions, 1, upper.tissue.trans.size.um),
#          shape_f != "Not Reported",
#          poly_f != "Not Reported",
#          !environment %in% c("Terrestrial", "Not Reported"),
#          Group != "Bacterium",
#          Group != "Plant",
#          effect.metric != "HONEC")
# 
# filtered.data.small.default_t3.4 <- filtered.data.small.default_t1.2 %>%
#   filter(risk.13 != 1,  bio_f %in% c("Organism", "Population"))
# 
# # get thresholds
# small.default.t1 <- SSD_function_t1(filtered.data = filtered.data.small.default_t1.2, hcxlcl = 5)
# small.default.t2 <- SSD_function_t2(filtered.data = filtered.data.small.default_t1.2, hcx = 5)
# small.default.t3 <- SSD_function_t3_4(filtered.data = filtered.data.small.default_t3.4, hcx = 5)
# small.default.t4 <- SSD_function_t3_4(filtered.data = filtered.data.small.default_t3.4, hcx = 10)
# 
# ####---- Food Dilution ------#####
# filtered.data.large.default_t1.2 <- aoc_risk_paper %>%
#   # remove algae, as food dilution MOE doesn't make sense for algae
#   filter(Group != "Algae") %>%
#   mutate(dose_new = particles.mL.food.dilution / (af.time * af.noec)) %>%
#   drop_na(dose_new) %>%
#   mutate(dose_new = dose_new * 1000) %>% #convert particles/mL to particles/L
#   filter(between(size.length.um.used.for.conversions, x1D_set, x2D_set),
#          poly_f != "Not Reported",
#          !environment %in% c("Terrestrial", "Not Reported"),
#          Group != "Bacterium",
#          Group != "Plant",
#          effect.metric != "HONEC")
# 
# filtered.data.large.default_t3.4 <- filtered.data.large.default_t1.2 %>%
#   filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))
# 
# # get thresholds
# large.default.t1 <- SSD_function_t1(filtered.data = filtered.data.large.default_t1.2, hcxlcl = 5)
# large.default.t2 <- SSD_function_t2(filtered.data = filtered.data.large.default_t1.2, hcx = 5)
# large.default.t3 <- SSD_function_t3_4(filtered.data = filtered.data.large.default_t3.4, hcx = 5)
# large.default.t4 <- SSD_function_t3_4(filtered.data = filtered.data.large.default_t3.4, hcx = 10)
# 
# base_thresholds <- tibble(
#   "Tier" = c('Tier1', 'Tier2', 'Tier3', 'Tier4'),
#   "Tissue Translocation (Default)" = c(small.default.t1, small.default.t2$hcx_est, small.default.t3$hcx_est, small.default.t4$hcx_est),
#   "Food Dilution (Default)" = c(large.default.t1, large.default.t2$hcx_est, large.default.t3$hcx_est, large.default.t4$hcx_est),  "Tissue Translocation (5th %)" = c(NA, small.default.t2$hcx05cl, small.default.t3$hcx05cl, small.default.t4$hcx05cl),
#   "Tissue Translocation (95th %)" = c(NA, small.default.t2$hcx95cl, small.default.t3$hcx95cl, small.default.t4$hcx95cl),
#   "Food Dilution (5th %)" = c(NA, large.default.t2$hcx05cl, large.default.t3$hcx05cl, large.default.t4$hcx05cl),
#   "Food Dilution (95th %)" = c(NA, large.default.t2$hcx95cl, large.default.t3$hcx95cl, large.default.t4$hcx95cl))
# 
# base_thresholds