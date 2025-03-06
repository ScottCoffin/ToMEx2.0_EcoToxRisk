#ToMEx 2.0 Data Tidying Script
#Date: 5/09/23
#Created By: Leah Thornton Hampton_ MODIFIED by Scott Coffin
#Description: Script to join validated data templates from ToMEx 2.0 exercise, re-structure to match ToMEx 1.0, bind to existing ToMEx 1.0 Database

#Load Packages
library(tidyverse)
library(readr)

#Set working directory
#setwd("~aq_mp_tox_shiny/")
#source("scripts/functions.R") # necessary for surface area, volume calculations

#### Extract Data from Submitted Templates ####
#### NOTE: ALREADY DONE IN aq_mp_tox_shiny REPO! JUST IMPORT tomex2_input.rds
#### tomex2_input.rds generated in ToMEx2.0_AquaticOrganisms/monte carlo/Master_Data_Tidying_Monte_Carlo and stored in ref data/ folder
#tomex2.0 <- read_rds("scripts/monte carlo/ref data/tomex2_input.rds")

#Read in ToMEx 1.0 Tidy Data sets
#aoc_setup <- readRDS("aoc_setup.RDS")

##### AOC SETUP #####

# R.ave.water.marine <- 0.77 # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
# R.ave.water.freshwater <- 0.67
# R.ave.sediment.marine <- 0.75
# R.ave.sediment.freshwater <- 0.70

ToMEx2.0fxn <- function(aoc_setup, #ToMEx 1.0 dataset
                        tomex2.0, #ToMEx 2.0 dataset
                        bodysize_addons, #datasets as input
                        R.ave.water.marine, R.ave.water.freshwater, R.ave.sediment.freshwater, R.ave.sediment.marine,
                        beta_log10_body_length, body_length_intercept){

tomex2.0_aoc_setup <- tomex2.0 %>%
   #Remove effect metrics when there are less than 3 treatments
   mutate(`effect metric` = ifelse(treatments < 3, NA_character_, `effect metric`)) %>% 
   #Replace NAs
   replace_na(list(shape = "Not Reported", polymer = "Not Reported", life.stage = "Not Reported")) %>% 
   #Add source column
   add_column(source = "ToMEx 2.0", .before = "doi") %>% 
   #Use only last name of first author
   mutate(authors = word(authors,1,sep = ",")) %>% 
   arrange(authors) %>% 
   base::transform(article=as.numeric(factor(doi))+162) %>% 
   relocate(article, .after = doi) %>%   
   #Add rowid column
   rowid_to_column() %>% 
   #Move screening scores up
   relocate(100:129, .after = article) %>% 
  ### Previously, all column names were renamed using base::transform() which replaced spaces with periods,...
  ### that package was updated to no longer do that, so we'll need to manually rename everything to be compatible
  ### with all downstream scripts.
  
  rename_with(~ gsub(" |[?]|[(]|[)]|-|/|\\^", ".",.x), everything()) %>% ### replaces all special chars with '.'
  
   #Change pass/fail descriptors to 2/0
   mutate(`tech.a1` = case_when(`tech.a1` == "Pass" ~ 2,
                                `tech.a1` == "Fail" ~ 0)) %>%
   mutate(`tech.a2` = case_when(`tech.a2` == "Pass" ~ 2,
                                `tech.a2` == "Fail" ~ 0)) %>%
   mutate(`tech.a3` = case_when(`tech.a3` == "Pass" ~ 2,
                                `tech.a3` == "Fail" ~ 0)) %>%
   mutate(`tech.a4` = case_when(`tech.a4` == "Pass" ~ 2,
                                `tech.a4` == "Fail" ~ 0)) %>%
   mutate(`tech.a5` = case_when(`tech.a5` == "Pass" ~ 2,
                                `tech.a5` == "Fail" ~ 0)) %>%
   mutate(`tech.a6` = case_when(`tech.a6` == "Pass" ~ 2,
                                `tech.a6` == "Fail" ~ 0)) %>%
   mutate(`risk.b1` = case_when(`risk.b1` == "Pass" ~ 2,
                                `risk.b1` == "Fail" ~ 0)) %>% 
   #Change column names of summed scores to match
   rename(technical.quality = `tech.total`,
          risk.quality = `risk.total`,
          total.quality = `overall.total`) %>% 
   #Add columns for red criteria failures
   mutate(tier_zero_tech_f = if_else((`tech.a1` == 0|`tech.a2` == 0|`tech.a3` == 0|`tech.a4` == 0|`tech.a5` == 0|`tech.a6` == 0), "Red Criteria Failed", "Red Criteria Passed")) %>% 
   relocate(tier_zero_tech_f , .after = total.quality) %>% 
   mutate(tier_zero_risk_f  = if_else((`risk.b1` == 0), "Red Criteria Failed", "Red Criteria Passed")) %>% 
   relocate(tier_zero_risk_f, .after = tier_zero_tech_f ) %>%
   replace_na(list(tier_zero_tech_f  = "Scoring Not Applicable", tier_zero_risk_f = "Scoring Not Applicable")) %>%
   mutate(tier_zero_tech_f  = factor(tier_zero_tech_f )) %>% 
   mutate(tier_zero_risk_f = factor(tier_zero_risk_f)) %>%
   #Fix life stages for algae and bacteria
   mutate(life.stage = if_else(organism.group == "Algae", "Not Reported", life.stage)) %>% 
   mutate(life.stage = if_else(organism.group == "Bacterium", "Not Reported", life.stage)) %>%  
   #Change Mytilus NA to Mytilus species
   mutate(species = if_else(species == "Mytilus NA", "Mytilus species", species)) %>% 
   #Factor species
   mutate(species = factor(species)) %>% 
   rename(species_f = species) %>% 
   #Factor organism group
   mutate(organism.group = factor(organism.group)) %>% 
   rename(org_f = organism.group) %>% 
   #Factor in vitro in vivo
   mutate(in.vitro.in.vivo = factor(case_when(
     in.vitro.in.vivo == "In Vivo" ~ "In Vivo",
     in.vitro.in.vivo == "In Vitro" ~ "In Vitro",
     in.vitro.in.vivo == "in vivo" ~ "In Vivo",
     in.vitro.in.vivo == "in vitro" ~ "In Vivo",
   ))) %>% 
   rename(vivo_f = in.vitro.in.vivo) %>%
   #Factor life stage
   mutate(life.stage = factor(life.stage)) %>%  
   rename(life_f = life.stage) %>%
   #Factor environment
   mutate(environment = factor(environment)) %>% 
   rename(env_f = environment) %>%
   #Factor experiment type
   mutate(experiment.type = factor(experiment.type)) %>% 
   rename(exp_type_f = experiment.type) %>%
   #Rename particle mix, reference particle, salinity, temp, exposure duration columns to match
   rename(mix = `particle.mix.`,
          reference.material = `reference.particle`,
          media.sal.ppt = `media.salinity..ppt.`,
          media.temp = `media.temp..mean.`,
          exposure.duration.d = `exposure.duration..days.`) %>%
   #Paste the range of treatments used within each study - there is often a range because this is different depending on the endpoint
   group_by(doi) %>% 
   mutate(treatment_range = ifelse(min(treatments) == max(treatments), paste0(treatments), paste0(min(treatments),"-",max(treatments)))) %>% 
   ungroup() %>% 
   relocate(treatment_range, .after = treatments) %>% 
   
  #Dosing restructuring #ONBOARDING CHECK - ADD UNITS AS NEEDED TO CASE_WHEN STATEMENTS#
   #Count - Nominal
   mutate(dose.particles.mL.nominal = case_when(
     nominal.dose...particles.units == "p/mL" ~ nominal.dose...particles,
     nominal.dose...particles.units == "particles/ml" ~ nominal.dose...particles,
     nominal.dose...particles.units == "particles/mL" ~ nominal.dose...particles,
     nominal.dose...particles.units == "particles/m3" ~ nominal.dose...particles/1000000,
     nominal.dose...particles.units == "particles/L" ~ nominal.dose...particles/1000,
     nominal.dose...particles.units == "particles/l" ~ nominal.dose...particles/1000,
     nominal.dose...particles.units == "L" ~ nominal.dose...particles/1000
     )) %>% 
  relocate(dose.particles.mL.nominal, .after = sample.size) %>% 
  #Count - Measured
  mutate(dose.particles.mL.measured = case_when(
    measured.dose...particles.units == "p/mL" ~ measured.dose...particles,
    measured.dose...particles.units == "particles/ml" ~ measured.dose...particles,
    measured.dose...particles.units == "particles/mL" ~ measured.dose...particles,
    measured.dose...particles.units == "particles/m3" ~ measured.dose...particles/1000000,
    measured.dose...particles.units == "particles/L" ~ measured.dose...particles/1000,
    measured.dose...particles.units == "particles/l" ~ measured.dose...particles/1000,
    measured.dose...particles.units == "L" ~ measured.dose...particles/1000
  )) %>% 
  relocate(dose.particles.mL.measured, .after = dose.particles.mL.nominal) %>% 
   #Mass - Nominal
   mutate(dose.mg.L.nominal = case_when(
     nominal.dose...mass.units == "g/L" ~ nominal.dose...mass*1000,
     nominal.dose...mass.units == "Kg/L" ~ nominal.dose...mass*1000000,
     nominal.dose...mass.units == "mg/L" ~ nominal.dose...mass,
     nominal.dose...mass.units == "mg/l" ~ nominal.dose...mass,
     nominal.dose...mass.units == "ug/mL" ~ nominal.dose...mass,
     nominal.dose...mass.units == "µg/mL" ~ nominal.dose...mass,
     nominal.dose...mass.units == "g/mL" ~ nominal.dose...mass*1000000,
     nominal.dose...mass.units == "mg/mL" ~ nominal.dose...mass*1000,
     nominal.dose...mass.units == "ug/L" ~ nominal.dose...mass/1000,
     nominal.dose...mass.units == "µg/L" ~ nominal.dose...mass/1000,
     nominal.dose...mass.units == "ug/l" ~ nominal.dose...mass/1000,
     nominal.dose...mass.units == "µg/l" ~ nominal.dose...mass/1000,
     nominal.dose...mass.units == "ng/L" ~ nominal.dose...mass/1000000,
     )) %>% 
  relocate(dose.mg.L.nominal, .after = dose.particles.mL.nominal) %>% 
  #Mass - Measured
  mutate(dose.mg.L.measured = case_when(
    measured.dose...mass.units == "g/L" ~ measured.dose...mass*1000,
    measured.dose...mass.units == "Kg/L" ~ measured.dose...mass*1000000,
    measured.dose...mass.units == "mg/L" ~ measured.dose...mass,
    measured.dose...mass.units == "mg/l" ~ measured.dose...mass,
    measured.dose...mass.units == "ug/mL" ~ measured.dose...mass,
    measured.dose...mass.units == "µg/mL" ~ measured.dose...mass,
    measured.dose...mass.units == "g/mL" ~ measured.dose...mass*1000000,
    measured.dose...mass.units == "mg/mL" ~ measured.dose...mass*1000,
    measured.dose...mass.units == "ug/L" ~ measured.dose...mass/1000,
    measured.dose...mass.units == "µg/L" ~ measured.dose...mass/1000,
    measured.dose...mass.units == "ug/l" ~ measured.dose...mass/1000,
    measured.dose...mass.units == "µg/l" ~ measured.dose...mass/1000,
    measured.dose...mass.units == "ng/L" ~ measured.dose...mass/1000000,
  )) %>% 
    relocate(dose.mg.L.measured, .after = dose.particles.mL.measured) %>% 
  #Mass - Nominal (sediment)
  mutate(nominal.dose.mg.kg.sediment = case_when(
    exposure.route == "sediment" & `nominal.dose...alternative.type.units` == "g/kg sediment" ~ `nominal.dose...alternative.type`/1000,
    exposure.route == "sediment" & `nominal.dose...alternative.type.units` == "mg/kg sediment" ~ `nominal.dose...alternative.type`,
    exposure.route == "sediment" & `nominal.dose...alternative.type.units` == "mg/kg sediment dry weight" ~ `nominal.dose...alternative.type`,
    exposure.route == "sediment" & `nominal.dose...alternative.type.units` == "mg/kg" ~ `nominal.dose...alternative.type`,
    exposure.route == "sediment" & `nominal.dose...alternative.type.units` == "g/kg (dw) sediment" ~ `nominal.dose...alternative.type`/1000,
    exposure.route == "sediment" & `nominal.dose...alternative.type.units` == "mg/Kg sediment" ~ `nominal.dose...alternative.type`,
  )) %>%
  mutate(`nominal.dose...alternative.type.units` = ifelse(!is.na(nominal.dose.mg.kg.sediment), NA_real_, `nominal.dose...alternative.type.units`)) %>% #Clear alternative dose column
  mutate(`nominal.dose...alternative.type` = ifelse(!is.na(nominal.dose.mg.kg.sediment), NA_real_, `nominal.dose...alternative.type`)) %>% #Clear alternative dose column
  #Count - Nominal (sediment)
  mutate(nominal.dose.particles.kg.sediment = case_when(
    exposure.route == "sediment" & `nominal.dose...alternative.type.units` == "particles/kg sediment dry weight" ~ `nominal.dose...alternative.type`
  )) %>%  
  #Mass - Measured (sediment)
  mutate(measured.dose.mg.kg.sediment = case_when(
    exposure.route == "sediment" & `measured.dose.alternative.units` == "mg/kg sediment" ~ `measured.dose.alternative`,
  )) %>%  #Place holder for measured mass doses
  #Count - Measured (sediment)
  mutate(measured.dose.particles.kg.sediment = case_when(
    exposure.route == "sediment" & `measured.dose.alternative.units` == "particles/kg sediment dry weight" ~ `measured.dose.alternative`,
    exposure.route == "sediment" & `measured.dose.alternative.units` == "particles/kg" ~ `measured.dose.alternative`,
  )) %>% 
  mutate(`measured.dose.alternative.units` = ifelse(!is.na(measured.dose.particles.kg.sediment), NA_real_, `measured.dose.alternative.units`)) %>% #Clear alternative dose column
  mutate(`measured.dose.alternative` = ifelse(!is.na(measured.dose.particles.kg.sediment), NA_real_, `measured.dose.alternative`)) %>% #Clear alternative dose column
  
  relocate(nominal.dose.mg.kg.sediment, .after = dose.mg.L.measured) %>% 
  relocate(nominal.dose.particles.kg.sediment, .after = nominal.dose.mg.kg.sediment) %>% 
  relocate(measured.dose.mg.kg.sediment, .after = nominal.dose.particles.kg.sediment) %>% 
  relocate(measured.dose.particles.kg.sediment, .after = measured.dose.mg.kg.sediment) %>% 
  #Create master columns for dose and count - measured doses preferred
  mutate(dose.particles.mL.master = if_else(!is.na(dose.particles.mL.measured), dose.particles.mL.measured, dose.particles.mL.nominal)) %>% 
  relocate(dose.particles.mL.master, .after = dose.mg.L.measured) %>% 
  mutate(dose.mg.L.master = if_else(!is.na(dose.mg.L.measured), dose.mg.L.measured, dose.mg.L.nominal)) %>% 
  relocate(dose.mg.L.master, .after = dose.particles.mL.master)  %>% 
  mutate(dose.mg.kg.sediment.master = if_else(!is.na(measured.dose.mg.kg.sediment), measured.dose.mg.kg.sediment, nominal.dose.mg.kg.sediment)) %>% 
  relocate(dose.mg.kg.sediment.master, .after = dose.mg.L.master)  %>% 
  mutate(dose.particles.kg.sediment.master = if_else(!is.na(measured.dose.particles.kg.sediment), measured.dose.particles.kg.sediment, nominal.dose.particles.kg.sediment)) %>% 
  relocate(dose.particles.kg.sediment.master, .after = dose.mg.kg.sediment.master)  %>%
  #Mark that doses were reported (converted doses are to be added later in script)
  mutate(dose.particles.mL.master.converted.reported = if_else(!is.na(dose.particles.mL.master), "reported", NA_character_)) %>% 
  relocate(dose.particles.mL.master.converted.reported, .after = dose.particles.mL.master) %>% 
  mutate(dose.mg.L.master.converted.reported = if_else(!is.na(dose.mg.L.master), "reported", NA_character_)) %>% 
  relocate(dose.mg.L.master.converted.reported, .after = dose.mg.L.master) %>% 
  mutate(dose.mg.kg.sediment.master.converted.reported = if_else(!is.na(dose.mg.kg.sediment.master), "reported", NA_character_)) %>% 
  relocate(dose.mg.kg.sediment.master.converted.reported, .after = dose.mg.kg.sediment.master) %>% 
  mutate(dose.particles.kg.sediment.master.converted.reported = if_else(!is.na(dose.particles.kg.sediment.master), "reported", NA_character_)) %>% 
  relocate(dose.particles.kg.sediment.master.converted.reported, .after = dose.particles.kg.sediment.master) %>% 
  #Chemical Dosing
  rename(chem.add.nominal = nominal.chemicals.added) %>% 
  mutate(chem.add.dose.mg.L.nominal = case_when(
    nominal.added.chemical.dose.units == "mg/L" ~ nominal.added.chemical.dose,
    nominal.added.chemical.dose.units == "ug/L" ~ nominal.added.chemical.dose/1000,
  )) %>% 
  relocate(chem.add.nominal, .after = dose.mg.L.master.converted.reported) %>% 
  relocate(chem.add.dose.mg.L.nominal, .after = chem.add.nominal) %>% 
  
  rename(chem.add.measured = measured.chemicals.added) %>% 
  mutate(chem.add.dose.mg.L.measured = case_when(
    measured.chemical.dose.units == "mg/L" ~ measured.chemical.dose,
    measured.chemical.dose.units == "ug/L" ~ measured.chemical.dose/1000,
  )) %>% 
  relocate(chem.add.measured, .after = chem.add.dose.mg.L.nominal) %>% 
  relocate(chem.add.dose.mg.L.measured, .after = chem.add.measured) %>%  
  #Factor effect
  mutate(effect = factor(effect)) %>% 
  rename(effect_f = effect) %>%
  #Factor effect metric
  mutate(effect.metric = factor(effect.metric)) %>%   
  #Add new columns for assessment factors
  mutate(af.time = case_when(
    org_f == "Algae" & exposure.duration.d < 3 ~ 10,
    org_f == "Algae" & exposure.duration.d >= 3 ~ 1,
    org_f == "Plant" & exposure.duration.d < 28 ~ 10,
    org_f == "Plant" & exposure.duration.d >= 28 ~ 1,
    org_f == "Crustacea" & exposure.duration.d < 21 ~ 10,
    org_f == "Crustacea" & exposure.duration.d >= 21 ~ 1,
    org_f == "Mollusca" & exposure.duration.d < 28 ~ 10,
    org_f == "Mollusca" & exposure.duration.d >= 28 ~ 1,
    org_f == "Fish" & exposure.duration.d < 28 ~ 10,
    org_f == "Fish" & exposure.duration.d >= 28 ~ 1,
    #Assessment factors outside of Wigger et al., 2020 DOI: 10.1002/ieam.4214
    exposure.duration.d < 21 ~ 10,
    exposure.duration.d >= 21 ~ 1)) %>%
  #Define acute and chronic exposures and factor
  mutate(acute.chronic_f = factor(
    if_else(af.time == 1, "Chronic", "Acute"))) %>% 
  #Rename recovery column 
  rename(`Recovery (Days)` = `recovery..days.`) %>% 
  relocate(`Recovery (Days)`, .after = exposure.duration.d) %>% 
  #Assign assessment factors for effect metrics Wigger et al., 2020 DOI: 10.1002/ieam.4214
  mutate(af.noec = case_when(
    effect.metric == "NOEC" ~ 1,
    effect.metric == "LOEC" ~ 2,
    effect.metric == "EC50" ~ 10,
    effect.metric == "EC10" ~ 1,
    effect.metric == "IC10" ~ 1,
    effect.metric == "HONEC" ~ 1,
    effect.metric == "IC50" ~ 10,
    effect.metric == "LC50" ~ 10)) %>% 
  #Factor and rename endpoint columns
  mutate(broad.endpoint.category = factor(broad.endpoint.category)) %>% 
  rename(lvl1_f = broad.endpoint.category) %>%
  mutate(specific.endpoint.category = factor(specific.endpoint.category)) %>% 
  rename(lvl2_f = specific.endpoint.category) %>% 
  mutate(endpoint = factor(endpoint)) %>% 
  rename(lvl3_f = endpoint) %>% 
  #Factor and rename biological level of organization
  mutate(level.of.biological.organization = factor(level.of.biological.organization)) %>% 
  rename(bio_f = level.of.biological.organization) %>%  
  #Rename target cell or tissue column
  rename(target.cell.tissue = target.cell.or.tissue) %>% 
  #Factor polymer
  mutate(polymer = factor(case_when(
    polymer == "Biopolymer" ~ "Biopolymer",
    polymer == "Latex" ~ "Latex",
    polymer == "Not Reported" ~ "Not Reported",
    polymer == "Polyamide" ~ "Polyamide",
    polymer == "Polycarbonate" ~ "Polycarbonate",
    polymer == "Polyethylene" ~ "Polyethylene",
    polymer == "Polyethylene Terephthalate" ~ "Polyethylene Terephthalate",
    polymer == "Polyethylene Vinyl Acetate" ~ "Polyethylene Vinyl Acetate",
    polymer == "Polyisoprene" ~ "Polyisoprene",
    polymer == "Polylactic Acid" ~ "Polylactic Acid",
    polymer == "Polymethylmethacrylate" ~ "Polymethylmethacrylate",
    polymer == "Polyoxymethylene" ~ "Polyoxymethylene",
    polymer == "Polypropylene" ~ "Polypropylene",
    polymer == "Polystyrene" ~ "Polystyrene",
    polymer == "PS" ~ "Polystyrene",
    polymer == "High Density Polyethylene" ~ "High Density Polyethylene",
    polymer == "High density polyethylene (HDPE)" ~ "High Density Polyethylene",
    polymer == "High density polyethylene" ~ "High Density Polyethylene",
    polymer == "Low Density Polyethylene" ~ "Low Density Polyethylene",
    polymer == "Low density polyethylene" ~ "Low Density Polyethylene",
    polymer == "LDPE" ~ "Low Density Polyethylene",
    polymer == "Medium Density Polyethylene" ~ "Medium Density Polyethylene",
    polymer == "Mix - See Original Study" ~ "Mix - See Original Study",
    polymer == "Poly-amidoamine (PAMAM)" ~ "Polyamidoamine",
    polymer == "Poly (Styrene-co-acrylonitrile)" ~ "Poly(Styrene-co-acrylonitrile)",
    polymer == "Polyamide 66" ~ "Polyamide 66",
    polymer == "Polyethylene (co-Vinyl acetate)" ~ "Polyethylene (co-Vinyl acetate)",
    polymer == "Polytetrafluoroethylene" ~ "Polytetrafluoroethylene",
    polymer == "Polyvinyl Acetate" ~ "Polyvinyl Acetate",
    polymer == "Polyvinylchloride" ~ "Polyvinylchloride",
    polymer == "polyvinyl chloride/ vinyl acetate co-polymer" ~ "Polyvinylchloride/vinylacetate co-polymer",
    polymer == "Polyurethane" ~ "Polyurethane",
    polymer == "pristine tire wear particles (P-TWP)" ~ "Tire Wear",
    polymer == "Sodium Polyacrylate" ~ "Sodium Polyacrylate",
    polymer == "Starch/PBAT/PLA" ~ "Starch/Polybutylene Adipate Terephthalate/Polylactic Acid"))) %>% 
  rename(poly_f = polymer) %>%  
  #Rename density column
  rename(density.g.cm3 = density..g.cm.3.) %>% 
  rename(density.reported.estimated = density..reported.estimated.) %>% 
  #Factor shape
  mutate(shape = factor(shape)) %>% 
  rename(shape_f = shape) %>%  
  #Rename zeta potential column
  rename(zetapotential.mV = zeta.potential..mv.,
         zetapotential.media = zeta.potential.media) %>% 
  #Rename size columns
  rename(
    size.length.min.mm.nominal = size.length.mm.nominal..minimum.,
    size.length.max.mm.nominal = size.length.mm.nominal..maximum.,
    size.length.min.mm.measured = size.length.mm.measured..minimum.,
    size.length.max.mm.measured = size.length.mm.measured..maximum.,
    
    size.width.min.mm.nominal = size.width.mm.nominal..minimum.,
    size.width.max.mm.nominal = size.width.mm.nominal..maximum.,
    size.width.min.mm.measured = size.width.mm.measured..minimum.,
    size.width.max.mm.measured = size.width.mm.measured..maximum.,
    
    size.height.min.mm.nominal = size.height.mm.nominal..minimum.,
    size.height.max.mm.nominal = size.height.mm.nominal..maximum.,
    size.height.min.mm.measured = size.height.mm.measured..minimum.,
    size.height.max.mm.measured = size.height.mm.measured..maximum.) %>%
  #Select size lengths to be used for conversions
  mutate(size.length.um.used.for.conversions = case_when(
    !is.na(size.length.mm.measured) ~ size.length.mm.measured*1000,
    !is.na(size.length.min.mm.measured) ~ ((size.length.max.mm.measured + size.length.min.mm.measured)/2)*1000,
    !is.na(size.length.mm.nominal) ~ size.length.mm.nominal*1000,
    !is.na(size.length.min.mm.nominal) ~ ((size.length.max.mm.nominal + size.length.min.mm.nominal)/2)*1000,
    
    !is.na(size.length.max.mm.measured) ~ size.length.max.mm.measured*1000,
    !is.na(size.width.max.mm.nominal) ~ size.width.max.mm.nominal*1000)) %>% 
  
  relocate(size.length.um.used.for.conversions, .after = zetapotential.media) %>% 
  #Select size widths to be used for conversions
  mutate(size.width.um.used.for.conversions = case_when(
    !is.na(size.width.mm.measured) ~ size.width.mm.measured*1000,
    !is.na(size.width.min.mm.measured) ~ ((size.width.max.mm.measured + size.width.min.mm.measured)/2)*1000,
    !is.na(size.width.mm.nominal) ~ size.width.mm.nominal*1000,
    !is.na(size.width.min.mm.nominal) ~ ((size.width.max.mm.nominal + size.width.min.mm.nominal)/2)*1000)) %>% 
  relocate(size.width.um.used.for.conversions, .after = size.length.um.used.for.conversions) %>% 
  #Add size category column 
  mutate(size_f = factor(case_when(
    size.length.um.used.for.conversions < 0.1 ~ "1nm < 100nm",
    size.length.um.used.for.conversions >= 0.1 & size.length.um.used.for.conversions < 1 ~ "100nm < 1µm",
    size.length.um.used.for.conversions >= 1 & size.length.um.used.for.conversions < 100 ~ "1µm < 100µm",
    size.length.um.used.for.conversions >= 100 & size.length.um.used.for.conversions < 1000 ~ "100µm < 1mm",
    size.length.um.used.for.conversions >= 1000 & size.length.um.used.for.conversions <= 5000 ~ "1mm < 5mm",
    is.na(size.length.um.used.for.conversions) ~ "Not Reported"),  
    levels = c("1nm < 100nm", "100nm < 1µm", "1µm < 100µm", "100µm < 1mm", "1mm < 5mm", "Not Reported"))) %>% # creates new column with nicer names and order by size levels.
  relocate(size_f, .after = size.width.um.used.for.conversions) %>% 
  ### assign average length to width ratio based on compartment
  mutate(R.ave = case_when(env_f == "Marine" & `exposure.route` == "water" ~ R.ave.water.marine,
                           env_f == "Marine" & `exposure.route` == "sediment" ~ R.ave.sediment.marine,
                           env_f == "Freshwater" & `exposure.route` == "water" ~ R.ave.water.freshwater,
                           env_f == "Freshwater" & `exposure.route` == "sediment" ~ R.ave.sediment.freshwater,
                           T ~ NA)) %>% # if doesn't meet these conditions, annotate as NA
  
  #Calculate particle surface area
  mutate(particle.surface.area.um2 = case_when(shape_f == "Sphere" ~ 4*pi*((size.length.um.used.for.conversions/2)^2),
                                               shape_f == "Fiber" & is.na(size.width.um.used.for.conversions) ~ SAfnx_fiber(width = 15, length = size.length.um.used.for.conversions), #assum 15 um width (kooi et al 2021)
                                               shape_f == "Fiber" & !is.na(size.width.um.used.for.conversions) ~ SAfnx_fiber(width = size.width.um.used.for.conversions, length = size.length.um.used.for.conversions), #if width is known
                                               shape_f == "Fragment" ~ SAfnx(a = size.length.um.used.for.conversions,
                                                                           b = R.ave * size.length.um.used.for.conversions,
                                                                           c = R.ave * 0.67 * size.length.um.used.for.conversions))) %>%
  relocate(particle.surface.area.um2, .after = size_f) %>% 
  #Calculate particle volume
  mutate(particle.volume.um3 = case_when(shape_f == "Sphere" ~ (4/3)*pi*((size.length.um.used.for.conversions/2)^3),
                                         shape_f == "Fiber" & is.na(size.width.um.used.for.conversions) ~ volumefnx_fiber(width = 15, length = size.length.um.used.for.conversions), #assume 15 um as width (kooi et al 2021)
                                         shape_f == "Fiber" & !is.na(size.width.um.used.for.conversions) ~ volumefnx_fiber(width = size.width.um.used.for.conversions, length = size.length.um.used.for.conversions), #if width reported
                                         shape_f == "Fragment" ~ volumefnx(R = R.ave, L = size.length.um.used.for.conversions))) %>% 
  relocate(particle.volume.um3, .after = particle.surface.area.um2) %>% 
  #Calculate particle mass
  mutate(mass.per.particle.mg = (particle.volume.um3*density.g.cm3)*0.000000001) %>% 
  relocate(mass.per.particle.mg, .after = particle.volume.um3) %>%
  ####
  # #calculate dose metrics accordingly
  # mutate(dose.surface.area.um2.mL.master = particle.surface.area.um2 * dose.particles.mL.master) %>% 
  # mutate(particle.surface.area.um2.mg = particle.surface.area.um2 / mass.per.particle.mg) %>% 
  
  # create label for polydispersity
  mutate(polydispersity = case_when(
    is.na(size.length.min.mm.nominal|size.length.min.mm.measured) ~ "monodisperse",
    !is.na(size.length.min.mm.nominal|size.length.min.mm.measured) ~ "polydisperse")) %>% 
  
  ####prioritize measured parameters for conversions ###
  # minima
  mutate(size.length.min.um.used.for.conversions = case_when(
    is.na(size.length.min.mm.measured) ~ size.length.min.mm.nominal * 1000,
    !is.na(size.length.min.mm.measured) ~ size.length.min.mm.measured * 1000)) %>% 
  mutate(size.width.min.um.used.for.conversions = case_when(
    shape_f == "Sphere" ~ size.length.min.um.used.for.conversions, #all dims same
    shape_f == "Fiber" ~ R.ave * size.length.min.um.used.for.conversions, #median holds for all particles (Kooi et al 2021)
    shape_f == "Not Reported" ~ R.ave * size.length.min.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
    shape_f == "Fragment" ~ R.ave * size.length.min.um.used.for.conversions)) %>% # average width to length ratio in the marine environment (kooi et al 2021)
  mutate(size.height.min.um.used.for.conversions = case_when(
    shape_f == "Sphere" ~ size.length.min.um.used.for.conversions, #all dims same
    shape_f == "Not Reported" ~ R.ave * 0.67 * size.length.min.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
    shape_f == "Fiber" ~  R.ave * size.length.min.um.used.for.conversions, #height same as width for fibers
    shape_f == "Fragment" ~ R.ave * 0.67 * size.length.min.um.used.for.conversions)) %>% # average width to length ratio in the marine environment AND average height to width ratio (kooi et al 2021)
  # maxima
  mutate(size.length.max.um.used.for.conversions = case_when(
    is.na(size.length.max.mm.measured) ~ size.length.max.mm.nominal * 1000,
    !is.na(size.length.max.mm.measured) ~ size.length.max.mm.measured * 1000)) %>% 
  mutate(size.width.max.um.used.for.conversions = case_when(
    shape_f == "Sphere" ~ size.length.max.um.used.for.conversions, #all dims same
    shape_f == "Fiber" ~ R.ave * size.length.max.um.used.for.conversions, #median holds for all particles (Kooi et al 2021) #there are no fibers
    shape_f == "Not Reported" ~ R.ave * size.length.max.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
    shape_f == "Fragment" ~ R.ave * size.length.max.um.used.for.conversions)) %>% # average width to length ratio in the marine environment (kooi et al 2021)
  mutate(size.height.max.um.used.for.conversions = case_when(
    shape_f == "Sphere" ~ size.length.max.um.used.for.conversions, #all dims same
    shape_f == "Not Reported" ~ R.ave * 0.67 * size.length.max.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
    shape_f == "Fiber" ~ R.ave * size.length.max.um.used.for.conversions, #hieght same as width
    shape_f == "Fragment" ~ R.ave * 0.67 * size.length.max.um.used.for.conversions)) %>%  # average width to length ratio in the marine environment AND average height to width ratio (kooi et al 2021)
  
  #calculate minimum and maximum surface area for polydisperse particles
  mutate(particle.surface.area.um2.min = SAfnx(a = size.length.min.um.used.for.conversions,
                                               b = size.width.min.um.used.for.conversions,
                                               c = size.height.min.um.used.for.conversions)) %>%
  mutate(particle.surface.area.um2.max = SAfnx(a = size.length.max.um.used.for.conversions,
                                               b = size.width.max.um.used.for.conversions,
                                               c = size.height.max.um.used.for.conversions)) %>% 
  #calculate minimum and maximum volume for polydisperse particles
  mutate(particle.volume.um3.min = volumefnx_poly(length = size.length.min.um.used.for.conversions,
                                                  width =  size.width.min.um.used.for.conversions)) %>% 
  mutate(particle.volume.um3.max = volumefnx_poly(length = size.length.max.um.used.for.conversions,
                                                  width = size.width.max.um.used.for.conversions)) %>% 
  #calculate minimum and maximum volume for polydisperse particles
  mutate(mass.per.particle.mg.min = massfnx_poly(length = size.length.min.um.used.for.conversions,
                                                 width = size.width.min.um.used.for.conversions,
                                                 p = density.g.cm3)) %>% #equation uses g/cm3
  mutate(mass.per.particle.mg.max = massfnx_poly(length = size.length.max.um.used.for.conversions,
                                                 width = size.width.max.um.used.for.conversions,
                                                 p = density.g.cm3)) %>%   #equation uses g/cm3
  
  #Mass (converted)
  mutate(dose.mg.L.master = ifelse(is.na(dose.mg.L.master), (dose.particles.mL.master*1000)*mass.per.particle.mg, dose.mg.L.master)) %>% 
  mutate(dose.mg.L.master.converted.reported = factor(ifelse((!is.na(dose.mg.L.master)&is.na(dose.mg.L.master.converted.reported)), "converted", dose.mg.L.master.converted.reported))) %>% 
      
  mutate(dose.mg.kg.sediment.master = ifelse(is.na(dose.mg.kg.sediment.master), (dose.particles.kg.sediment.master)*mass.per.particle.mg, dose.mg.kg.sediment.master)) %>% 
  mutate(dose.mg.kg.sediment.master.converted.reported = factor(ifelse((!is.na(dose.mg.kg.sediment.master)&is.na(dose.mg.kg.sediment.master.converted.reported)), "converted", dose.mg.kg.sediment.master.converted.reported))) %>% 
  
  #Count (converted)
  mutate(dose.particles.mL.master = ifelse(is.na(dose.particles.mL.master),(dose.mg.L.master/mass.per.particle.mg)/1000,dose.particles.mL.master)) %>% 
  mutate(dose.particles.mL.master.converted.reported = factor(ifelse((!is.na(dose.particles.mL.master)&is.na(dose.particles.mL.master.converted.reported)), "converted", dose.particles.mL.master.converted.reported))) %>% 

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

  #Factor particle weathering
  mutate(weathered.or.biofouled. = factor(weathered.or.biofouled.)) %>% 
  #Rename columns to match/look nicer
  rename(weather.biofoul_f = weathered.or.biofouled.,
         size.valid = size.validated.,
         shape.valid = shape.validated.,
         polymer.valid = polymer.validated.,
         sodium.azide = sodium.azide.present.,
         contaminant.screen = screened.for.chemical.contamination.,
         clean.method = particle.cleaning.,
         sol.rinse = solvent.rinse,
         background.plastics = background.contamination.monitored.,
         concentration.valid = concentration.validated.,
         uptake.valid = uptake.validated.,
         uptake.valid.method = uptake.validation.method,
         fed = organisms.fed.,
         `Nominal Dose Alternative Category` = `nominal.dose.alternative.category`,
         `Nominal Dose Alternative Type` = `nominal.dose...alternative.type`,
         `Nominal Dose Alternative Type Units` = `nominal.dose...alternative.type.units`,
         `Measured Dose Alternative Category` = `measured.dose.alternative.category`,
         `Measured Dose Alternative Type` = `measured.dose.alternative`,
         `Measured Dose Alternative Type Units` = `measured.dose.alternative.units`) 
  
  
#Add-in body size metrics from ToMEx 1.0

#Create summary data frame from ToMEx 1.0
  bodysize_summary <- aoc_setup %>%
    group_by(species_f, body.length.cm, body.size.source, max.size.ingest.mm, max.size.ingest.um) %>%
    summarise(.groups = "drop") 

  # bodysize_addons <- read_csv("scripts/monte carlo/ref data/gape_size.csv",
  #                             show_col_types = FALSE) #copied from aq_mp_tox_shiny main folder
    
  bodysize_addons <- bodysize_addons %>% 
    mutate(species_f = as.factor(species_f)) %>% 
    #annotate whether max size ingest was estimated or reported (all estiamted here)
    mutate(max.size.ingest.reported.estimated = "estimated",
           #calculate maximum ingestible size (if not already in database)
           max.size.ingest.mm = 10^(beta_log10_body_length * log10(body.length.cm * 10) - body_length_intercept),#(Jâms, et al 2020 Nature paper)correction for cm to mm
           max.size.ingest.um = 1000 * max.size.ingest.mm)

  bodysize_summary <- bind_rows(bodysize_summary,bodysize_addons)
  
#Join summary to tidy ToMEx 2.0 data frame

tomex2.0_aoc_setup <- left_join(tomex2.0_aoc_setup, bodysize_summary, by = c("species_f"),
                                relationship = "many-to-many")

#Re-structure alternative dosing columns in aoc-setup
aoc_setup <- aoc_setup %>% 
  #Remove effect metrics when there are less than 3 treatments
  mutate(effect.metric = as.character(effect.metric)) %>% 
  #Remove effect metrics when there are less than 3 treatments
  mutate(effect.metric = factor(ifelse(treatments < 3, NA_character_, effect.metric))) %>% 
  #Nominal Alternative Doses
  # mutate(dose.mg.kg.sed.nominal = nominal.dose.mg.kg.sediment,
  #        dose.particles.kg.sed.nominal = nominal.dose.particles.kg.sediment,
  #        dose.mg.kg.sed.measured = measured.dose.mg.kg.sediment,
  #        dose.particles.kg.sed.measured = measured.dose.particles.kg.sediment
  #        ) %>% 
  mutate(`Nominal Dose Alternative Type` = case_when(
    # !is.na(dose.mg.kg.sed.nominal) ~ dose.mg.kg.sed.nominal,
    !is.na(dose.mg.kg.food.nominal) ~ dose.mg.kg.food.nominal,
    !is.na(dose.particles.kg.food.nominal) ~ dose.particles.kg.food.nominal,
    !is.na(dose.percent.sed.nominal) ~ dose.percent.sed.nominal,
    !is.na(dose.particles.m2.nominal) ~ dose.particles.m2.nominal,
    !is.na(dose.percent.food.nominal) ~ dose.percent.food.nominal,
    # !is.na(dose.particles.kg.sed.nominal) ~ dose.particles.kg.sed.nominal,
  )) %>%  
  mutate(`Nominal Dose Alternative Type Units` = case_when(
    # !is.na(dose.mg.kg.sed.nominal) ~ "mg/kg sediment",
    !is.na(dose.mg.kg.food.nominal) ~ "mg/kg food",
    !is.na(dose.particles.kg.food.nominal) ~ "particles/kg food",
    !is.na(dose.percent.sed.nominal) ~ "percent mass (%) sediment",
    !is.na(dose.particles.m2.nominal) ~ "particles/m^2",
    !is.na(dose.percent.food.nominal) ~ "percent mass (%) food",
    # !is.na(dose.particles.kg.sed.nominal) ~ "particles/kg sediment",
  )) %>% 
  #Measured Alternative Doses
  mutate(`Measured Dose Alternative Type` = case_when(
    # !is.na(dose.mg.kg.sed.measured) ~ dose.mg.kg.sed.measured,
    !is.na(dose.mg.kg.food.measured) ~ dose.mg.kg.food.measured,
    !is.na(dose.particles.kg.food.measured) ~ dose.particles.kg.food.measured,
    !is.na(dose.particles.kg.food.min.measured) ~ dose.particles.kg.food.min.measured,
    !is.na(dose.particles.kg.food.max.measured) ~ dose.particles.kg.food.max.measured,
    !is.na(dose.percent.sed.measured) ~ dose.percent.sed.measured,
    !is.na(dose.particles.m2.measured) ~ dose.particles.m2.measured,
  )) %>%  
  mutate(`Measured Dose Alternative Type Units` = case_when(
    # !is.na(dose.mg.kg.sed.measured) ~ "mg/kg sediment",
    !is.na(dose.mg.kg.food.measured) ~ "mg/kg food",
    !is.na(dose.particles.kg.food.measured) ~ "particles/kg food",
    !is.na(dose.particles.kg.food.min.measured) ~ "particles/kg food minimum",
    !is.na(dose.particles.kg.food.max.measured) ~ "particles/kg food maximum",
    !is.na(dose.percent.sed.measured) ~ "percent mass (%) sediment",
    !is.na(dose.particles.m2.measured) ~ "particles/m^2",
  )) %>% 
  rename(nominal.dose.mg.kg.sediment = dose.mg.kg.sed.nominal,
         measured.dose.mg.kg.sediment = dose.mg.kg.sed.measured,
         nominal.dose.particles.kg.sediment = dose.particles.kg.sed.nominal)
  
#Join tomex2.0_aoc_setup to aoc_setup (from ToMEx 1.0)

#Get names of relevant columns from tomex 2.0
names <- tomex2.0_aoc_setup %>%
  select(-c(nominal.dose...mass, nominal.dose...mass.units, nominal.dose...particles, nominal.dose...particles.units, 
            nominal.added.chemical.dose, nominal.added.chemical.dose.units, 
            measured.dose...mass, measured.dose...mass.units, 
            measured.dose...particles, measured.dose...particles.units, 
            measured.chemical.dose, measured.chemical.dose.units)) %>% 
  colnames() 

#Select columns from each data frame
tomex2.0_aoc_setup <- tomex2.0_aoc_setup %>%
  select(all_of(names))

aoc_setup <- aoc_setup %>%
  mutate(sample.size = as.character(sample.size)) %>% 
  mutate(replicates = as.character(replicates)) %>% 
  mutate(media.sal.ppt = as.character(media.sal.ppt)) %>% 
  mutate(media.ph = as.character(media.ph)) %>% 
  mutate(media.temp = as.character(media.temp)) %>% 
  mutate(media.temp.min = as.character(media.temp.min)) %>% 
  mutate(media.temp.max = as.character(media.temp.max)) %>% 
  mutate(zetapotential.mV = as.character(zetapotential.mV)) %>% 
  mutate(chem.add.dose.mg.L.measured = suppressWarnings(as.numeric(aoc_setup$chem.add.dose.mg.L.measured))) %>% #warning message expected
  add_column(chem.add.measured = NA_character_,
             `Nominal Dose Alternative Category` = NA_character_,
             `Measured Dose Alternative Category` = NA_character_,
             `Recovery (Days)` = NA_real_,
             measured.dose.particles.kg.sediment = NA_real_,)

#remove unwanted character strings from DOI column
aoc_setup$doi <- gsub('https://dx.doi.org/','',aoc_setup$doi)
aoc_setup$doi <- gsub('https://doi.org/','',aoc_setup$doi)
aoc_setup$doi <- gsub('doi.org/','',aoc_setup$doi)
aoc_setup$doi <- gsub('https://','',aoc_setup$doi)

aoc_setup <- aoc_setup %>%
  select(all_of(names)) #If this throws an error, aoc_setup doesn't have all the column names as tomex2.0 df

#Join rows
tomex2.0_aoc_setup_final <- bind_rows(aoc_setup, tomex2.0_aoc_setup)

tomex2.0_aoc_setup_final %>% 
  mutate(source = as.factor(source))

##### Variable Management - Ana's code ####

#life stage
tomex2.0_aoc_setup_final$life_f =
  fct_collapse(tomex2.0_aoc_setup_final$life_f,
               "adult" = "Adult",
               "early" = c("Early", "Early, F0 generation", "Early, F1 generation", "Early, F2 generation", "Instar II", "Instar III"),
               "Not Reported" = c("early/juvenile/adult", "Not Reported"),
               "juvenile" = "Juvenile")

tomex2.0_aoc_setup_final$life_f = fct_na_value_to_level(tomex2.0_aoc_setup_final$life_f, "Not Reported")

#sex
tomex2.0_aoc_setup_final$sex = fct_collapse(tomex2.0_aoc_setup_final$sex,
                       "F" = "f",
                       "M" = "m",
                       "M_F" = c("M,F",  "M/F"))

tomex2.0_aoc_setup_final$sex = fct_na_value_to_level(tomex2.0_aoc_setup_final$sex, "Not Reported")

#exposure route
tomex2.0_aoc_setup_final$exposure.route = fct_collapse(tomex2.0_aoc_setup_final$exposure.route,
                                  "parental_water" = "coparental exposure+water",
                                  "water_food" = "water+food",
                                  "parental" =  c("Co-Parental Exposure", "coparental exposure", "maternal exposure","paternal exposure"))

tomex2.0_aoc_setup_final$exposure.route = fct_na_value_to_level(tomex2.0_aoc_setup_final$exposure.route, "Not Reported")

#particle mix
tomex2.0_aoc_setup_final$mix = fct_collapse(tomex2.0_aoc_setup_final$mix,
                                 "No"= c("N", "No"),
                                 "Yes"= c("Y", "Yes"))

#new column - dissolved organic matter
tomex2.0_aoc_setup_final$DOM_present = NA

tomex2.0_aoc_setup_final = tomex2.0_aoc_setup_final %>% group_by(exposure.media) %>%
  mutate(DOM_present = case_when(
    exposure.media %in% c("brackish water", "filtered aerated natural seawater", "filtered autoclaved natural seawater",  "filtered natural seawater",
                          "filtered seawater", "filtered sterilized seawater", "filtered uv treated seawater", "freshwater", "Hemolymph serum (mussel)",
                          "high nitrogen and phosphorus lake water", "marine sediment seawater", "natural brackish water", "natural river water",
                          "natural seawater", "seawater", "Seawater", "sediment", "sediment filtered seawater", "sediment freshwater",
                          "sediment seawater", "sterile natural seawater", "uv treated filtered seawater", "uv treated seawater") ~ "Yes",
    chem.add.nominal %in% c("harbor effluent", "seawater", "sewage") ~ "Yes")) %>% 
  ungroup()
   
tomex2.0_aoc_setup_final$DOM_present = fct_na_value_to_level(tomex2.0_aoc_setup_final$DOM_present, "No")

#exposure media
tomex2.0_aoc_setup_final$exposure.media = fct_collapse(tomex2.0_aoc_setup_final$exposure.media,
                                  "artificial_medium" = c("0.3x Danieau's solution", "10% Hoagland Solution", "30percent danieaus medium",
                                                          "AAP medium", "adam medium", "alga gro freshwater medium", "artificial freshwater",
                                                          "artificial medium", "artificial seawater", "basal medium", "bg11 medium", 
                                                          "BG11 medium", "Bold's Basal Medium", "culture medium", "cyanobacterial growth medium - AA/8+N medium",
                                                          "distilled water containing 10% Hoagland solution", "dulbecco's modified eagle's medium", "e3 embryo medium",
                                                          "egg water", "Elendt m4 medium", "Elendt M4 medium", "elendt m7 medium", "embryo culture medium", "f2 medium",
                                                          "filter-sterilized artificial seawater", "filtered artificial seawater", "filtered artificial water", "hanks balanced salt soution",
                                                          "hutners medium", "Immune cell exposure media", "instant ocean", "iso medium", "jarworski media", "M4 medium", "m7 medium",
                                                          "M7 medium", "oecd 202 culture media", "oecd media", "osmosis water", "reverse osmosis water", "RT medium", "SE medium",
                                                          "steinburg medium", "sterilized artificial seawater", "synthetic freshwater", "Synthetic seawater", "tris acetate phosphate medium",
                                                          "trisacetatephosphatemedium", "trout culture medium (TCM)", "WC medium", "zebrafish water", "AA/8+N"),
                                  "water" = c("3 parts tap 1 part softened water", "astm hard water", "circlulated water", "copper free water", "culture water",
                                              "Culture water", "dechlorinated tap water", "dechlorinated water", "deionized water", "diluted mineral water",
                                              "double distilled water", "filtered dechlorinated water", "hard reconstituted water", "moderately hard water",
                                              "Naturally dechlorinated water", "sediment dutch standard water", "standard salt water", "tap deionized water",
                                              "tap water", "uv sterilized water", "uv treated dechlorinated tap water", "uv treated water", "uv-sterilized milli-q water",
                                              "uv-sterilized well aerated dechlorinated tap water"),
                                  "freshwater"= c("culture fresh water", "high nitrogen and phosphorus lake water", "natural river water",
                                                  "sediment freshwater"),
                                  "other" =  c("brackish water", "Hemolymph serum (mussel)", "natural brackish water"),
                                  "seawater" = c("filtered aerated natural seawater", "filtered autoclaved natural seawater", "filtered natural seawater", "filtered seawater",
                                                 "filtered sterilized seawater", "filtered uv treated seawater", "marine sediment seawater", "natural seawater", "seawater", "Seawater",
                                                 "sediment artificial seawater", "sediment filtered seawater", "sediment seawater", "sterile natural seawater", "uv treated filtered seawater",
                                                 "uv treated seawater"))
tomex2.0_aoc_setup_final$exposure.media = fct_na_value_to_level(tomex2.0_aoc_setup_final$exposure.media, "Not Reported")

#solvent
tomex2.0_aoc_setup_final$solvent = fct_na_value_to_level(tomex2.0_aoc_setup_final$solvent, "Not Reported")

#detergent
tomex2.0_aoc_setup_final$detergent = fct_collapse(tomex2.0_aoc_setup_final$detergent,
                             "cetyl_alcohol" = "cetyl alcohol",
                             "SDS"= c("SDS", "unknown, but similar to sodium dodecyl sulfate (SDS)"),
                             "triton_x-100" = "Triton X-100", 
                             "tween80"= c("Tween80", "tween80"))

tomex2.0_aoc_setup_final$detergent = fct_na_value_to_level(tomex2.0_aoc_setup_final$detergent, "Not Reported")

#chemicals added
tomex2.0_aoc_setup_final$chem.add.nominal = fct_collapse(tomex2.0_aoc_setup_final$chem.add.nominal,
                                   "silver_nitrate" = "AgNO3",
                                   "silver" = "AgNP",
                                   "bap" = "BaP",
                                   "benzophenone" = c("BP3", "BP4", "benzophenone-3"),
                                   "cadmium" = "Cadmium", 
                                   "chlorpyrifos" = "Chlorpyrifos", 
                                   "glyphosate_formulation" = "glyphosate-based herbicide (Roundup)",
                                   "phenantherene" = "phenanthrene",
                                   "esfenvalerate" = "pyrethroid pesticide (esfenvalerate)",
                                   "mix" = "various - see paper")

tomex2.0_aoc_setup_final$chem.add.nominal = fct_na_value_to_level(tomex2.0_aoc_setup_final$chem.add.nominal, "Not Reported")

#direction
tomex2.0_aoc_setup_final$direction = fct_collapse(tomex2.0_aoc_setup_final$direction,
                             "decrease"= "Decreasing",
                             "increase" = c("Increase", "Increasing"),
                             "increase.decrease" = c("increase and decrease", "increase.decrease", "increase/decrease"),
                             "Not Reported" = c("Na"))

tomex2.0_aoc_setup_final$direction = fct_na_value_to_level(tomex2.0_aoc_setup_final$direction, "Not Reported")

#target cell or tissue
tomex2.0_aoc_setup_final$target.cell.tissue  = fct_collapse(tomex2.0_aoc_setup_final$target.cell.tissue,
                                          "digestive.gland" = "digestive.gland:total",
                                          "gills" = c("gill", "Gills"),
                                          "intestine" = c("digestive.tract", "distal intestine", "gut", "intestines", "proximal intestine"),
                                          "whole.body" = c("dstage.larvae", "embryo", "gut.whole.body"),
                                          "kidney" = "kidney.neutrophil",
                                          "liver" = "Liver")

tomex2.0_aoc_setup_final$target.cell.tissue = fct_na_value_to_level(tomex2.0_aoc_setup_final$target.cell.tissue , "Not Reported")

#polymer
tomex2.0_aoc_setup_final$poly_f = fct_collapse(tomex2.0_aoc_setup_final$poly_f,
"Polyamide" = c("Polyamide", "Polyamide 66"),
"Polyethylene vinyl acetate" = c("Polyethylene (co-Vinyl acetate)", "Polyethylene Vinyl Acetate"),
"Not Reported" = "Not Reported",
"Polystyrene" = "Latex")

tomex2.0_aoc_setup_final$poly_f = fct_na_value_to_level(tomex2.0_aoc_setup_final$poly_f, "Not Reported")

# exclude Polyamidoamine -> not MP
tomex2.0_aoc_setup_final <- tomex2.0_aoc_setup_final[tomex2.0_aoc_setup_final$poly_f != "Polyamidoamine", ]

#particle source
tomex2.0_aoc_setup_final$particle.source = fct_collapse(tomex2.0_aoc_setup_final$particle.source,
                                   "environmental" = c("Environmental",  "Field collected"),
                                   "not_reported" = c("N",  "No"),
                                   "commercial" = "Commercial") 

#sodium azide
tomex2.0_aoc_setup_final$sodium.azide = fct_collapse(tomex2.0_aoc_setup_final$sodium.azide,
                                         "No" =  c("N",  "No"),
                                         "Yes" = c("Y",  "Yes"))

#cleaning
tomex2.0_aoc_setup_final$clean.method = fct_collapse(tomex2.0_aoc_setup_final$clean.method,
                                      "not_cleaned" = c("N",  "No", "Not Cleaned"),
                                      "centrifugal_filter" = c("centrifugal filter",  "dialysis"),
                                      "rinse" = c("centrifugation rinse",  "rinse", "rinse dialyzed"),
                                      "ultrasonic" = "ultrasonic extraction")

#solvent
tomex2.0_aoc_setup_final$sol.rinse = fct_collapse(tomex2.0_aoc_setup_final$sol.rinse,
                                 "none" = c("N",  "No"),
                                 "distilled water" = "distillled water", 
                                 "ultrapure water" = "ultrapurewater")

tomex2.0_aoc_setup_final$sol.rinse = fct_na_value_to_level(tomex2.0_aoc_setup_final$sol.rinse, "Not Reported")

#particle behavior
tomex2.0_aoc_setup_final$particle.behavior = fct_collapse(tomex2.0_aoc_setup_final$particle.behavior,
                                     "none" = c("aeration",  "aeration,pdi quantified", "air pump", "beads did not float air pump in tanks",
                                                "cultures agitated", "detergent", "dynamic light scatter", "dynamic light scatter aerated",
                                                "gentle resuspension", "homogenized cement mill", "homogenized paintmixer", "images captured in test medium",
                                                "microscopic observation", "mixed thoroughly in bulk sediment", "mixing testing solution twice daily",
                                                "particles conditioned in a solution of Triton X-100", "particles conditioned in river water",
                                                "particles in food pellets", "particles mixed into sediment", "particles mixed into top layer of sediment",
                                                "particles weathered to ensure suspension in water column", "polydisperisty index quantified",
                                                "quantified particles on sediment surface", "resuspended daily", "sedimentation rate quantified",
                                                "Sedimentation rate quantified", "shaken by hand twice a day", "shaker", "shaker,mechnical rotator",
                                                "sonicated", "sonicated aeration", "sonication aeration", "stirring", "stock solution on rolling mixer",
                                                "surfactant added", "vortex sonicated", "vortexed dynamic light scatter", "water circulation",
                                                "water column concentrations monitored"),
                                     "adhering_to_surface" = "adhere surface chorion",
                                     "not_settling" = "aeration no settling observed",
                                     "aggregates" = c("aggregates observed by dynamic light scattering",  "aggregates observed via dynamic light scattering", "aggregates observed via TEM",
                                                      "aggregates observed visually", "poly disperisity quantified aggregates observed", "SEM aggregates observed", "some aggregates",
                                                      "sonicated aggregates observed", "sonication aeration aggregates observed", "TEM aggregates"),
                                     "uniformly_dispersed" = c("confirmed homogeneity in food microscopy", "homogenous mixture within food microscopy", "methodunclear noaggregation",
                                                               "microscopic observation no aggregation", "particles described as neutrally bouyant dispersion quantified", "particles mixed into food confirmed by microscopy",
                                                               "sonicated no aggregates observed", "sonicated no aggregation observed via dynamic light scattering", "sonication no aggregates observed",
                                                               "suspended in water column", "visual observation no aggregates"),
                                     "floating" = "floated on medium surface",
                                     "Not Evaluated" = c("No", "not evaluated", "Not Evaluated"),
                                     "settling" = c("particles allowed to sink", "quantified settling", "settled on bottom", "sink to bottom"))

tomex2.0_aoc_setup_final$particle.behavior = fct_na_value_to_level(tomex2.0_aoc_setup_final$particle.behavior, "Not Evaluated")

#organisms fed
tomex2.0_aoc_setup_final$fed = fct_collapse(tomex2.0_aoc_setup_final$fed,
                                  "Yes" = c("0.5 mg CL-1 M. minutum combined w/MPs",  "0.5 mg CL-1 mix of M. minutum and Cryptomonas sp. Combined w/MPs",
                                            "16 blue mussels", "2 blue mussels", "2 mg CL-1 M. minutum combined w/MPs",
                                            "2 mg CL-1 mix of M. minutum and Cryptomonas sp. Combined w/MPs", "3*10^4 cells T. suecica",
                                            "32 blue mussels", "4 blue mussels", "6*10^4 cells T. suecica", "6*10^5 cells T. suecica",
                                            "64 blue mussels", "8 blue mussels", "High food - food content 1.2% dry sediment weight",
                                            "Low food - food content 0.15% dry sediment weight", "Y", "Yes"),
                                  "No" = c("N",  "No"))

tomex2.0_aoc_setup_final$fed = fct_na_value_to_level(tomex2.0_aoc_setup_final$fed, "No")

#weathering
tomex2.0_aoc_setup_final$weather.biofoul_f = fct_na_value_to_level(tomex2.0_aoc_setup_final$weather.biofoul_f, "No")

#add 0s for chemical dose
tomex2.0_aoc_setup_final$chem.add.dose.mg.L.nominal = ifelse(is.na(tomex2.0_aoc_setup_final$chem.add.dose.mg.L.nominal), 0, tomex2.0_aoc_setup_final$chem.add.dose.mg.L.nominal)

#Fill width with length for rows where shape is "sphere"
tomex2.0_aoc_setup_final$size.width.um.used.for.conversions[tomex2.0_aoc_setup_final$shape_f == "Sphere" & !is.na(tomex2.0_aoc_setup_final$size.length.um.used.for.conversions)] = 
  tomex2.0_aoc_setup_final$size.length.um.used.for.conversions[tomex2.0_aoc_setup_final$shape_f == "Sphere" & !is.na(tomex2.0_aoc_setup_final$size.length.um.used.for.conversions)]

#Add 1 as sample size for Algae, Bacterium, Cyanobacteria, Plant, and Mixed
tomex2.0_aoc_setup_final <- tomex2.0_aoc_setup_final %>%
  mutate(sample.size = ifelse(org_f %in% c("Algae", "Bacterium", "Cyanobacteria", "Dinoflagellate"), 1, sample.size))

#Calculate the median density for each polymer
median_density_by_polymer = tomex2.0_aoc_setup_final %>%
  group_by(poly_f) %>%
  summarize(median_density = median(density.g.cm3, na.rm = TRUE)) %>% 
  ungroup()

#Join the mean temperature back to the original data
tomex2.0_aoc_setup_final = left_join(tomex2.0_aoc_setup_final, median_density_by_polymer, by = "poly_f")

#Fill NA values of density.g.cm.3. with the calculated median per polymer type
tomex2.0_aoc_setup_final = tomex2.0_aoc_setup_final %>%
  mutate(density.g.cm3 = ifelse(is.na(density.g.cm3), median_density, density.g.cm3))

#Remove the temporary column median_density
tomex2.0_aoc_setup_final = tomex2.0_aoc_setup_final %>% select(-median_density)

#Identify unique Polymers with NA values in density.g.cm3
na_polymer = unique(tomex2.0_aoc_setup_final$poly_f[is.na(tomex2.0_aoc_setup_final$density.g.cm3)])

# Convert to character and sort the unique species alphabetically
na_polymer = data.frame(poly_f = sort(as.character(na_polymer)))
na_polymer # >>> check the density in the literature

#
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "High Density Polyethylene"] = "0.953"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "Poly(Styrene-co-acrylonitrile)"] = "1.08"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "polyethylene_vinyl_acetate"] = "0.953"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "Polyisoprene"] = "0.905"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "Polytetrafluoroethylene"] = "2.2"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "Polyurethane"] = "1.19"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "Polyvinyl Acetate"] = "1.19"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "Polyvinylchloride/vinylacetate co-polymer"] = "1.38"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "Sodium Polyacrylate"] = "1.32"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "Starch/Polybutylene Adipate Terephthalate/Polylactic Acid"] = "1.05"
tomex2.0_aoc_setup_final$density.g.cm3[tomex2.0_aoc_setup_final$poly_f == "Tire Wear"] = "1.45"

# Restructure temp column
tomex2.0_aoc_setup_final$media.temp = as.numeric(tomex2.0_aoc_setup_final$media.temp)

# Calculate the median temperature for each Species
median_temp_by_group = tomex2.0_aoc_setup_final %>%
  group_by(species_f) %>%
  summarize(median_temp = median(media.temp, na.rm = TRUE)) %>% 
  ungroup()

# Join the mean temperature back to the original data
tomex2.0_aoc_setup_final = left_join(tomex2.0_aoc_setup_final, median_temp_by_group, by = "species_f")

# Fill NA values of media.temp with the calculated median temperature
tomex2.0_aoc_setup_final = tomex2.0_aoc_setup_final %>%
  mutate(media.temp = ifelse(is.na(media.temp), median_temp, media.temp))

# Remove the temporary column median_temp
tomex2.0_aoc_setup_final = tomex2.0_aoc_setup_final %>% select(-median_temp)

# Calculate the mean temperature for each organism group
median_temp_by_group = tomex2.0_aoc_setup_final %>%
  group_by(org_f) %>%
  summarize(median_temp = median(media.temp, na.rm = TRUE)) %>% 
  ungroup()

# Join the mean temperature back to the original data
tomex2.0_aoc_setup_final <- left_join(tomex2.0_aoc_setup_final, median_temp_by_group, by = "org_f")

# Fill NA values of media.temp with the calculated mean temperature
tomex2.0_aoc_setup_final <- tomex2.0_aoc_setup_final %>%
  mutate(media.temp = ifelse(is.na(media.temp), median_temp, media.temp))

# Populate NA values for quality scores where studies are not "particle only"
tomex2.0_aoc_setup_final_test <- tomex2.0_aoc_setup_final %>% 
  mutate(technical.quality = case_when(exp_type_f == "Particle Only" ~ technical.quality,
                                       exp_type_f == "Leachate" ~ NA_real_,
                                       exp_type_f == "Chemical Co-Exposure" ~ NA_real_,
                                       exp_type_f == "Chemical Transfer" ~ NA_real_)) %>% 
  mutate(risk.quality = case_when(exp_type_f == "Particle Only" ~ risk.quality,
                                  exp_type_f == "Leachate" ~ NA_real_,
                                  exp_type_f == "Chemical Co-Exposure" ~ NA_real_,
                                  exp_type_f == "Chemical Transfer" ~ NA_real_)) %>%
  mutate(total.quality = case_when(exp_type_f == "Particle Only" ~ total.quality,
                                   exp_type_f == "Leachate" ~ NA_real_,
                                   exp_type_f == "Chemical Co-Exposure" ~ NA_real_,
                                   exp_type_f == "Chemical Transfer" ~ NA_real_))

#Save RDS file
#saveRDS(tomex2.0_aoc_setup_final, file = "aoc_setup_tomex2.RDS")

##### AOC ENDPOINT #####

#Endpoint Categorization setup
tomex2.0_aoc_endpoint_final <- tomex2.0_aoc_setup_final %>% 
  group_by(lvl1_f,lvl2_f,lvl3_f,bio_f) %>% 
  summarise(.groups = "drop")

#Save RDS file
#saveRDS(tomex2.0_aoc_endpoint_final, file = "aoc_endpoint_tomex2.RDS")

##### AOC SEARCH #####

tomex2.0_aoc_search <- tomex2.0_aoc_setup_final %>% 
  dplyr::select(source, doi, authors, year, tier_zero_tech_f, tier_zero_risk_f, species_f, org_f, env_f, life_f, vivo_f, sex, body.length.cm, max.size.ingest.mm,
                #experimental parameters
                exp_type_f, exposure.route, mix, negative.control, reference.material, exposure.media, solvent, detergent,
                media.ph, media.sal.ppt, media.temp, media.temp.min, media.temp.max, exposure.duration.d, `Recovery (Days)`, acute.chronic_f,
                treatments, replicates, sample.size, dosing.frequency, chem.add.nominal, chem.add.dose.mg.L.nominal, chem.add.dose.mg.L.measured,
                #master/alternative doses
                dose.particles.mL.master, dose.particles.mL.master.converted.reported, dose.mg.L.master, dose.mg.L.master.converted.reported,
                dose.um3.mL.master, dose.um2.mL.master, dose.um2.ug.mL.master, 
                
                dose.particles.kg.sediment.master, dose.particles.kg.sediment.master.converted.reported, dose.mg.kg.sediment.master, dose.mg.kg.sediment.master.converted.reported,
                dose.um3.kg.sediment.master, dose.um2.kg.sediment.master, dose.um2.ug.kg.sediment.master,
                
                `Nominal Dose Alternative Category`, `Nominal Dose Alternative Type`, `Nominal Dose Alternative Type Units`,
                `Measured Dose Alternative Category`, `Measured Dose Alternative Type`, `Measured Dose Alternative Type Units`,
                #biological effects
                effect_f, direction, lvl1_f, lvl2_f, lvl3_f, bio_f, target.cell.tissue, effect.metric, af.time, af.noec,
                #particle characteristics
                poly_f, shape_f, density.g.cm3, density.reported.estimated, charge, zetapotential.mV, zetapotential.media, functional.group,
                size.length.um.used.for.conversions, size.width.um.used.for.conversions, size_f, particle.surface.area.um2, particle.volume.um3,
                mass.per.particle.mg, weather.biofoul_f,
                #quality
                size.valid, polymer.valid, shape.valid, particle.source, sodium.azide, contaminant.screen, clean.method, sol.rinse, background.plastics,
                concentration.valid, particle.behavior, uptake.valid, uptake.valid.method, tissue.distribution, fed,
                #scores
                tech.a1, tech.a2, tech.a3, tech.a4, tech.a5, tech.a6, tech.1, tech.2, tech.3, tech.4, tech.5,
                tech.6, tech.7, tech.8, tech.9, tech.10, tech.11, tech.12, risk.b1, risk.13, risk.14, risk.15, risk.16, risk.17, risk.18, risk.19, risk.20,
                )

#Turn all character strings into factors if they aren't already so they are searchable via dropdown
tomex2.0_aoc_search[sapply(tomex2.0_aoc_search, is.character)] <- lapply(tomex2.0_aoc_search[sapply(tomex2.0_aoc_search, is.character)], as.factor)

tomex2.0_aoc_search_final <- tomex2.0_aoc_search %>% 
  dplyr::rename('Source' = source,'DOI' = doi,'Authors' = authors, 'Year' = year, 'Technical "Red Criteria"' = tier_zero_tech_f, 
                'Risk Assessment "Red Criteria"' = tier_zero_risk_f,'Species' = species_f, 
                'Organism Group' = org_f, 'Environment' = env_f, 'Life Stage' = life_f, 'In vitro/in vivo' = vivo_f,
                'Sex' = sex, 'Estimated Body Length (cm)' = body.length.cm, 
                'Estimated Maximum Ingestible Size (mm)' = max.size.ingest.mm, 'Experiment Type' = exp_type_f,
                'Exposure Route' = exposure.route, 'Particle Mix?' = mix, 'Negative Control' = negative.control, 
                'Reference Particle' = reference.material, 'Exposure Media' = exposure.media,
                'Solvent' = solvent, 'Detergent' = detergent, 'pH' = media.ph, 'Salinity (ppt)' = media.sal.ppt, 
                'Temperature (Avg)' = media.temp, 'Temperature (Min)' = media.temp.min,
                'Temperature (Max)' = media.temp.max, 'Exposure Duration (days)' = exposure.duration.d, 
                'Acute/Chronic' = acute.chronic_f, 'Number of Doses' = treatments, 'Replicates' = replicates,
                'Sample Size' = sample.size, 'Dosing Frequency' = dosing.frequency, 'Chemicals Added' = chem.add.nominal, 
                'Added Chemical Dose (mg/L, nominal)' = chem.add.dose.mg.L.nominal,
                'Added Chemical Dose (mg/L, measured)' = chem.add.dose.mg.L.measured, 
                'Effect' = effect_f, 'Direction' = direction, 'Broad Endpoint Category' = lvl1_f, 
                'Specific Endpoint Category' = lvl2_f,'Endpoint' = lvl3_f, 
                'Level of Biological Organization' = bio_f, 'Target Cell or Tissue' = target.cell.tissue, 
                'Effect Metric' = effect.metric, 'AF Time' = af.time,
                'AF NOEC' = af.noec, 'Polymer' = poly_f, 'Shape' = shape_f, 'Density (g/cm^3)' = density.g.cm3, 
                'Density, reported or estimated' = density.reported.estimated, 'Charge' = charge,
                'Zeta Potential (mV)' = zetapotential.mV, 'Size Category' = size_f,'Zeta Potential Media' = zetapotential.media, 
                'Functional Group' = functional.group, 'Particle Length (μm)' = size.length.um.used.for.conversions,
                'Particle Width (μm)' = size.width.um.used.for.conversions, 'Particle Surface Area (μm^2)' = particle.surface.area.um2, 
                'Particle Volume (μm^3)' = particle.volume.um3, 'Particle Mass (mg)'= mass.per.particle.mg,
                'Weathered or Biofouled?' = weather.biofoul_f, 'Size Validated?' = size.valid, 
                'Polymer Validated?' = polymer.valid, 'Shape Validated' = shape.valid, 'Particle Source' = particle.source,
                'Sodium Azide Present?' = sodium.azide,
                'Screened for Chemical Contamination?' = contaminant.screen, 'Particle Cleaning?' = clean.method, 
                'Solvent Rinse'= sol.rinse, 'Background Contamination Monitored?' = background.plastics,
                'Concentration Validated?' = concentration.valid, 
                'Particle Behavior' = particle.behavior, 'Uptake Validated?' = uptake.valid, 'Uptake Validation Method' = uptake.valid.method,
                'Tissue Distribution' = tissue.distribution, 'Organisms Fed?' = fed, 
                'Test Medium Score' = tech.a1, 'Administration Route Score' = tech.a2, 'Test Species Score'= tech.a3, 'Sample Size Score'= tech.a4, 
                'Control Group Score'= tech.a5, 'Exposure Duration Score'= tech.a6, 'Particle Size Score'= tech.1,
                'Particle Shape Score'= tech.2, 'Polymer Type Score'= tech.3, 'Source of Microplastics Score'= tech.4, 'Data Reporting Score'= tech.5, 
                'Chemical Purity Score'= tech.6, 'Laboratory Preparation Score'= tech.7, 'Background Contamination Score'= tech.8,
                'Exposure Verification Score'= tech.9, 'Exposure Homogeneity Score'= tech.10, 'Exposure Assessment Score'= tech.11, 
                'Replication Score'= tech.12, 'Number of Treatments Score'= risk.b1, 'Endpoints Score'= risk.13, 'Food Availability Score'= risk.14,
                'Effect Thresholds Score'= risk.15, 'Dose Response Score'= risk.16, 'Concentration Range Score'= risk.17, 'Aging and Biofouling Score'= risk.18, 
                'Microplastic Diversity Score' = risk.19, 'Exposure Time Score' = risk.20,
                #water
                "particles/mL water (master)" = dose.particles.mL.master, "particles/mL water (master), reported or converted" = dose.particles.mL.master.converted.reported,
                "μg/mL water (master)" = dose.mg.L.master, "μg/mL water (master), reported or converted" = dose.mg.L.master.converted.reported,
                "μm^3/mL water (master)" = dose.um3.mL.master, "μm^2/mL water (master)" = dose.um2.mL.master, "μm^2/ug/mL water (master)" = dose.um2.ug.mL.master,
                #sediment
                "particles/kg sediment dry weight (master)" = dose.particles.kg.sediment.master, "particles/kg sediment dry weight (master), reported or converted" = dose.particles.kg.sediment.master.converted.reported, 
                "mg/kg sediment dry weight (master)" = dose.mg.kg.sediment.master, "mg/kg sediment dry weight (master), reported or converted" = dose.mg.kg.sediment.master.converted.reported,
                "μm^3/kg sediment dry weight (master)" = dose.um3.kg.sediment.master, "μm^2/kg sediment dry weight (master)" = dose.um2.kg.sediment.master, "μm^2/ug/kg sediment dry weight (master)" = dose.um2.ug.kg.sediment.master)

#Save RDS file
#saveRDS(tomex2.0_aoc_search_final, file = "aoc_search_tomex2.RDS")

##### AOC SSD (AOC_Z) #####

tomex2.0_aoc_z_final <- tomex2.0_aoc_setup_final %>%
  rename(Species = species_f) %>%  #must make value 'Species" (uppercase)
  rename(Group = org_f)

#tomex2.0_aoc_z_final$rowid <- sapply(seq_len(nrow(tomex2.0_aoc_z_final)), function(x) uuid::UUIDgenerate())

#Save RDS file
#saveRDS(tomex2.0_aoc_z_final, file = "aoc_z_tomex2.RDS")

##### AOC QUALITY #####

tomex2.0_aoc_quality_final <- tomex2.0_aoc_setup_final %>%
  filter(tier_zero_tech_f != "Scoring Not Applicable") %>% 
  filter(tier_zero_risk_f != "Scoring Not Applicable") %>% 
  mutate(Study = paste0(authors, " (", year,")")) %>%
  mutate(Study_plus = as.factor(paste0(authors, " (", year,")", " (",doi,")"))) %>%
  group_by(Study, Study_plus, doi, treatment_range, tech.a1, tech.a2, tech.a3, tech.a4, tech.a5, tech.a6, tech.1, tech.2, tech.3, tech.4, tech.5,
           tech.6, tech.7, tech.8, tech.9, tech.10, tech.11, tech.12, risk.13, risk.14, risk.15, risk.16, risk.17, risk.18, risk.19, risk.20,
           lvl1_f, lvl2_f, bio_f, effect_f, life_f, poly_f, shape_f, size_f, species_f, env_f, org_f, acute.chronic_f, tier_zero_tech_f, tier_zero_risk_f) %>% 
  summarise(.groups = "drop")

tomex2.0_aoc_quality_final <- tomex2.0_aoc_quality_final %>% 
  pivot_longer(!c(Study, Study_plus, doi, treatment_range, lvl1_f, lvl2_f, bio_f, effect_f, life_f, poly_f, shape_f, size_f, species_f, env_f, org_f, acute.chronic_f, tier_zero_tech_f, tier_zero_risk_f),
               names_to ="Criteria", 
               values_to ="Score") %>%  
  #Assign descriptions to numerical scores
  mutate(Score_f = factor(case_when(Score == 0 ~ "Inadequate",
                                    Score == 1 ~ "Adequate with Restrictions",
                                    Score == 2 ~ "Adequate"))) %>% 
  #Assign each criteria to appropriate category
  mutate(Category = case_when(Criteria == "tech.a1" ~ "Technical",
                              Criteria == "tech.a2" ~ "Technical",
                              Criteria == "tech.a3" ~ "Technical",
                              Criteria == "tech.a4" ~ "Technical",
                              Criteria == "tech.a5" ~ "Technical",
                              Criteria == "tech.a6" ~ "Technical",
                              Criteria == "tech.1" ~ "Technical",
                              Criteria == "tech.2" ~ "Technical",
                              Criteria == "tech.3" ~ "Technical",
                              Criteria == "tech.4" ~ "Technical",
                              Criteria == "tech.5" ~ "Technical",
                              Criteria == "tech.6" ~ "Technical",
                              Criteria == "tech.7" ~ "Technical",
                              Criteria == "tech.8" ~ "Technical",
                              Criteria == "tech.9" ~ "Technical",
                              Criteria == "tech.10" ~ "Technical",
                              Criteria == "tech.11" ~ "Technical",
                              Criteria == "tech.12" ~ "Technical",
                              Criteria == "risk.13" ~ "Risk Assessment",
                              Criteria == "risk.13" ~ "Risk Assessment",
                              Criteria == "risk.14" ~ "Risk Assessment",
                              Criteria == "risk.15" ~ "Risk Assessment",
                              Criteria == "risk.16" ~ "Risk Assessment",
                              Criteria == "risk.17" ~ "Risk Assessment",
                              Criteria == "risk.18" ~ "Risk Assessment",
                              Criteria == "risk.19" ~ "Risk Assessment",
                              Criteria == "risk.20" ~ "Risk Assessment")) %>%
  #Set order of categories so they plot in correct order
  mutate(Category_f = factor(Category, levels = c("Technical", "Risk Assessment"))) %>% 
  #Assign descriptions to each criteria
  mutate(Criteria = case_when(Criteria == "tech.a1" ~ "Test Medium*",
                              Criteria == "tech.a2" ~ "Administration Route*",
                              Criteria == "tech.a3" ~ "Test Species*",
                              Criteria == "tech.a4" ~ "Sample Size*",
                              Criteria == "tech.a5" ~ "Control Group*",
                              Criteria == "tech.a6" ~ "Exposure Duration*",
                              Criteria == "tech.1" ~ "Particle Size*",
                              Criteria == "tech.2" ~ "Particle Shape*",
                              Criteria == "tech.3" ~ "Polymer Type*",
                              Criteria == "tech.4" ~ "Source of Microplastics*",
                              Criteria == "tech.5" ~ "Data Reporting*",
                              Criteria == "tech.6" ~ "Chemical Purity",
                              Criteria == "tech.7" ~ "Laboratory Preparation",
                              Criteria == "tech.8" ~ "Background Contamination",
                              Criteria == "tech.9" ~ "Exposure Verification",
                              Criteria == "tech.10" ~ "Exposure Homogeneity",
                              Criteria == "tech.11" ~ "Exposure Assessment",
                              Criteria == "tech.12" ~ "Replication",
                              Criteria == "risk.13" ~ "Endpoints",
                              Criteria == "risk.14" ~ "Food Availability",
                              Criteria == "risk.15" ~ "Effect Thresholds",
                              Criteria == "risk.16" ~ "Dose Response",
                              Criteria == "risk.17" ~ "Concentration Range",
                              Criteria == "risk.18" ~ "Aging and Biofouling",
                              Criteria == "risk.19" ~ "Microplastic Diversity",
                              Criteria == "risk.20" ~ "Exposure Time")) %>% 
  #Create factor for criteria and set order - need to be in reverse order here to plot correctly
  mutate(Criteria_f = factor(Criteria, levels = c("Exposure Time", "Microplastic Diversity", "Aging and Biofouling", "Concentration Range", "Dose Response",
                                                  "Effect Thresholds", "Food Availability", "Endpoints", "Replication", "Exposure Assessment", "Exposure Homogeneity",
                                                  "Exposure Verification", "Background Contamination", "Laboratory Preparation","Chemical Purity","Data Reporting*",
                                                  "Source of Microplastics*","Polymer Type*","Particle Shape*","Particle Size*","Exposure Duration*","Control Group*",
                                                  "Sample Size*", "Test Species*", "Administration Route*","Test Medium*")))

#Save RDS file
#saveRDS(tomex2.0_aoc_quality_final, file = "aoc_quality_tomex2.RDS")

return(tomex2.0_aoc_z_final)
} #end function (DO NOT DELETE)


#quick stats
# library(tidyverse)
# 
# pubs <- as.data.frame(unique(aoc_setup_tomex2$doi))
# 
# study_types <- aoc_setup_tomex2 %>%
#   group_by(doi, exp_type_f) %>% 
#   summarise()
#   
# acute_chronic <- aoc_setup_tomex2 %>% 
#   filter(acute.chronic_f == "Chronic")
#   
# vivo <- aoc_setup_tomex2 %>% 
#   filter(vivo_f == "In Vitro")
# 
# species <- aoc_setup_tomex2 %>%
#   filter(env_f == "Freshwater") %>% 
#   group_by(species_f) %>% 
#   summarise()
# 
# effect_metrics <- aoc_setup_tomex2 %>% 
#   # filter(source != "ToMEx 2.0") %>% 
#   filter(effect.metric %in% c("EC50", "LC50", "EC10", "IC50", "EC20", "LC20")) %>% 
#   group_by(doi, env_f, org_f, species_f, effect.metric, lvl1_f, lvl2_f, lvl3_f) %>% 
#   summarise()

# not_tidy <- tomex2.0_aoc_setup_final %>% 
#   filter(!is.na(`Issue Flag`)) #857 lines
# 
# not_tidy_studies <- not_tidy %>% 
#   group_by(doi,authors,`Issue Flag`) %>% 
#   summarise() #24 studies with issue flags
