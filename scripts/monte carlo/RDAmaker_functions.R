##### RDA MAKER ####
# preps dataframe for use in RShiny app #

library(tidyverse) #General everything

source("scripts/functions.R") # necessary for surface area, volume calculations

# R.ave.water.marine <- 0.77 # average length to width ratio of microplastics in marine environment (Kooi et al. 2021)
# R.ave.water.freshwater <- 0.67
# R.ave.sediment.marine <- 0.75
# R.ave.sediment.freshwater <- 0.70


ToMEx1.0fxn <- function(R.ave.water.marine, R.ave.water.freshwater, R.ave.sediment.freshwater, R.ave.sediment.marine,
                        beta_log10_body_length, body_length_intercept){
##### Read in Data ####
#aoc <- read_csv("AquaticOrganisms_Clean_final.csv", guess_max = 10000) %>% rowid_to_column()
aoc <- read_csv("scripts/monte carlo/ref data/AquaticOrganisms_Clean_final.csv", guess_max = 10000) %>% rowid_to_column()

# Master dataset for scatterplots - for Heili's tab.
aoc_v1 <- aoc %>% # start with original dataset
  #Remove 26C temperature treatment data from Jaimukar et al. 2018 (necessary to get same thresholds as in manuscript)
  filter(!(article == 42 & media.temp == 26)) %>% 
  # full dataset filters.
  mutate(effect_f = factor(case_when(effect == "Y" ~ "Yes",
                                     effect == "N" ~ "No"),
                           levels = c("No", "Yes"))) %>%
  # removing NAs to make data set nicer
  replace_na(list(size.category = 0, shape = "Not Reported", polymer = "Not Reported", life.stage = "Not Reported", chem.exp.typ.nominal = "Particle Only",
                  tech.tier.zero = "Scoring Not Applicable", risk.tier.zero = "Scoring Not Applicable"))

aoc_setup <- aoc_v1 %>% # start with original dataset
  mutate(size_f = factor(case_when(
    size.category == 1 ~ "1nm < 100nm",
    size.category == 2 ~ "100nm < 1µm",
    size.category == 3 ~ "1µm < 100µm",
    size.category == 4 ~ "100µm < 1mm",
    size.category == 5 ~ "1mm < 5mm",
    size.category == 0 ~ "Not Reported"),
    levels = c("1nm < 100nm", "100nm < 1µm", "1µm < 100µm", "100µm < 1mm", "1mm < 5mm", "Not Reported"))) %>% # creates new column with nicer names and order by size levels.
  # shape category data tidying.
  mutate(shape_f = factor(case_when(
    shape == "fiber" ~ "Fiber",
    shape == "fragment" ~ "Fragment",
    shape == "sphere" ~ "Sphere",
    shape == "Not Reported" ~ "Not Reported"),
    levels = c("Fiber", "Fragment", "Sphere", "Not Reported"))) %>% # order our different shapes.
  # polymer category data tidying.
  mutate(poly_f = factor(case_when(
    polymer == "BIO" ~ "Biopolymer",
    polymer == "EVA" ~ "Polyethylene Vinyl Acetate",
    polymer == "LTX" ~ "Latex",
    polymer == "PA" ~ "Polyamide",
    polymer == "PE" ~ "Polyethylene",
    polymer == "PC" ~ "Polycarbonate",
    polymer == "PET" ~ "Polyethylene Terephthalate",
    polymer == "PI" ~ "Polyisoprene",
    polymer == "PMMA" ~ "Polymethylmethacrylate",
    polymer == "PP" ~ "Polypropylene",
    polymer == "PS" ~ "Polystyrene",
    polymer == "PUR" ~ "Polyurethane",
    polymer == "PVC" ~ "Polyvinylchloride",
    polymer == "PLA" ~ "Polylactic Acid",
    polymer == "Not Reported" ~ "Not Reported"))) %>%
  # taxonomic category data tidying.
  mutate(org_f = factor(organism.group, levels = c("Algae", "Annelida", "Bacterium", "Cnidaria", "Crustacea",
                                                   "Echinoderm", "Fish", "Insect", "Mollusca", "Plant", "Rotifera", "Mixed"))) %>% # order our different organisms.
  mutate(lvl1_f = factor(case_when(lvl1 == "alimentary.excretory" ~ "Alimentary, Excretory",
                                   lvl1 == "behavioral.sense.neuro" ~ "Behavioral, Sensory, Neurological",
                                   lvl1 == "circulatory.respiratory" ~ "Circulatory, Respiratory",
                                   lvl1 == "community" ~ "Community",
                                   lvl1 == "fitness" ~ "Fitness",
                                   lvl1 == "immune" ~ "Immune",
                                   lvl1 == "metabolism" ~ "Metabolism",
                                   lvl1 == "microbiome" ~ "Microbiome",
                                   lvl1 == "stress" ~ "Stress"))) %>% # creates new column with nicer names.
  # Level 2 Data tidying
  mutate(lvl2_f = factor(case_when(lvl2 == "abundance"~"Abundance",
                                   lvl2 == "actinobacteria" ~ "Actinobacteria",
                                   lvl2 == "aggressivity"~"Agressivity",
                                   lvl2 == "ammonia.excretion" ~ "Ammonia Excretion",
                                   lvl2 == "bacteroidetes"~ "Bacteriodetes",
                                   lvl2 == "blood"~"Blood",
                                   lvl2 == "body.condition"~"Body Condition",
                                   lvl2 == "boldness"~"Boldness",
                                   lvl2 == "brain.histo"~"Brain Histological Abnormalities",
                                   lvl2 == "burrowing"~"Burrowing",
                                   lvl2 == "carb.metabolism"~"Carb Metabolism",
                                   lvl2 == "chemokines.cytokines"~"Chemokines",
                                   lvl2 == "circulatory"~"Circulatory",
                                   lvl2 == "detoxification"~"Detoxification",
                                   lvl2 == "development"~"Development",
                                   lvl2 == "digestion"~"Digestion",
                                   lvl2 == "digestive.enzymes"~"Digestive Enzymes",
                                   lvl2 == "digestive.tract.histo"~"Digestive Tract Histological Abnormalities",
                                   lvl2 == "diversity"~ "Diversity",
                                   lvl2 == "feeding"~ "Feeding",
                                   lvl2 == "firmicutes"~ "Firmicutes",
                                   lvl2 == "gall.bladder.histo" ~ "Gall Bladder Histological Abnormalities",
                                   lvl2 == "gen.metabolism"~ "General Metabolism",
                                   lvl2 == "gill.histo"~ "Gill Histological Abnormalities",
                                   lvl2 == "gonad.histo"~"Gonad Histological Abnormalities",
                                   lvl2 == "growth"~ "Growth",
                                   lvl2 == "immune.cells"~"Immune Cells",
                                   lvl2 == "immune.other"~"Immune Other ",
                                   lvl2 == "intestinal.permeability"~"Intestinal Permeability",
                                   lvl2 == "kidney.histo"~"Kidney Histological Abnormalities",
                                   lvl2 == "lipid.metabolism"~"Lipid Metabolism",
                                   lvl2 == "liver.histo"~"Liver Histological Abnormalities",
                                   lvl2 == "liver.kidney.products" ~ "Liver and Kidney Products",
                                   lvl2 == "locomotion"~"Locomotion",
                                   lvl2 == "mortality"~"Mortality",
                                   lvl2 == "nervous.system"~"Nervous System",
                                   lvl2 == "oxidative.stress"~"Oxidative Stress",
                                   lvl2 == "photosynthesis"~ "Photosynthesis",
                                   lvl2 == "proteobacteria"~"Protebacteria",
                                   lvl2 == "reproduction"~"Reproduction",
                                   lvl2 == "respiration"~"Respiration",
                                   lvl2 == "sexhormones"~"Sex Hormones",
                                   lvl2 == "shoaling"~"Shoaling",
                                   lvl2 == "stress"~"Stress",
                                   lvl2 == "vision.system"~"Vision System"))) %>% #Renames for widget
  #Level 3 endpoint tidying
  mutate(lvl3_f = factor(case_when(lvl3 == "11bhsd.mrnaexpression"~"11β HSD mRNA expression",
                                   lvl3 == "17bhsd.mrnaexpression"~"17β HSD mRNA expression",
                                   lvl3 == "1433e.mrnaexpression"~"1433e mRNA expression",
                                   lvl3 == "3bhsd.mrnaexpression"~"3β HSD mRNA expression",
                                   lvl3 == "5ht1.mrnaexpression"~"5ht1 mRNA expression",
                                   lvl3 == "aadd45a.mrnaexpression"~"aadd45a mRNA expression",
                                   lvl3 == "abcb1.mrnaexpression"~"abcb1 mRNA expression",
                                   lvl3 == "abcc.mrnaexpression"~"abcc mRNA expression",
                                   lvl3 == "abcc5.mrnaexpression"~"abcc5 mRNA expression",
                                   lvl3 == "abcd.mrnaexpression"~"abcd mRNA expression",
                                   lvl3 == "abnormality.occurrence"~"Abnormality Occurrence",
                                   lvl3 == "abs.rc"~"Absorption Flux/Reaction Center",
                                   lvl3 == "absorption.efficiency"~"Absorption Efficiency",
                                   lvl3 == "abundance"~"Abundance",
                                   lvl3 == "acc1.mrnaexpression"~"acc1 mRNA expression",
                                   lvl3 == "ace.diversity.index"~"ACE Diversity Index",
                                   lvl3 == "ache.activity"~"AcHE Activity",
                                   lvl3 == "ache.mrnaexpression"~"AcHE mRNA expression",
                                   lvl3 == "acid.phosphatase.activity"~"Acid Phosphatase Activity",
                                   lvl3 == "aco.mrnaexpression"~"aco mRNAexpression",
                                   lvl3 == "acrosomal.integrity"~"Acrosomal Integrity",
                                   lvl3 == "actinobacteria"~"Actinobacteria",
                                   lvl3 == "activity"~"Activity",
                                   lvl3 == "activity.during.feeding"~"Activity During Feeding",
                                   lvl3 == "aeromonas"~"Aeromonas",
                                   lvl3 == "agal.consumption"~"Agal Consumption",
                                   lvl3 == "agression.index"~"Aggression Index",
                                   lvl3 == "ahra.mrnaexpression"~"ahrα mRNA expression",
                                   lvl3 == "ahrb.mrnaexpression"~"ahrβ mRNA expression",
                                   lvl3 == "akt.mrnaexpression"~"akt mRNA expression",
                                   lvl3 == "albumin.con"~"Albumin Concentration",
                                   lvl3 == "albumin.globulin.ratio"~"Albumin to Globulin Ratio",
                                   lvl3 == "algal.density"~"Algal Density",
                                   lvl3 == "alkaline.phosphatase.activity"~"Alkaline Phosphatase Activity",
                                   lvl3 == "aproteobacteria"~"α proteobacteria",
                                   lvl3 == "alp.activity"~"Alkaline Phosphatase Activity",
                                   lvl3 == "alp.con"~"Alkaline Phosphatase Concentration",
                                   lvl3 == "alt.activity"~"Alanine Aminotransferase Activity",
                                   lvl3 == "alt.con"~"Alanine Aminotransferase Concentration",
                                   lvl3 == "amino.acid.trans"~"Amino Acid Transport",
                                   lvl3 == "ammonia.con"~"Ammonia Concentration",
                                   lvl3 == "ampka.mrnaexpression"~"ampka mRNA expression",
                                   lvl3 == "ampkb.mrnaexpression"~"ampkb mRNA expression",
                                   lvl3 == "ampkc.mrnaexpression"~"ampkc mRNA expression",
                                   lvl3 == "amylase.mrnaexpression"~"Amylanse mRNA expression",
                                   lvl3 == "antioxidant.act"~"Antioxidant Activity",
                                   lvl3 == "aox.activity"~"Alternative Oxidase Activity",
                                   lvl3 == "apo.mrnaexpression"~"Apolipoprotein mRNA expression",
                                   lvl3 == "ar.mrnaexpression"~"Androgen Receptor mRNA expression",
                                   lvl3 == "ara.mrnaexpression"~"Androgen Receptor α mRNA expression",
                                   lvl3 == "arb.mrnaexpression"~"Androgen Receptor β mRNA expression",
                                   lvl3 == "area"~"Area",
                                   lvl3 == "area.traveled"~"Area Traveled",
                                   lvl3 == "aspartate.aminotransferase.con"~"Aspartate Aminotransferase Concentration",
                                   lvl3 == "assemblage"~"Assemblage",
                                   lvl3 == "assimilation.efficiency"~"Assimilation Efficiency",
                                   lvl3 == "ast.activity"~"Aspartate Aminotransferase Activity",
                                   lvl3 == "ast.con"~"Aspartate Aminotransferase Concentration",
                                   lvl3 == "attachment"~"Attachment",
                                   lvl3 == "atubulin.mrnaexpression"~"α tubulin mRNA expression",
                                   lvl3 == "atubulin.proteinexpression"~"α tubulin Protein expression",
                                   lvl3 == "bacteroidetes"~"Bacteroidetes",
                                   lvl3 == "bcl2.mrnaexpression"~"bcl2 mRNA expression",
                                   lvl3 == "bfcod.activity"~"BFCOD Activity",
                                   lvl3 == "bile.acid.con"~"Bile Acid Concentration",
                                   lvl3 == "bilirubin.con"~"Bilirubin Concentration",
                                   lvl3 == "biomass"~"Biomass",
                                   lvl3 == "bioturbation"~"Bioturbation",
                                   lvl3 == "body.length"~"Body Length",
                                   lvl3 == "body.mass"~"Body Mass",
                                   lvl3 == "body.width"~"Body Width",
                                   lvl3 == "bproteobacteria"~"βproteobacteria",
                                   lvl3 == "brain.histo"~"Brain Histology",
                                   lvl3 == "brain.mass.index"~"Brain Mass Index",
                                   lvl3 == "brain.water.content"~"Brain Water Content",
                                   lvl3 == "burrow.length"~"Burrow Length",
                                   lvl3 == "burrow.litter"~"Burrow Litter",
                                   lvl3 == "burrow.no"~"Burrow Number",
                                   lvl3 == "burrow.vol"~"Burrow Volume",
                                   lvl3 == "burrow.walldensity"~"Burrow Wall Density",
                                   lvl3 == "burrow.wallweight"~"Burrow Wall Weight",
                                   lvl3 == "burrowing"~"Burrowing",
                                   lvl3 == "ca.mrnaexpression"~"ca mRNA expression",
                                   lvl3 == "calcification"~"Calcification",
                                   lvl3 == "calcium.con"~"Calcium Concentration",
                                   lvl3 == "capture.rate"~"Capture Rate",
                                   lvl3 == "carbohydrate.reserves"~"Carbohydrate Reserves",
                                   lvl3 == "cardiolipin.con"~"Cardiolipin Concentration",
                                   lvl3 == "cas8.mrnaexpression"~"cas8 mRNA expression",
                                   lvl3 == "casp3.mrnaexpression"~"casp3 mRNA expression",
                                   lvl3 == "casp373.mrnaexpression"~"casp373 mRNA expression",
                                   lvl3 == "casp3a.mrnaexpression"~"casp3a mRNA expression",
                                   lvl3 == "casp8.mrnaexpression"~"casp8 mRNA expression",
                                   lvl3 == "casp9.mrnaexpression"~"casp9 mRNA expression",
                                   lvl3 == "caspase3.7.mrnaexpression"~"caspase 3/7 mRNA expression",
                                   lvl3 == "cat.activity"~"Catalase Activity",
                                   lvl3 == "cat.mrnaexpression"~"Catalase mRNA expression",
                                   lvl3 == "cell.autoflourescence"~"Cell Autoflourescence",
                                   lvl3 == "cell.density"~"Cell Density",
                                   lvl3 == "cell.granularity"~"Cell Granularity",
                                   lvl3 == "cell.morphology"~"Cell Morphology",
                                   lvl3 == "cell.necrosis"~"Cell Necrosis",
                                   lvl3 == "cell.size"~"Cell Size",
                                   lvl3 == "cell.viability"~"Cell Viability",
                                   lvl3 == "cellular.energy.allocation"~"Cellular Energy Allocation",
                                   lvl3 == "chg.mrnaexpression"~"Choriogenin mRNA expression",
                                   lvl3 == "chgh.mrnaexpression"~"Choriogenin H mRNA expression",
                                   lvl3 == "chgl.mrnaexpression"~"Choriogenin L mRNA expression",
                                   lvl3 == "chlorophyll"~"Chlorophyll",
                                   lvl3 == "cholesterol.con"~"Cholesterol Concentration",
                                   lvl3 == "chrna1.mrnaexpression"~"chrna1 mRNA expression",
                                   lvl3 == "chymotrypsin.activity"~"Chymotrypsin Activity",
                                   lvl3 == "clearance.rate"~"Clearance Rate",
                                   lvl3 == "clutch.num"~"Number of Clutches",
                                   lvl3 == "clutch.size"~"Clutch Size",
                                   lvl3 == "clutch.size.first"~"First Clutch Size",
                                   lvl3 == "clutch.size.firstthree"~"First Three Clutch Sizes",
                                   lvl3 == "co2.consumption"~"CO2 Consumption",
                                   lvl3 == "condition.index"~"Condition Index",
                                   lvl3 == "congestion"~"Congestion",
                                   lvl3 == "copper.con.hepatopancreas"~"Hepatopancreas Copper Concentration",
                                   lvl3 == "cox.activity"~"Cyclooxygenase Activity",
                                   lvl3 == "coxIV.mrnaexpression"~"Cyclooxygenase IV mRNA expression",
                                   lvl3 == "cpt1.mrnaexpression"~"CPT1 mRNA expression",
                                   lvl3 == "creatine.kinase.con"~"Creatine Kinase Concentration",
                                   lvl3 == "crude.ash"~"Crude Ash",
                                   lvl3 == "crude.protein"~"Curde Protein",
                                   lvl3 == "cs.mrnaexpression"~"Chitin Synthase mRNA expression",
                                   lvl3 == "ctsl.mrnaexpression"~"Cathepsin-L mRNA expression",
                                   lvl3 == "cuznsod.mrnaexpression"~"CuZn Superoxide Dismutase mRNA expression",
                                   lvl3 == "cxcr5.mrnaexpression"~"cxcr5 mRNA expression",
                                   lvl3 == "cyp11.mrnaexpression"~"cyp11 mRNA expression",
                                   lvl3 == "cyp11a2.mrnaexpression"~"cyp11a2 mRNA expression",
                                   lvl3 == "cyp17a1.mrnaexpression"~"cyp17a1 mRNA expression",
                                   lvl3 == "cyp19a.mrnaexpression"~"cyp19a mRNA expression",
                                   lvl3 == "cyp19a2.mrnaexpression"~"cyp19a2 mRNA expression",
                                   lvl3 == "cyp19b.mrnaexpression"~"cyp19b mRNA expression",
                                   lvl3 == "cyp1a.mrnaexpression"~"cyp1a  mRNA expression",
                                   lvl3 == "cyp1a1.mrnaexpression"~"cyp1a1 mRNA expression",
                                   lvl3 == "cyp32.mrnaexpression"~"cyp32 mRNA expression",
                                   lvl3 == "cyp450.proteinexpression"~"cyp450 mRNA expression",
                                   lvl3 == "cyp4an1.mrnaexpression"~"cyp4an1 mRNA expression",
                                   lvl3 == "cyp4c33.mrnaexpression"~"cyp4c33 mRNA expression",
                                   lvl3 == "cyp4c34.mrnaexpression"~"cyp4c34 mRNA expression",
                                   lvl3 == "cytoplasmic.membrane.potential"~"Cytoplasmic Membrane Potential",
                                   lvl3 == "cytotoxicity"~"Cytotoxicity",
                                   lvl3 == "d.lactate"~"D Lactate Concentration",
                                   lvl3 == "dao.activity"~"Diamine Oxidase Activity",
                                   lvl3 == "degranulation"~"Degranulation",
                                   lvl3 == "development"~"Development",
                                   lvl3 == "developmental.abnormalities"~"Developmental Abnormalities",
                                   lvl3 == "dgat.mrnaexpression"~"Diglyceride acyltransferase mRNA expression",
                                   lvl3 == "digestive.tract.histo"~"Digestive Tract Histology",
                                   lvl3 == "dio.rc"~"Energy Dissipated/Reaction Center",
                                   lvl3 == "distance.between.fish"~"Distance Between Fish",
                                   lvl3 == "distance.during.feeding"~"Distance During Feeding",
                                   lvl3 == "diversity"~"Diversity",
                                   lvl3 == "dlac.activity"~"DLac Activity",
                                   lvl3 == "dna.damage"~"DNA Damage",
                                   lvl3 == "dopamine.con"~"Dopamine Concentration",
                                   lvl3 == "egg.size"~"Egg Size",
                                   lvl3 == "embryotoxicity"~"Embryotoxicity",
                                   lvl3 == "emergence"~"Emergence",
                                   lvl3 == "endopeptidase.activity"~"Endopeptidase Activity",
                                   lvl3 == "energy.balance"~"Energy Balance",
                                   lvl3 == "energy.reserves"~"Energy Reserves",
                                   lvl3 == "ep.mrnaexpression"~"ep mRNA expression",
                                   lvl3 == "era.mrnaexpression"~"Estrogen Receptor α mRNA expression",
                                   lvl3 == "erb.mrnaexpression"~"Estrogen Receptor β mRNA expression",
                                   lvl3 == "erk.mrnaexpression"~"erk mRNA expression",
                                   lvl3 == "erod.activity"~"EROD Activity",
                                   lvl3 == "esr2.mrnaexpression"~"esr2 mRNA expression",
                                   lvl3 == "esterase.activity"~"Esterase Activity",
                                   lvl3 == "estradiol.level"~"Estradiol Levels",
                                   lvl3 == "estradiol.testosterone.ratio"~"Estradiol to Testosterone Ratio",
                                   lvl3 == "eto.abs"~"ETO/ABS",
                                   lvl3 == "eto.rc"~"ETO/RC",
                                   lvl3 == "exopeptidase.activity"~"Exopeptidase Activity",
                                   lvl3 == "exploration"~"Exploration",
                                   lvl3 == "f0.fm"~"f0/fm",
                                   lvl3 == "f0.fv"~"f0/fv",
                                   lvl3 == "fabp6.mrnaexpression"~"fabp6 mRNA expression",
                                   lvl3 == "fadd.mrnaexpression"~"fadd mRNA expression",
                                   lvl3 == "fas.mrnaexpression"~"fas mRNA expression",
                                   lvl3 == "fatty.vac"~"Fatty Vacuoles",
                                   lvl3 == "fecundity"~"Fecundity",
                                   lvl3 == "feeding.activity"~"Feeding Activity",
                                   lvl3 == "feeding.rate"~"Feeding Rate",
                                   lvl3 == "feeding.time"~"Feeding Time",
                                   lvl3 == "fertilization"~"Fertilization",
                                   lvl3 == "fertilization.rate"~"Fertilization Rate",
                                   lvl3 == "fertilization.yield"~"Fertilization Yield",
                                   lvl3 == "fibrosis"~"Fiborsis",
                                   lvl3 == "filter.feeding"~"Filter Feeding",
                                   lvl3 == "firmicutes"~"Firmicutes",
                                   lvl3 == "food.consumption"~"Food Consumption",
                                   lvl3 == "foraging.activity"~"Foraging Activity",
                                   lvl3 == "foraging.time"~"Foraging Time",
                                   lvl3 == "foxl2.mrnaexpression"~"fox12 mRNA expression",
                                   lvl3 == "free.fatty.acid.con"~"Free Fatty Acid Concentration",
                                   lvl3 == "fshb.mrnaexpression"~"fshβ mRNA expression",
                                   lvl3 == "fshr.mrnaexpression"~"fshr mRNA expression",
                                   lvl3 == "ftzf1.mrnaexpression"~"ftzf1 mRNA expression",
                                   lvl3 == "fv'.fm'"~"fv'/fm",
                                   lvl3 == "fv.f0"~"fv/f0",
                                   lvl3 == "fv.fm"~"fv/fm",
                                   lvl3 == "fv.fm.sm.tfm"~"fv/fm/sm/tfm",
                                   lvl3 == "gadph.mrnaexpression"~"gadhph mRNA expression",
                                   lvl3 == "gall.bladder.histo"~"Gall Bladder Histology",
                                   lvl3 == "gallery.vol"~"Gallery Volume",
                                   lvl3 == "gcl.mrnaexpression"~"gcl mRNA expression",
                                   lvl3 == "gcl1.mrnaexpression"~"gcl1 mRNA expression",
                                   lvl3 == "gclcs.mrnaexpression"~"gclcs mRNA expression",
                                   lvl3 == "gclcu.mrnaexpression"~"gclcu mRNA expression",
                                   lvl3 == "genera"~"Genera",
                                   lvl3 == "gfap.mrnaexpression"~"gfap mRNA expression",
                                   lvl3 == "ggt.activity"~"Gamma-Glutamyl Transferase Activity",
                                   lvl3 == "ggt.con"~"Gamma-Glutamyl Transferase Concentration",
                                   lvl3 == "gill.histo"~"Gill Histology",
                                   lvl3 == "gk.mrnaexpression"~"gk mRNA expression",
                                   lvl3 == "globulin.con"~"Globulin Concentration",
                                   lvl3 == "glucose.con"~"Glucose Concentration",
                                   lvl3 == "glutamate.con"~"Glutamate Concentration",
                                   lvl3 == "glutathione.con"~"Glutathione Concentration",
                                   lvl3 == "glycogen.con"~"Glycogen Concentration",
                                   lvl3 == "glycogen.deplet"~"Glycogen Depletion",
                                   lvl3 == "gnrh.mrnaexpression"~"gnrh mRNA expression",
                                   lvl3 == "gnrhr.mrnaexpression"~"gnrhr mRNA expression",
                                   lvl3 == "gonad.development.index"~"Gonad Development Index",
                                   lvl3 == "gonad.histo"~"Gonad Histology",
                                   lvl3 == "gonad.regression"~"Gonad Regression",
                                   lvl3 == "gonad.somatic.index"~"Gonad Somatic Index",
                                   lvl3 == "gordonia"~"Gordonia",
                                   lvl3 == "got.activity"~"Glutamic Oxaloacetic Transaminase Activity",
                                   lvl3 == "gpt.activity"~"Glutamic Pyruvic Transaminase Activity",
                                   lvl3 == "gpx.activity"~"Glutathione Peroxidase Activity",
                                   lvl3 == "gpx.mrnaexpression"~"Glutathione Peroxidase mRNA expression",
                                   lvl3 == "gpx4b.mrnaexpression"~"gpx4b mRNA expression",
                                   lvl3 == "gr.activity"~"Glutathione Reductase Activity",
                                   lvl3 == "gr.mrnaexpression"~"Glutathione Reductase mRNA expression",
                                   lvl3 == "granulocyte.con"~"Granulocyte Concentration",
                                   lvl3 == "granulocyte.hyalinocyte.ratio"~"Granulocyte to Hyalinocyte Ratio",
                                   lvl3 == "granulocyte.oxidative.activity"~"Granulocyte Oxidative Activity",
                                   lvl3 == "granulocyte.size"~"Granulocyte Size",
                                   lvl3 == "gross.energy"~"Gross Energy",
                                   lvl3 == "growth"~"Growth",
                                   lvl3 == "growth.rate"~"Growth Rate",
                                   lvl3 == "gs.mrnaexpression"~"Glutathione Synthetase mRNA expression",
                                   lvl3 == "gsh-px.activity"~"Glutathione Peroxidase Activity",
                                   lvl3 == "gsh.activity"~"Glutathione Peroxidase Activity",
                                   lvl3 == "gsh.con"~"Glutathione Peroxidase Concentration",
                                   lvl3 == "gsp.activity"~"Glutathione Peroxidase Activity",
                                   lvl3 == "gst.activity"~"Glutathione S-Transferase Activity",
                                   lvl3 == "gst.mrnaexpression"~"Glutathione S-Transferase mRNA expression",
                                   lvl3 == "gst4.mrnaexpression"~"gst4 mRNA expression",
                                   lvl3 == "gstd.mrnaexpression"~"gstd mRNA expression",
                                   lvl3 == "gstp1.mrnaexpression"~"Glutathione S-Transferase mRNA expression",
                                   lvl3 == "gt.con"~"Glutathione Concentration",
                                   lvl3 == "gtha.mrnaexpression"~"gtha mRNA expression",
                                   lvl3 == "gusb.mrnaexpression"~"gusb mRNA expression",
                                   lvl3 == "gut.mucus.vol"~"Gut Mucus Volume",
                                   lvl3 == "gyrus.size"~"Gyrus Size",
                                   lvl3 == "h2o2.con"~"Hydorgen Peroxide Concentration",
                                   lvl3 == "h3.mrnaexpression"~"h3 mRNA expression",
                                   lvl3 == "haemocyanin.con"~"Haemocyanin Concentration",
                                   lvl3 == "hemocyte.con"~"Hemocyte Concentration",
                                   lvl3 == "haemolymph.con"~"Haemolymph Concentration",
                                   lvl3 == "hatching.success"~"Hatching Success",
                                   lvl3 == "hdl.con"~"HDL Concentration",
                                   lvl3 == "heart.rate"~"Heart Rate",
                                   lvl3 == "hepatic.somatic.index"~"Hepatic Somatic Index",
                                   lvl3 == "hepatopancreas.water.con"~"Hepatopancreas Water Concentration",
                                   lvl3 == "hex.mrnaexpression"~"hex mRNA expression",
                                   lvl3 == "hk.mrnaexpression"~"hk mRNA expression",
                                   lvl3 == "hk1.mrnaexpression"~"hk1 mRNA expression",
                                   lvl3 == "ho1.mrnaexpression"~"ho1 mRNA expression",
                                   lvl3 == "hsp70.mrnaexpression"~"hsp70 mRNA expression",
                                   lvl3 == "hsp90.mrnaexpression"~"hsp90 mRNA expression",
                                   lvl3 == "hyalinocyte.size"~"Hyalinocyte Size",
                                   lvl3 == "hylaniocyte.oxidative.activity"~"Hylaniocyte Oxidative Activity",
                                   lvl3 == "id.mrnaexpression"~"Isocitrate Dehydrogenase mRNA expression",
                                   lvl3 == "idh.activity"~"Isocitrate Dehydrogenase Activity",
                                   lvl3 == "idp.mrnaexpression"~"idp mRNA expression",
                                   lvl3 == "ifn.mrnaexpression"~"Interferon mRNA expression",
                                   lvl3 == "ifny.mrnaexpression"~"Interferon y mRNA expression",
                                   lvl3 == "igm.con"~"Immunoglobulin M Concentration",
                                   lvl3 == "igm.mrnaexpression"~"Immunoglobulin M mRNA expression",
                                   lvl3 == "il10.mrnaexpression"~"Interleukin 10 mRNA expression",
                                   lvl3 == "il17.mrnaexpression"~"Interleukin 17 mRNA expression",
                                   lvl3 == "il1a.mrnaexpression"~"Interleukin 1α mRNA expression",
                                   lvl3 == "il1b.mrnaexpression"~"Interleukin 1β mRNA expression",
                                   lvl3 == "il413a.mrnaexpression"~"Interleukin 413α mRNA expression",
                                   lvl3 == "il6.mrnaexpression"~"Interleukin 6 mRNA expression",
                                   lvl3 == "il8.mrnaexpression"~"Interluekin 8 mRNA expression",
                                   lvl3 == "ila.con"~"Interleukin 1α Concentration",
                                   lvl3 == "immobilization"~"Immobilization",
                                   lvl3 == "inflam.infiltrate"~"Inflammatory Filtrate",
                                   lvl3 == "ingestion.rate"~"Ingestion Rate",
                                   lvl3 == "inhibitory.rate"~"Inhibitory Rate",
                                   lvl3 == "instar.develop.time"~"Instar Development Time",
                                   lvl3 == "intermolt.time"~"Intermolt Time",
                                   lvl3 == "intestinal.ca.con"~"Intestinal Calcium Concentration",
                                   lvl3 == "ion.exch.ca"~"Ion Exchange (Calcium)",
                                   lvl3 == "ion.exch.k"~"Ion Exchange (Potassium)",
                                   lvl3 == "ion.exch.na"~"Ion Exchange (Sodium)",
                                   lvl3 == "ion.flow"~"Ion Flow",
                                   lvl3 == "ion.permeability"~"Ion Permeability",
                                   lvl3 == "ion.trans"~"Ion Transport",
                                   lvl3 == "jaw.tissue.damage"~"Jaw Tissue Damage",
                                   lvl3 == "jnk.mrnaexpression"~"jnk mRNA expression",
                                   lvl3 == "jump.freq"~"Jump Frequency",
                                   lvl3 == "jump.height"~"Jump Height",
                                   lvl3 == "kidney.histo"~"Kidney Histology",
                                   lvl3 == "labilization.period"~"Labilization Period",
                                   lvl3 == "lactate.con"~"Lactate Concentration",
                                   lvl3 == "ldh.activity"~"Lactate Dehydrogenase Activity",
                                   lvl3 == "ldh.con"~"Lactate Dehydrogenase Concentration",
                                   lvl3 == "ldl.con"~"LDL Concentration",
                                   lvl3 == "leaf.growth.rate"~"Leaf Growth Rate",
                                   lvl3 == "length"~"Length",
                                   lvl3 == "leukocyte.con"~"Leukocyte Concentration",
                                   lvl3 == "lhb.mrnaexpression"~"Leutinizing Hormone β mRNA expression",
                                   lvl3 == "lhr.mrnaexpression"~"Leutinizing Hormone Receptor mRNA expression",
                                   lvl3 == "lifespan"~"Lifespan",
                                   lvl3 == "light.capture"~"Light Capture",
                                   lvl3 == "light.dark.activity"~"Light/Dark Activity",
                                   lvl3 == "line.crossing"~"Line Crossing",
                                   lvl3 == "lipase.activity"~"Lipase Activity",
                                   lvl3 == "lipase.con"~"Lipase Concentration",
                                   lvl3 == "lipid.content"~"Lipid Content",
                                   lvl3 == "lipid.membrane.con"~"Lipid Membrane Concentration",
                                   lvl3 == "lipid.perox"~"Lipid Peroxidation",
                                   lvl3 == "lipid.storage.con"~"Lipid Storage Concentration",
                                   lvl3 == "lipofuscin.accum"~"Lipofuscin Accumulation",
                                   lvl3 == "liver.histo"~"Liver Histology",
                                   lvl3 == "liver.protein"~"Liver Protein",
                                   lvl3 == "lymphocyte.con"~"Lymphocyte Concentration",
                                   lvl3 == "lys.mrnaexpression"~"Lysozyme mRNA expression",
                                   lvl3 == "lysosomal.membrane.stab"~"Lysosomal Membrane Stability",
                                   lvl3 == "lysozyme.activity"~"Lysozyme Activity",
                                   lvl3 == "manf.mrnaexpression"~"manf mRNA expression",
                                   lvl3 == "mao.activity"~"Monoamine Oxidase Activity",
                                   lvl3 == "mapk.phosphorylation"~"MAPK Phosphorylation",
                                   lvl3 == "mbp.mrnaexpression"~"Myelin Basic Protein mRNA expression",
                                   lvl3 == "mbp.proteinexpression"~"Myelin Basic Protein Expression",
                                   lvl3 == "mda.con"~"Malondialdehyde Concentration",
                                   lvl3 == "meer1.mrnaexpression"~"meer 1 mRNA expression",
                                   lvl3 == "meer2.mrnaexpression"~"meer 2 mRNA expression",
                                   lvl3 == "mek.mrnaexpression"~"mek mRNA expression",
                                   lvl3 == "metamorphosis.yield"~"Metamorphosis Yield",
                                   lvl3 == "mgnrh.mrnaexpression"~"mgnrh mRNA expression",
                                   lvl3 == "min.sat.irradiance"~"Minimum Saturation Irradiance",
                                   lvl3 == "mirconuclei.freq"~"Micronuclei Frequency",
                                   lvl3 == "mitochondrial.mem.integrity"~"Mitochondrial Membrane Integrity",
                                   lvl3 == "mnsod.mrnaexpression"~"Mn Superoxide Dismutase mRNA expression",
                                   lvl3 == "moisture"~"Moisture",
                                   lvl3 == "molecular.permeability"~"Molecular Permeability",
                                   lvl3 == "molt.period.length"~"Molt Period Length",
                                   lvl3 == "mortality"~"Mortality",
                                   lvl3 == "mt10.mrnaexpression"~"mt10 mRNA expression",
                                   lvl3 == "mt20.mrnaexpression"~"mt20 mRNA expression",
                                   lvl3 == "mta.mrnaexpression"~"mta mRNA expression",
                                   lvl3 == "mtb.mrnaexpression"~"mtb mRNA expression",
                                   lvl3 == "mxr.efflux"~"MXR Efflux",
                                   lvl3 == "mytc.mrnaexpression"~"mytc mRNA expression",
                                   lvl3 == "myticina.mrnaexpression"~"myticin α mRNA expression",
                                   lvl3 == "mytilinb.mrnaexpression"~"mytilin β mRNA expression",
                                   lvl3 == "nauplius.to.adult"~"Nauplius to Adult",
                                   lvl3 == "nauplius.to.copepodid"~"Nauplius to Copepodid",
                                   lvl3 == "neutrophil.extracell.trap"~"Neutrophil Extracellular Trap",
                                   lvl3 == "ngn1.mrnaexpression"~"ngn1 mRNA expression",
                                   lvl3 == "nitrogen.conversion"~"Nitrogen Conversion",
                                   lvl3 == "no.prod"~"Nitric Oxide Production",
                                   lvl3 == "nrf2.mrnaexpression"~"nrf2 mRNA expression",
                                   lvl3 == "num.indiv.nosnails"~"Number of Individuals (No Snails)",
                                   lvl3 == "num.indiv.snails"~"Number of Individuals (With Snails)",
                                   lvl3 == "num.offspring"~"Number of Offspring",
                                   lvl3 == "number.of.molts"~"Number of Molts",
                                   lvl3 == "number.of.roots"~"Number of Roots",
                                   lvl3 == "occludin.mrnaexpression"~"occludin mRNA expression",
                                   lvl3 == "offspring.size"~"Offspring Size",
                                   lvl3 == "oocyte.diameter"~"Oocyte Diameter",
                                   lvl3 == "oocyte.num"~"Oocyte Number",
                                   lvl3 == "oxygen.consumption"~"Oxygen Consumption",
                                   lvl3 == "p38.mrnaexpression"~"p38 mRNA expression",
                                   lvl3 == "p38mapk.mrnaexpression"~"p38mapk mRNA expression",
                                   lvl3 == "p53.mrnaexpression"~"p53 mRNA expression",
                                   lvl3 == "pc.pe.ratio"~"PC/PE Ratio",
                                   lvl3 == "pcc.level"~"Protein Carbonylation Content",
                                   lvl3 == "pche.activity"~"Pseudocholinesterase Activity",
                                   lvl3 == "pepckc.mrnaexpression"~"pepckc mRNA expression",
                                   lvl3 == "percent.dstage.larvae"~"Percentage of D-stage Larvae",
                                   lvl3 == "percent.dveliger"~"Percentage of Veliger Larvae",
                                   lvl3 == "percent.motile.sperm"~"Percentage of Motile Sperm",
                                   lvl3 == "percent.tank.used"~"Percentage of Tank Used",
                                   lvl3 == "percent.trocophore.larvae"~"Percentage of Trocophore Larvae",
                                   lvl3 == "pericardial.area"~"Pericardial Area",
                                   lvl3 == "perox.activity"~"Peroxidase Activity",
                                   lvl3 == "perox.mrnaexpression"~"Peroxidase mRNA expression",
                                   lvl3 == "petb.mrnaexpression"~"petb mRNA expression",
                                   lvl3 == "pglycoprotein.activity"~"Pglycoprotein Activity",
                                   lvl3 == "pgp.mrnaexpression"~"pgp mRNA expression",
                                   lvl3 == "pgrp.mrnaexpression"~"pgrp mRNA expression",
                                   lvl3 == "phagocytosis"~"Phagocytosis",
                                   lvl3 == "phenoloxydase.activity"~"Phenoloxydase Activity",
                                   lvl3 == "phosphorus.con"~"Phosphorus Concentration",
                                   lvl3 == "phosphatase.activity"~"Phosphatase Activity",
                                   lvl3 == "phosphatidylcholine.con"~"Phosphatidylcholine Concentration",
                                   lvl3 == "phosphatidylethanolamine.con"~"Phosphatidylethanolamine Concentration",
                                   lvl3 == "phosphatidylinositol.ceramide.amino.ethylphosphonate.con"~"PCAE Concentration",
                                   lvl3 == "phosphatidylserine.con"~"Phosphatidylserine Concentration",
                                   lvl3 == "photochemical.efficiency"~"Photochemical Efficiency",
                                   lvl3 == "pk.mrnaexpression"~"pk mRNA expression",
                                   lvl3 == "polyp.activity"~"Polyp Activity",
                                   lvl3 == "population.size"~"Population Size",
                                   lvl3 == "potassium.con"~"Potassium Concentration",
                                   lvl3 == "ppara.mrnaexpression"~"pparα mRNA expression",
                                   lvl3 == "ppary.mrnaexpression"~"ppary mRNA expression",
                                   lvl3 == "prdx1.mrnaexpression"~"prdx1 mRNA expression",
                                   lvl3 == "prdx2.mrnaexpression"~"prdx2 mRNA expression",
                                   lvl3 == "prdx3.mrnaexpression"~"prdx3 mRNA expression",
                                   lvl3 == "prdx5.mrnaexpression"~"prdx5 mRNA expression",
                                   lvl3 == "predation.risk"~"Predation Risk",
                                   lvl3 == "predatory.performance"~"Predatory Performance",
                                   lvl3 == "pregnancy.rate"~"Pregnancy Rate",
                                   lvl3 == "protein.carbonyl.con"~"Protein Carbonyl Concentration",
                                   lvl3 == "protein.con"~"Protein Concentration",
                                   lvl3 == "protein.reserves"~"Protein Reserves",
                                   lvl3 == "proteobacteria"~"proteobacteria",
                                   lvl3 == "ps2"~"Photosystem II",
                                   lvl3 == "psaa.mrnaexpression"~"psaa mRNA expression",
                                   lvl3 == "psab.mrnaexpression"~"psab mRNA expression",
                                   lvl3 == "psba.mrnaexpression"~"psba mRNA expression",
                                   lvl3 == "psbl.mrnaexpression"~"psbl mRNA expression",
                                   lvl3 == "pseudomonas"~"pseudomonas",
                                   lvl3 == "rapid.light.resp.curve"~"Rapid Light Response Curve",
                                   lvl3 == "rbc.con"~"Red Blood Cell Concentration",
                                   lvl3 == "rbcl.mrnaexpression"~"rbcl mRNA expression",
                                   lvl3 == "rc.abs"~"RC/ABS",
                                   lvl3 == "redox.state"~"Redox State",
                                   lvl3 == "relative.e.transport.rate"~"Relative e Transport State",
                                   lvl3 == "relocate"~"Relocation",
                                   lvl3 == "reproduction.factor"~"Reproduction Factor",
                                   lvl3 == "reproduction.inhibition"~"Reproduction Inhibition",
                                   lvl3 == "reproductive.freq"~"Reproductive Frequency",
                                   lvl3 == "resistance.time"~"Resistance Time",
                                   lvl3 == "respiration.rate"~"Respiration Rate",
                                   lvl3 == "respiratory.burst"~"Respiratory Burst",
                                   lvl3 == "retrmax"~"Maximal Election Transport Rate",
                                   lvl3 == "root.length"~"Root Length",
                                   lvl3 == "ros.prod"~"Reactive Oxygen Species Production",
                                   lvl3 == "s100a1.mrnaexpression"~"s100a1 mRNA expression",
                                   lvl3 == "saa.mrnaexpression"~"saa mRNA expression",
                                   lvl3 == "sd.mrnaexpression"~"sd mRNA expression",
                                   lvl3 == "segment.regeneration"~"Segment Regeneration",
                                   lvl3 == "serotonin.con"~"Serotonin Concentration",
                                   lvl3 == "settling.rate"~"Settling Rate",
                                   lvl3 == "sex.ratio"~"Sex Ratio",
                                   lvl3 == "sexual.maturity"~"Sexual Maturity",
                                   lvl3 == "sgst.mrnaexpression"~"sgst mRNA expression",
                                   lvl3 == "shannon.diversity.index"~"Shannon Diversity Index",
                                   lvl3 == "shell.growth"~"Shell Growth",
                                   lvl3 == "shell.length"~"Shell Length",
                                   lvl3 == "shell.weight"~"Shell Weight",
                                   lvl3 == "size"~"Size",
                                   lvl3 == "sod.activity"~"Superoxide Dismutase Activity",
                                   lvl3 == "sod.mrnaexpression"~"sod mRNA expression",
                                   lvl3 == "sod1.mrnaexpression"~"sod1 mRNA expression",
                                   lvl3 == "sod3.mrnaexpression"~"sod3 mRNA expression",
                                   lvl3 == "sodium.con"~"Sodium Concentration",
                                   lvl3 == "species.diversity"~"Species Diversity",
                                   lvl3 == "species.richness"~"Species Richness",
                                   lvl3 == "specific.growth.rate"~"Specific Growth Rate",
                                   lvl3 == "sperm.velocity"~"Sperm Velocity",
                                   lvl3 == "srebp1a.mrnaexpression"~"srebp1a mRNA expression",
                                   lvl3 == "star.mrnaexpression"~"star mRNA expression",
                                   lvl3 == "sterol.con"~"Sterol Concentration",
                                   lvl3 == "surface.area"~"Surface Area",
                                   lvl3 == "swim.distance"~"Swim Distance",
                                   lvl3 == "swim.speed"~"Swim Speed",
                                   lvl3 == "swimbladder.inflation"~"Swimbladder Inflation",
                                   lvl3 == "synapsin.mrnaexpression"~"synapsin mRNA expression",
                                   lvl3 == "taxa.num"~"Number of Taxa",
                                   lvl3 == "testosterone.level"~"Testosterone Levels",
                                   lvl3 == "tgfb.mrnaexpression"~"tgfb mRNA expression",
                                   lvl3 == "thiol.con"~"Thiol Concentration",
                                   lvl3 == "thrombocyte.con"~"Thrombocyte Concentration",
                                   lvl3 == "time.to.hatch"~"Time to Hatch",
                                   lvl3 == "time.to.maturation"~"Time to Maturation",
                                   lvl3 == "time.to.nauplius"~"Time to Nauplius",
                                   lvl3 == "time.to.offspring"~"Time to Offspring",
                                   lvl3 == "tlr.mrnaexpression"~"tlr mRNA expression",
                                   lvl3 == "tnfa.mrnaexpression"~"tnfα mRNA expression",
                                   lvl3 == "tnfsf13b.mrnaexpression"~"tnfsf13b mRNA expression",
                                   lvl3 == "total.distance.traveled"~"Total Distance Travelled",
                                   lvl3 == "tp53.con"~"tp53 Concentration",
                                   lvl3 == "tp53.mrnaexpression"~"tp53 mRNA expression",
                                   lvl3 == "tph2.mrnaexpression"~"tph2 mRNA expression",
                                   lvl3 == "tricellulin.mrnaexpression"~"tricellulin mRNA expression",
                                   lvl3 == "triglyceride.cholesterol.ratio"~"Tryglyceride Cholesterol Ratio",
                                   lvl3 == "triglycerides.con"~"Triglyceride Concentration",
                                   lvl3 == "tro.rc"~"TRO/RC",
                                   lvl3 == "trypsin.activity"~"Trypsin Activity",
                                   lvl3 == "tumor.promotion"~"Tumor Promotion",
                                   lvl3 == "turn.angle.light.dark.activity"~"Turn Angle, Light/Dark Activity",
                                   lvl3 == "ucp1.mrnaexpression"~"ucp1 mRNA expression",
                                   lvl3 == "ugst.mrnaexpression"~"ugst mRNA expression",
                                   lvl3 == "uricacid.con"~"Uric Acid Concentration",
                                   lvl3 == "viserosomatic.index"~"Viserosomatic Index",
                                   lvl3 == "volume"~"Volume",
                                   lvl3 == "vtg.mrnaexpression"~"vtg mRNA expression",
                                   lvl3 == "vtg.proteinexpression"~"vtg protein expression",
                                   lvl3 == "vtg1.mrnaexpression"~"vtg1 mRNA expression",
                                   lvl3 == "vtg2.mrnaexpression"~"vtg2 mRNA expression",
                                   lvl3 == "weight.change"~"Weight Change",
                                   lvl3 == "wgst.mrnaexpression"~"wgst mRNA expression",
                                   lvl3 == "wnt16.mrnaexpression"~"wnt16 mRNA expression",
                                   lvl3 == "wnt4.mrnaexpression"~"wnt4 mRNA expression",
                                   lvl3 == "yproteobacteria"~"y proteobacteria",
                                   lvl3 == "zfblue.mrnaexpression"~"zfblue mRNA expression",
                                   lvl3 == "zfrho.mrnaexpression"~"zfrho mRNA expression",
                                   lvl3 == "zo1.mrnaexpression"~"zo1 mRNA expression",
                                   lvl3 == "zoozanthelle.density"~"Zoozanthelle Density"))) %>%
  mutate(bio_f = factor(case_when( #Bio Org Data Tidying
    bio.org == "subcell"~"Subcell",
    bio.org == "cell"~"Cell",
    bio.org == "tissue" ~ "Tissue",
    bio.org == "organism"~"Organism",
    bio.org == "population"~ "Population")))%>%
  mutate(vivo_f = factor(case_when(invitro.invivo == "invivo"~"In Vivo",
                                   invitro.invivo == "invitro"~"In Vitro")))%>% ##Renames for widget (Not using a widget right now, but saving for human health database)
  mutate(life_f = factor(case_when(life.stage == "Early"~"Early",
                                   life.stage == "Juvenile"~"Juvenile",
                                   life.stage == "Adult"~"Adult",
                                   life.stage == "Not Reported"~"Not Reported")))%>% #Renames for widget
  mutate(env_f = factor(case_when(environment == "Freshwater"~"Freshwater",
                                  environment == "Marine" ~ "Marine",
                                  environment == "Terrestrial" ~ "Terrestrial"))) %>%
  mutate(species = if_else(genus == "Mytilus" & is.na(species), "species", species)) %>% 
  mutate(species_f = as.factor(paste(genus,species))) %>% 
  mutate(dose.mg.L.master.converted.reported = factor(dose.mg.L.master.converted.reported)) %>%
  mutate(dose.particles.mL.master.converted.reported = factor(dose.particles.mL.master.converted.reported)) %>% 
  #a few of the mortality LC50's are misclassified as EC50. fixing here
  mutate(effect.metric = case_when(
    effect.metric == "NOEC" ~ "NOEC",
    effect.metric == "LOEC" ~ "LOEC",
    effect.metric == "EC50" & effect.score == "6" ~ "LC50", #reclassify
    effect.metric == "EC50" ~ "EC50",
    effect.metric == "EC10" ~ "EC10",
    effect.metric == "IC10" ~ "IC10",
    effect.metric == "HONEC" ~ "HONEC",
    effect.metric == "IC50" ~ "IC50",
    effect.metric == "LC50" ~ "LC50"
  )) %>% 
  mutate(effect.metric = factor(effect.metric)) %>% #factorize
  mutate(acute.chronic_f = factor(case_when(af.time == 10 ~ "Acute",
                                            af.time == 1 ~ "Chronic"))) %>% #factorize assessment factor time into chronic/acute
  mutate(tier_zero_tech_f = factor(case_when(tech.tier.zero == "Fail" ~ "Red Criteria Failed",
                                             tech.tier.zero == "Pass" ~ "Red Criteria Passed",
                                             tech.tier.zero == "Scoring Not Applicable" ~ "Scoring Not Applicable"))) %>% 
  mutate(tier_zero_risk_f = factor(case_when(risk.tier.zero == "Fail" ~ "Red Criteria Failed",
                                             risk.tier.zero == "Pass" ~ "Red Criteria Passed",
                                             risk.tier.zero == "Scoring Not Applicable" ~ "Scoring Not Applicable"))) %>%
  mutate(weather.biofoul_f = factor(case_when(weather.biofoul == "Y" ~ "Yes",
                                              weather.biofoul == "N" ~ "No"))) %>%
  #Paste the range of treatments used within each study - there is often a range because this is different depending on the endpoint
  group_by(doi) %>% 
  mutate(treatment_range = ifelse(min(treatments) == max(treatments), paste0(treatments), paste0(min(treatments),"-",max(treatments)))) %>% 
  ungroup() %>% 
  #annotate whether max size ingest was estimated or reported
  mutate(max.size.ingest.reported.estimated = case_when(
    is.na(max.size.ingest.mm) ~ "estimated",  # If the value is NA, label it "estimated"
    TRUE ~ "reported"                   # Otherwise, label it "converted"
  )) %>% 
  #calculate maximum ingestible size (if not already in database)
  mutate(max.size.ingest.mm = ifelse(is.na(max.size.ingest.mm), 
                                     10^(beta_log10_body_length * log10(body.length.cm) - body_length_intercept) * 10,  #(Jamm et al 2020 Nature paper)correction for cm to mm
                                     max.size.ingest.mm)) %>%  # if already present, just use that
  mutate(max.size.ingest.um = 1000 * max.size.ingest.mm) %>%  #makes it less confusing below
  #Make factor for experiment type
  mutate(exp_type_f = factor(case_when(leachate.only == "Y" ~ "Leachate",
                                       chem.exp.typ.nominal == "Particle Only" ~ "Particle Only",
                                       chem.exp.typ.nominal == "co.exp" ~ "Chemical Co-Exposure",
                                       chem.exp.typ.nominal == "sorbed" ~ "Chemical Transfer"))) %>%
  
  #### Recalculation of surface area and volume based on shape ####

## Assign average length to width ratio (R.Ave) based on compartment (values from Kooi et al. 2021 SI)
  mutate(R.ave = case_when(environment == "Marine" & exposure.route == "water" ~ R.ave.water.marine,
                           environment == "Marine" & exposure.route == "sediment" ~ R.ave.sediment.marine,
                           environment == "Freshwater" & exposure.route == "water" ~ R.ave.water.freshwater,
                           environment == "Freshwater" & exposure.route == "sediment" ~ R.ave.sediment.freshwater,
                           T ~ NA)) %>% # if doesn't meet these conditions, annotate as NA
  #calculate surface area based on shape
#calculate surface area based on shape
mutate(particle.surface.area.um2 = case_when(shape == "sphere" ~ particle.surface.area.um2,
                                             shape == "fiber" & is.na(size.width.um.used.for.conversions) ~ SAfnx_fiber(width = 15, length = size.length.um.used.for.conversions), #assum 15 um width (kooi et al 2021)
                                             shape == "fiber" & !is.na(size.width.um.used.for.conversions) ~ SAfnx_fiber(width = size.width.um.used.for.conversions, length = size.length.um.used.for.conversions), #if width is known
                                             shape == "fragment" ~ SAfnx(a = size.length.um.used.for.conversions,
                                                                         b = R.ave * size.length.um.used.for.conversions,
                                                                         c = R.ave * 0.67 * size.length.um.used.for.conversions))) %>% 
  mutate(particle.volume.um3 = case_when(shape == "sphere" ~ particle.volume.um3, #sphere volume is correct in excel
                                         shape == "fiber" & is.na(size.width.um.used.for.conversions) ~ volumefnx_fiber(width = 15, length = size.length.um.used.for.conversions), #assume 15 um as width (kooi et al 2021)
                                         shape == "fiber" & !is.na(size.width.um.used.for.conversions) ~ volumefnx_fiber(width = size.width.um.used.for.conversions, length = size.length.um.used.for.conversions), #if width reported
                                         shape == "fragment" ~ volumefnx(R = R.ave, L = size.length.um.used.for.conversions))) %>% 
  
  #calcualte dose metrics accordingly
  mutate(dose.surface.area.um2.mL.master = particle.surface.area.um2 * dose.particles.mL.master) %>% 
  mutate(particle.surface.area.um2.mg = particle.surface.area.um2 / mass.per.particle.mg) %>% 
  
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
    shape == "sphere" ~ size.length.min.um.used.for.conversions, #all dims same
    shape == "fiber" ~ R.ave * size.length.min.um.used.for.conversions, #median holds for all particles (Kooi et al 2021)
    shape == "Not Reported" ~ R.ave * size.length.min.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
    shape == "fragment" ~ R.ave * size.length.min.um.used.for.conversions)) %>% # average width to length ratio in the marine environment (kooi et al 2021)
  mutate(size.height.min.um.used.for.conversions = case_when(
    shape == "sphere" ~ size.length.min.um.used.for.conversions, #all dims same
    shape == "Not Reported" ~ R.ave * 0.67 * size.length.min.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
    shape == "fiber" ~  R.ave * size.length.min.um.used.for.conversions, #height same as width for fibers
    shape == "fragment" ~ R.ave * 0.67 * size.length.min.um.used.for.conversions)) %>% # average width to length ratio in the marine environment AND average height to width ratio (kooi et al 2021)
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
    shape == "Not Reported" ~ R.ave * 0.67 * size.length.max.um.used.for.conversions, # average width to length ratio in the marine environment (kooi et al 2021)
    shape == "fiber" ~ R.ave * size.length.max.um.used.for.conversions, #hieght same as width
    shape == "fragment" ~ R.ave * 0.67 * size.length.max.um.used.for.conversions)) %>%  # average width to length ratio in the marine environment AND average height to width ratio (kooi et al 2021)
  
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
                                                 p = density.g.cm3)) %>% #equation usess g/cm3
  mutate(mass.per.particle.mg.max = massfnx_poly(length = size.length.max.um.used.for.conversions,
                                                 width = size.width.max.um.used.for.conversions,
                                                 p = density.g.cm3)) %>%   #equation usess g/cm3
  
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
  
  #Additional tidying for nicer values
  mutate(authors = gsub(".", " & ", as.character(authors), fixed = TRUE)) %>% 
  mutate(exposure.media = gsub(".", " ", as.character(exposure.media), fixed = TRUE)) %>%
  mutate(detergent = gsub(".", " ", as.character(detergent), fixed = TRUE)) %>%
  mutate(chem.add.nominal = gsub(".", " ", as.character(chem.add.nominal), fixed = TRUE)) %>%
  mutate(exposure.route = gsub(".", " ", as.character(exposure.route), fixed = TRUE)) %>% 
  mutate(sol.rinse = gsub(".", " ", as.character(sol.rinse), fixed = TRUE)) %>%
  mutate(sol.rinse = if_else(sol.rinse == "N", "No", sol.rinse)) %>% 
  mutate(uptake.valid.method = gsub(".", " ", as.character(uptake.valid.method), fixed = TRUE)) %>% 
  mutate(clean.method = gsub(".", " ", as.character(clean.method), fixed = TRUE)) %>% 
  mutate(clean.method = if_else(clean.method == "N", "Not Cleaned", clean.method)) %>% 
  mutate(particle.behavior = gsub(".", " ", as.character(particle.behavior), fixed = TRUE)) %>% 
  mutate(particle.behavior = if_else(particle.behavior == "N", "Not Evaluated", particle.behavior)) %>%
  mutate(tissue.distribution = gsub(".", " ", as.character(tissue.distribution), fixed = TRUE))

#### Overview Setup Continued ####

#Endpoint Categorization setup
aoc_endpoint <- aoc_setup %>% 
  group_by(lvl1_f,lvl2_f,lvl3_f,bio_f) %>% 
  summarise()


#### Search Setup ####

aoc_search <- aoc_setup %>%
  #general
  dplyr::select(doi, authors, year, tier_zero_tech_f, tier_zero_risk_f, species_f, org_f, env_f, life_f, vivo_f, sex, body.length.cm, max.size.ingest.mm,
                #experimental parameters
                exp_type_f, exposure.route, mix, negative.control, reference.material, exposure.media, solvent, detergent,
                media.ph, media.sal.ppt, media.temp, media.temp.min, media.temp.max, exposure.duration.d, acute.chronic_f,
                treatments, replicates, sample.size, dosing.frequency, chem.add.nominal, chem.add.dose.mg.L.nominal, chem.add.dose.mg.L.measured,
                #master doses
                dose.particles.mL.master, dose.particles.mL.master.converted.reported, dose.mg.L.master, dose.mg.L.master.converted.reported,
                dose.um3.mL.master, dose.um2.mL.master, dose.um2.ug.mL.master,
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
                #reported doses
                dose.particles.mL.nominal, dose.particles.mL.min.nominal,
                dose.particles.mL.max.nominal, dose.mg.L.nominal, dose.mg.kg.sed.nominal, dose.mg.kg.food.nominal,
                dose.particles.kg.food.nominal, dose.percent.sed.nominal, dose.particles.m2.nominal, dose.mg.m2.nominal, dose.percent.food.nominal,
                dose.particles.kg.sed.nominal, dose.particles.mL.measured, dose.mg.L.measured, dose.mg.kg.sed.measured, dose.mg.kg.food.measured,
                dose.particles.kg.food.measured, dose.particles.kg.food.min.measured, dose.particles.kg.food.max.measured,
                dose.percent.sed.measured, dose.particles.m2.measured) %>%
  #rename 'master' dose columns
  rename("particles/mL (master)" = dose.particles.mL.master, "particles/mL (master), reported or converted" = dose.particles.mL.master.converted.reported,
         "μg/mL (master)" = dose.mg.L.master, "μ/mL (master), reported or converted" = dose.mg.L.master.converted.reported,
         "μm^3/mL (master)" = dose.um3.mL.master, "μm^2/mL (master)" = dose.um2.mL.master, "μm/ug/mL (master)" = dose.um2.ug.mL.master) %>% 
  #rename other dose columns
  rename("particles/mL - nominal" = dose.particles.mL.nominal, "particles/mL (min) - nominal" = dose.particles.mL.min.nominal,
         "particles/mL (max) - nominal" = dose.particles.mL.max.nominal, "mg/L - nominal" = dose.mg.L.nominal, "mg/kg sediment - nominal" = dose.mg.kg.sed.nominal, 
         "mg/kg food - nominal" = dose.mg.kg.food.nominal,
         "particles/kg food - nominal" = dose.particles.kg.food.nominal, "% sediment - nominal" = dose.percent.sed.nominal, 
         "particles/m^2 - nominal" = dose.particles.m2.nominal, "mg/m^2 - nominal" = dose.mg.m2.nominal, "% food - nominal" = dose.percent.food.nominal,
         "particles/kg sediment - nominal" = dose.particles.kg.sed.nominal, 
         
         "particles/mL - measured" = dose.particles.mL.measured, "mg/L - measured" = dose.mg.L.measured, "mg/kg sediment - measured" = dose.mg.kg.sed.measured, 
         "mg/kg food - measured" = dose.mg.kg.food.measured,
         "particles/kg food - measured" = dose.particles.kg.food.measured, "particles/kg food (min) - measured" = dose.particles.kg.food.min.measured, 
         "particles/kg food (max) - measured" = dose.particles.kg.food.max.measured,
         "% sediment - measured" = dose.percent.sed.measured, "particles/m^2 - measured" = dose.particles.m2.measured)

#Turn all character strings into factors if they aren't already so they are searchable via dropdown
aoc_search[sapply(aoc_search, is.character)] <- lapply(aoc_search[sapply(aoc_search, is.character)], as.factor)

#Rename table columns

aoc_search <- aoc_search %>% 
  dplyr::rename('DOI' = doi,'Authors' = authors, 'Year' = year, 'Technical "Red Criteria"' = tier_zero_tech_f, 
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
                'Added Chemical Dose (nominal)' = chem.add.dose.mg.L.nominal,
                'Added Chemical Dose (measured)' = chem.add.dose.mg.L.measured, 
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
                'Microplastic Diversity Score' = risk.19, 'Exposure Time Score' = risk.20)

#### SSD Setup ####

# Master dataset for SSDs
aoc_z <- aoc_setup %>% 
  # environment category data tidying.
  mutate(environment.noNA = replace_na(environment, "Not Reported")) %>% # replaces NA to better relabel.
  mutate(env_f = factor(environment.noNA, levels = c("Marine", "Freshwater", "Terrestrial", "Not Reported"))) 
# final cleanup and factoring  
aoc_z$Species <- as.factor(paste(aoc_setup$genus,aoc_setup$species)) #must make value 'Species" (uppercase)
aoc_z$Group <- as.factor(aoc_z$organism.group) #must make value "Group"
aoc_z$Group <- fct_explicit_na(aoc_z$Group) #makes sure that species get counted even if they're missing a group


#### Screening Setup ####

aoc_quality <- aoc_setup %>%
  filter(tier_zero_tech_f != "Scoring Not Applicable") %>% 
  filter(tier_zero_risk_f != "Scoring Not Applicable") %>% 
  mutate(Study = paste0(authors, " (", year,")")) %>%
  mutate(Study_plus = as.factor(paste0(authors, " (", year,")", " (",doi,")"))) %>%
  distinct(Study, Study_plus, doi, treatment_range, tech.a1, tech.a2, tech.a3, tech.a4, tech.a5, tech.a6, tech.1, tech.2, tech.3, tech.4, tech.5,
           tech.6, tech.7, tech.8, tech.9, tech.10, tech.11, tech.12, risk.13, risk.14, risk.15, risk.16, risk.17, risk.18, risk.19, risk.20,
           lvl1_f, lvl2_f, bio_f, effect_f, life_f, poly_f, shape_f, size_f, species_f, env_f, org_f, acute.chronic_f, tier_zero_tech_f, tier_zero_risk_f) %>%   
  
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


###### SAVE DATA #####
# saveRDS(aoc, file = "aoc.RDS")
# saveRDS(aoc_endpoint, file = "aoc_endpoint.RDS")
# saveRDS(aoc_quality, file = "aoc_quality.RDS")
# saveRDS(aoc_search, file = "aoc_search.RDS")
# saveRDS(aoc_setup, file = "aoc_setup.RDS")
# saveRDS(aoc_v1, file = "aoc_v1.RDS")
# saveRDS(aoc_z, file = "aoc_z.RDS")

return(aoc_setup)
} #end function
