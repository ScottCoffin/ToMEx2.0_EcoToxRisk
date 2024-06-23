# aim: prepare data from Reisser et al. 2014 doi:10.5061/dryad.mt0j5
# Primary article: https://doi.org/10.1371/journal.pone.0080466


# Description of the variables of plastics_info.csv matrix (see its heading/first line)
# 
# Plastic --> plastic ID number
# NetStation --> net station identification number (where 3 consecutive net tows were undertaken); goes from 1 to 57
# Replicate --> indicate if it is from the first (1), second (2), or third (3) net tow undertaken at the net station
# Type --> type of plastic (e.g. hard, soft, styrofoam)
# Length --> plastic size (in mm)
# Colour --> plastic colour
# PolymerType --> plastic polymer as determined by FT-IR spectra

library(dplyr)
library(tidyr)

dat= read.csv("Environmental_Data_Collection/00_Extraction_finished/ReisserEtAl2014/plastics_info.csv", stringsAsFactors = TRUE)
str(dat)
head(dat)

dat = dat %>% mutate(across(NetStation, as.factor))

# include only particles <= 5 mm
dat = dat[dat$Length <= 5,]
dat = droplevels(dat)
str(dat)
summary(dat) #no NA values :)

# polymer types
levels(dat$PolymerType)

# number of samples per station
summary(dat$NetStation)

# According to the primary article, there are seven separate sampling events. 
# Aggregate data according to these events to get 1 value per sampling event for the final NMDS plot

dat$SamplingTrip = NA
dat$SamplingTrip[dat$NetStation %in% 1:4] = 1
dat$SamplingTrip[dat$NetStation %in% 5:14] = 2
dat$SamplingTrip[dat$NetStation %in% 15:23] = 3
dat$SamplingTrip[dat$NetStation %in% 24:27] = 4
dat$SamplingTrip[dat$NetStation %in% 28:39] = 5
dat$SamplingTrip[dat$NetStation %in% 40:47] = 6
dat$SamplingTrip[dat$NetStation %in% 48:57] = 7
dat = dat %>% mutate(across(SamplingTrip, as.factor))

# calculate summary stats

lengths = dat %>%  group_by(SamplingTrip) %>% summarise(mean.length = mean(Length), median.length = median(Length)) 

polymers = dat %>%  group_by(SamplingTrip, PolymerType) %>% summarise(n = n()) %>%  mutate(prop = n/sum(n))  %>% pivot_wider(names_from = PolymerType, values_from = prop) %>% 
  summarise(across(everything(), function(x) sum(x, na.rm = TRUE)))

(final = merge(lengths, polymers, by = "SamplingTrip"))                  
final$other = final$`NA` + final$`Ethylene Vinyl Acetate`
final
