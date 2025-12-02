# load packages
library(plyr)


# import data
dat = readRDS("data/output/aoc_final.RDS")
dat = droplevels(dat)
str(dat)

# create particle and treatment ID
dat$particleID = paste0(dat$doi, dat$poly_f, dat$shape_f, dat$charge, dat$size.length.um.used.for.conversions, 
                        dat$functional.group, dat$weather.biofoul_f)
dat$treatmentID = paste0(dat$particleID, dat$Species, dat$life_f, dat$vivo_f, dat$sex, dat$exposure.route)

# To avoid over-representation of particles for which many endpoints have been measured within the same experimental treatment
# prepare dataset with only one entry per unique treatment-----
dat = dat[!duplicated(dat$treatmentID), ] 

# Check factor level names
levels(dat$poly_f)
levels(dat$shape_f)
levels(dat$sodium.azide)
dat$sodium.azide <- plyr::revalue(dat$sodium.azide, c("unknown" = "Unknown"))
levels(dat$DOM_present)
levels(as.factor(dat$charge))
dat$charge = plyr::revalue(dat$charge, c("negative" = "Negative", "positive" = "Positive"))
levels(as.factor(dat$functional.group))

# Create polymer groups: all polymers that were tested only once are group "Others"
others.list = as.data.frame(table(dat$poly_f))$Var1[as.data.frame(table(dat$poly_f))$Freq == 1]
dat$poly_group = as.character(dat$poly_f)
dat$poly_group[dat$poly_f %in% others.list] = "Other"
dat$poly_group[dat$poly_group == "Mix - See Original Study"] = "Polymer mix"
dat$poly_group = as.factor(dat$poly_group)

# Replace NAs in charge and functional.group with "Unknown" and "None", respectively
dat$charge[is.na(dat$charge)] = "Unknown"
dat$functional.group[is.na(dat$functional.group)] = "None"

# save data
saveRDS(dat, "data/output/characteristics and NMDS/prepared_data.RDS")

### END ###