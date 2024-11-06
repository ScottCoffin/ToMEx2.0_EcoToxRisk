# import data
dat = readRDS("../data/output/aoc_final.RDS")
dat = droplevels(dat)
str(dat)

# create particle and treatment ID
dat$particleID = paste0(dat$doi, dat$poly_f, dat$shape_f, dat$charge, dat$size.length.um.used.for.conversions, 
                        dat$functional.group, dat$weather.biofoul_f)
dat$treatmentID = paste0(dat$particleID, dat$Species, dat$life_f, dat$vivo_f, dat$sex, dat$exposure.route)

# To avoid over-representation of particles for which many endpoints have been measured within the same experimental treatment
# prepare dataset with only one entry per unique treatment-----
dat = dat[!duplicated(dat$treatmentID), ] 

# save data
saveRDS(dat, "Data/prepared_data.RDS")

### END ###