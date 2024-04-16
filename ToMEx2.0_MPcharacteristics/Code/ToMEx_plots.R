
library(dplyr)
source("Code/Helper_boxplot_Kooi2021.R")


dat = readRDS("Data/RDS/aoc_z_tomex2.RDS")


# import data from Kooi et al. 2021:
dat.lengths = read.csv("Data/KooiEtAl2021_particle_length.csv")

str(dat)
summary(dat)

# !! split between freshwater and marine! # compare to aquatic env data ??? # $ env_f 
# The reason is that the composition of microplastics in the environment will differ between these two compartments

# each particle type only once in the dataset > create particle ID
# The problem is that in some experiments, several endpoints were measured in the same organisms in the same setup using the same particles
# if we do not account for these, particle properties from particles with many endpoint measurements would be overrepresented
# We need a unique experiment ID with 1 or more unique particle IDs
# This allows for the same particle to be used in different experiments within the same publication. These are counted repeatedly
# > summary on experiment-particle level!

# create unique particle-experiment-ID
# experiment is defined by: doi, Species, life_f, vivo_f, sex, exposure.route
# particle is defined by: doi, poly_f, shape_f, charge, size.length.um.used.for.conversions, functional.group, weather.biofoul_f

dat$particleID = paste0(dat$doi, dat$poly_f, dat$shape_f, dat$charge, dat$size.length.um.used.for.conversions, 
                        dat$functional.group, dat$weather.biofoul_f)
dat$treatmentID = paste0(dat$particleID, dat$Species, dat$life_f, dat$vivo_f, dat$sex, dat$exposure.route)

# number of unique treatments:
length(unique(dat$treatmentID))

# prepare dataset with only 1 entry per unique treatment-----

treatdat = dat[!duplicated(dat$treatmentID), ] # takes only the first entry per treatmentID
dim(treatdat)


# We need summary stats and plot:
# - Polymer = poly_f
# - density.g.cm3
# - size.length.um.used.for.conversions
# - size_f 
# - shape_f
# - charge
# - zetapotential.mV

# - dose.particles.mL.master
# - dose.mg.L.master 


## Aim: Is the composition of MP properties in ToMEx different to those in natural environmental samples?
# this could be looked at in terms of: number of datapoints with specific properties in ToMEx compared to number of particles in env. samples 
# >> need per particle data
# or in terms of: concentrations of particles with specific properties in ToMEx compared to concentrations of same props in env. samples 
# >> need per sample frequency distribution of props


# Polymer-----

## number of studies per polymer type
polymer = treatdat %>% group_by(poly_f) %>% summarize(count = n())
barplot(polymer$count)

polymer = polymer[complete.cases(polymer),]
polymer = polymer %>% filter(poly_f != "Not Reported") %>% arrange(desc(count))
polymer = as.data.frame(polymer)

# calcualte proportions
polymer$props = polymer$count/sum(polymer$count)
polymer

# stacked barchart of proportions
polymer2 = matrix((polymer[,3]),ncol = 1)
rownames(polymer2) = (polymer[,1])
colnames(polymer2) = "TomEx2.0"

par(mar = c(3,3,1,20))
barplot(polymer2, legend=rownames(polymer2), col = hcl.colors(30, "RdYlGn"),
        args.legend = list(x = 1.2, cex = 0.8, xjust = 0))

#### How to deal with that many different polymer types? Can we lump some of them together?
### This will be also an issue when trying to compare them to environmental samples, if the polymer names do not match

# size_f-----
# size_class = treatdat %>%  group_by(size_f) %>%summarize(count = n())
# 
# opar = par(mar = c(10,3,1,1))
# barplot(size_class$count, names.arg = size_class$size_f, las =2)
# par(opar)

# Size numeric -------
# boxplot comparing ToMEx2.0 measurements of particle length with data provided by Kooi et al. 2021 on environmental samples

boxplot_comp(treatdat$size.length.um.used.for.conversions, tomex_color = "darkorange2",
             dat.lengths, c("palegreen", "cornsilk3", "lightblue3", "skyblue","lightblue1", "skyblue3"))

#### We can do similar plots for particle width, ratio of width to length
#### + we can add more values from additional studies for which raw data are available


# shape_f ----
shape = treatdat %>%  group_by(shape_f) %>%summarize(count = n())
(shape= as.data.frame(shape))

# remove not reported
shape = shape[-4,]

# calcualte proportions
shape$props = shape$count/sum(shape$count)
shape

# stacked barchart of proportions
shape2 = matrix((shape[,3]),ncol = 1)
rownames(shape2) = (shape[,1])
colnames(shape2) = "TomEx2.0"
barplot(shape2, legend=rownames(shape2), col = c("grey60", "royalblue", "darkolivegreen3"))

#### Add data from environmental samples next to this for comparison.
#### Kooi et al. 2021 unfortunately do not provide their data


### NOT THE REST ####

# charge
(charge = treatdat %>%  group_by(charge) %>%summarize(count = n()))
charge$charge[is.na(charge$charge)] = "not reported"

opar = par(mar = c(10,3,1,1))
barplot(charge$count, names.arg = charge$charge, las =2)
par(opar)

# zeta potential
length(treatdat$zetapotential.mV)
hist(as.numeric(treatdat$zetapotential.mV)) # there are still +- values in the data


# how many datapoints have missing dose measures?
sum(is.na(dat$dose.mg.L.master))
sum(is.na(dat$dose.particles.mL.master))


