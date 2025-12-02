library(dplyr)

env = read.csv("data/input/characteristics and NMDS/data_comp_to_env_compiled.csv", stringsAsFactors = TRUE)
# remove ToMEx datapoints
env = env[-c(1:4),]
# read in data prepared in 01_General_data_prep_tomex.R
tomex = readRDS("data/output/characteristics and NMDS/prepared_data.RDS")

str(dat)
str(tomex)

summary(env)

## Percent PS, PTFE and PE/PET/Polyester etc.
median(env$PE_PET_Polyester, na.rm = TRUE)
median(env$PP_prop, na.rm = TRUE)
median(env$PS_prop, na.rm = TRUE)
median(env$PA_prop, na.rm = TRUE)
median(env$PUR_prop, na.rm = TRUE)
median(env$PVC, na.rm = TRUE)
median(env$PTFE, na.rm = TRUE)

mean(env$PE_PET_Polyester, na.rm = TRUE)
mean(env$PP_prop, na.rm = TRUE)
mean(env$PS_prop, na.rm = TRUE)
mean(env$PA_prop, na.rm = TRUE)
mean(env$PUR_prop, na.rm = TRUE)
mean(env$PVC, na.rm = TRUE)
mean(env$PTFE, na.rm = TRUE)


# TomEx2.0
tomex %>%
  group_by(poly_group) %>%
  dplyr::summarize(count = n()) %>%
  mutate(freq = count / sum(count))

## Percent spheres
# Environment
median(env$spheres_prop, na.rm = TRUE)
range(env$spheres_prop, na.rm = TRUE)
# only freshwater environmental samples
median(env$spheres_prop[env$fresh_marine_binary == "freshwater"], na.rm = TRUE)
# TomEx2.0
tomex %>%
  group_by(shape_f) %>%
  dplyr::summarize(count = n()) %>%
  mutate(freq = count / sum(count))

## Percent fibers
# Environment
median(env$fibers_prop, na.rm = TRUE)
range(env$fibers_prop, na.rm = TRUE)
# only freshwater environmental samples
median(env$fibers_prop[env$fresh_marine_binary == "freshwater"], na.rm = TRUE)
# only marine environmental samples
median(env$fibers_prop[env$fresh_marine_binary == "marine"], na.rm = TRUE)
# TomEx2.0
tomex %>%
  group_by(shape_f) %>%
  dplyr::summarize(count = n()) %>%
  mutate(freq = count / sum(count))


## Particle length
# Environment
summary(env$length_mean_um)
sd(env$length_mean_um, na.rm = TRUE)
# ToMEx 2.0
summary(tomex$size.length.um.used.for.conversions)
sd(tomex$size.length.um.used.for.conversions, na.rm = TRUE)
# plot
boxplot(env$length_mean_um, tomex$size.length.um.used.for.conversions)

## Particle width
# Environment
summary(env$width_mean_um)
sd(env$width_mean_um, na.rm = TRUE)
# ToMEx 2.0
summary(tomex$size.width.um.used.for.conversions)
sd(tomex$size.width.um.used.for.conversions, na.rm = TRUE)
# only freshwater Tomex
summary(tomex$size.width.um.used.for.conversions[
  tomex$environment == "Freshwater"
])

# plot
boxplot(env$width_mean_um, tomex$size.width.um.used.for.conversions)

