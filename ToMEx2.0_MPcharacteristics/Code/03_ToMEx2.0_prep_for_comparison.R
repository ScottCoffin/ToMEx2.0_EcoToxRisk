# load packages
detach("package:plyr")
library(dplyr)

# import data
dat = readRDS("ToMEx2.0_MPcharacteristics/Data/prepared_data.RDS")

# subset variables useful for this analysis
dat = dat %>% select(environment, 
               poly_f, shape_f, particle.surface.area.um2, size.length.um.used.for.conversions, size.width.um.used.for.conversions,
               density.g.cm3)

# aggregate polymer levels
levels(dat$poly_f)
dat$poly_f = as.character(dat$poly_f)
dat$poly_f[dat$poly_f %in% c("Biopolymer", "Not Reported", "Polycarbonate", "Polyethylene vinyl acetate",
                             "Polyisoprene", "Polylactic Acid", "Polymethylmethacrylate", "Mix - See Original Study",
                             "Poly(Styrene-co-acrylonitrile)", "Polyamidoamine", "Polyoxymethylene", 
                             "Polyvinyl Acetate", "Polyvinylchloride/vinylacetate co-polymer", "Sodium Polyacrylate",
                             "Starch/Polybutylene Adipate Terephthalate/Polylactic Acid", "Tire Wear")] = "Others"
dat$poly_f[dat$poly_f %in% c("High Density Polyethylene", "Low Density Polyethylene", "Medium Density Polyethylene", 
                             "Polyethylene Terephthalate")] = "Polyethylene"
dat$poly_f = as.factor(dat$poly_f)
levels(dat$poly_f)
 

# split freshwater from marine
fresh = dat[dat$environment == "Freshwater",]
marine = dat[dat$environment == "Marine",]


# prepare summary stats for comparison of MP characteristics with environmental samples

# rel prop of different polymer types
# aggregate groups to following categories: 

# mean length
# median length
# mean width
# median width
# mean density


fresh$density.g.cm3 = as.numeric(fresh$density.g.cm3)
summary(fresh)
marine$density.g.cm3 = as.numeric(marine$density.g.cm3)
summary(marine)

# polymer type proportion
poly.fresh  = fresh %>% count(poly_f)
poly.fresh$prop = round(poly.fresh$n/sum(poly.fresh$n), digits = 3)
print(poly.fresh, n =100)

poly.marine  = marine %>% count(poly_f)
poly.marine$prop = round(poly.marine$n/sum(poly.marine$n), digits = 3)
print(poly.marine, n =100)

# one hot shape proportion

shape.fresh  = fresh %>% count(shape_f)
shape.fresh$prop = round(shape.fresh$n/sum(shape.fresh$n), digits = 3)
shape.fresh

shape.marine  = marine %>% count(shape_f)
shape.marine$prop = round(shape.marine$n/sum(shape.marine$n), digits = 3)
shape.marine

### Summaries for quality filtered  -----

dat = readRDS("ToMEx2.0_MPcharacteristics/Data/prepared_data.RDS")

# Select only high quality data
dat = dat[dat$tier_zero_risk_f == "Red Criteria Passed",]
dat = dat[dat$tier_zero_tech_f == "Red Criteria Passed",]

# subset variables useful for this analysis
dat = dat %>% select(environment, 
                     poly_f, shape_f, particle.surface.area.um2, size.length.um.used.for.conversions, size.width.um.used.for.conversions,
                     density.g.cm3)

# aggregate polymer levels
levels(dat$poly_f)
dat$poly_f = as.character(dat$poly_f)
dat$poly_f[dat$poly_f %in% c("Biopolymer", "Not Reported", "Polycarbonate", "Polyethylene vinyl acetate",
                             "Polyisoprene", "Polylactic Acid", "Polymethylmethacrylate", "Mix - See Original Study",
                             "Poly(Styrene-co-acrylonitrile)", "Polyamidoamine", "Polyoxymethylene", 
                             "Polyvinyl Acetate", "Polyvinylchloride/vinylacetate co-polymer", "Sodium Polyacrylate",
                             "Starch/Polybutylene Adipate Terephthalate/Polylactic Acid", "Tire Wear")] = "Others"
dat$poly_f[dat$poly_f %in% c("High Density Polyethylene", "Low Density Polyethylene", "Medium Density Polyethylene", 
                             "Polyethylene Terephthalate")] = "Polyethylene"
dat$poly_f = as.factor(dat$poly_f)
levels(dat$poly_f)


# split freshwater from marine
fresh = dat[dat$environment == "Freshwater",]
marine = dat[dat$environment == "Marine",]


# prepare summary stats for comparison of MP characteristics with environmental samples

# rel prop of different polymer types
# aggregate groups to following categories: 

# mean length
# median length
# mean width
# median width
# mean density


fresh$density.g.cm3 = as.numeric(fresh$density.g.cm3)
summary(fresh)
marine$density.g.cm3 = as.numeric(marine$density.g.cm3)
summary(marine)

# polymer type proportion
poly.fresh  = fresh %>% count(poly_f)
poly.fresh$prop = round(poly.fresh$n/sum(poly.fresh$n), digits = 3)
print(poly.fresh, n =100)

poly.marine  = marine %>% count(poly_f)
poly.marine$prop = round(poly.marine$n/sum(poly.marine$n), digits = 3)
print(poly.marine, n =100)

# one hot shape proportion

shape.fresh  = fresh %>% count(shape_f)
shape.fresh$prop = round(shape.fresh$n/sum(shape.fresh$n), digits = 3)
shape.fresh

shape.marine  = marine %>% count(shape_f)
shape.marine$prop = round(shape.marine$n/sum(shape.marine$n), digits = 3)
shape.marine
