library(metaDigitise)
library(dplyr)

data = metaDigitise("Figures_to_extract/", summary = FALSE)

get_bar_lengths = function(y) {
  out = rep(NA, length(y))
  for(i in 2:length(y)){
    out[i] = y[i] - y[i-1]
  }
  out
}

# SoltaniEtAl2022 ----

# polymer types surface
x = data$scatterplot$SoltaniEtAl2022_Fig6A.png
x$barheight = get_bar_lengths(x$y)
x[,c(1,3,9)]

# polymer types sediment
x = data$scatterplot$SoltaniEtAl2022_Fig6B.png
x$barheight = get_bar_lengths(x$y)
x[,c(1,3,9)]

# size surface
x = data$scatterplot$SoltaniEtAl2022_FigS3A.png
x$barheight = get_bar_lengths(x$y)
x$land.use = round(x$x,0)
x$mean.group.size = rep(c(1000, 750, 375, 175, 50), 4)
x$barheight[x$mean.group.size == 1000] = x$y[x$mean.group.size == 1000]
x$barheight[x$barheight < 0] = 0

calc.mean.size = function(p) {
  p[1]*1000 + p[2]*750 + p[3]*375 + p[4]*175 + p[5]*50
}

mean.sizes = x %>% group_by(land.use) %>% summarize(mean.size = calc.mean.size(barheight))
mean.sizes 

# size sediment
x = data$scatterplot$SoltaniEtAl2022_FigS3B.png
x$barheight = get_bar_lengths(x$y)
x$land.use = round(x$x,0)
x$mean.group.size = rep(c(1000, 750, 375, 175, 50), 4)
x$barheight[x$mean.group.size == 1000] = x$y[x$mean.group.size == 1000]
x$barheight[x$barheight < 0] = 0

calc.mean.size = function(p) {
  p[1]*1000 + p[2]*750 + p[3]*375 + p[4]*175 + p[5]*50
}

mean.sizes = x %>% group_by(land.use) %>% summarize(mean.size = calc.mean.size(barheight))
mean.sizes 


# KooiEtAl2021 ----
x = data$scatterplot$KooiEtAl2021_Fig3.png
class(x$id)
# get heights for the different bars; the first value in each group is the 0 point and is discarded
result = data.frame(x$id, get_bar_lengths(x$y))
result$polymer_type = c(NA, "PVC", "Polyester", "Acrylates,PU,varnish", "Polycaprolactone", "Rupper Type 3", "Other", "Nitrile","PP","PE-chlorinated","PE",
                        NA, "PVC","Polyester","PA","Acrylates...","Polycaprolactone","Rupper Type 3","Other","PP","PE chlorinated","PE",
                        NA, "PVC","Polyester","Polychloropene","PA","Acrylates...","Polycaprolactone","Rubber Type 3","Other","Nitrile","PP","PE chlorinated","PE",
                        NA, "Polyester", "PA", "Acrylates...","Polycaprolactone","Rubber type 3", "Other","Nitrile","Ethylene-vinyl-acetate","PP", "PEchlorinated","PE")

# DevereuxEtAl2020 ----
## Shapes
x = data$scatterplot$DevereuxEtAl2020_Fig4A.png
x = rbind(0,x)
x$barheight = get_bar_lengths(x$y)
x

## Calculate mean particle lengths
x = data$scatterplot$DevereuxEtAl2020_Fig4C.png
x =x[,c(1,3)]
x$id = as.character(x$id)
x = rbind(c("29-4-5",0),x)
x$date = sapply(strsplit(x$id, "-"), `[[`, 1)
x = x[-39,]
x$average.group.size = as.numeric(rep(c(4.5,3.5, 2.5,1.5,0.75,0.25),9))
x$y = as.numeric(x$y)
x$barheight = get_bar_lengths(x$y)

x$barheight[x$average.group.size == 4.5] = x$y[x$average.group.size == 4.5]
x

table(x$date) # each date has data for all 6 size classes

calc.mean.size = function(p) {
  p[1]*4.5 + p[2]*3.5 + p[3]*2.5 + p[4]*1.5 + p[5]*0.75 + p[6]*0.25
}

mean.sizes = x %>% group_by(date) %>% summarize(mean.size = calc.mean.size(barheight))
mean.sizes 


## Polymer types
x = data$scatterplot$DevereuxEtAl2020_Fig4D.png
x$barheight = get_bar_lengths(x$y)
x[,c(1,3,9)]

# CorcoranEtAl2020 -----

dat1 = data$mean_error$CorcoranEtAl2020_Fig4A.png
dat2 = data$mean_error$CorcoranEtAl2020_Fig4B.png
dat = rbind(dat1, dat2)
dat$n.particles = round(dat$mean, 0)

# calculate proportions of shapes
dat %>% group_by(variable) %>% summarise(sum = sum(n.particles)) %>% mutate(freq = sum/sum(sum))

# calculate proportions of polymer types
dat %>% group_by(id) %>% summarise(sum = sum(n.particles))  %>% mutate(freq = sum/sum(sum))
