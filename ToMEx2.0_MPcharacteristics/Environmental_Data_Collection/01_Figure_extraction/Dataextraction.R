library(metaDigitise)


data = metaDigitise("Figures_to_extract/", summary = FALSE)


library(dplyr)
get_bar_lengths = function(y) {
  out = rep(NA, length(y))
  for(i in 2:length(y)){
    out[i] = y[i] - y[i-1]
  }
  out
}

x = data$scatterplot$KooiEtAl2021_Fig3.png
class(x$id)
# get heights for the different bars; the first value in each group is the 0 point and is discarded
result = data.frame(x$id, get_bar_lengths(x$y))
result$polymer_type = c(NA, "PVC", "Polyester", "Acrylates,PU,varnish", "Polycaprolactone", "Rupper Type 3", "Other", "Nitrile","PP","PE-chlorinated","PE",
                        NA, "PVC","Polyester","PA","Acrylates...","Polycaprolactone","Rupper Type 3","Other","PP","PE chlorinated","PE",
                        NA, "PVC","Polyester","Polychloropene","PA","Acrylates...","Polycaprolactone","Rubber Type 3","Other","Nitrile","PP","PE chlorinated","PE",
                        NA, "Polyester", "PA", "Acrylates...","Polycaprolactone","Rubber type 3", "Other","Nitrile","Ethylene-vinyl-acetate","PP", "PEchlorinated","PE")
