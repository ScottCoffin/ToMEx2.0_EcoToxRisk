# Load packages ---
library(vegan)

# Import data ---

dat = read.csv("Data/data_comp_to_env_compiled.csv", stringsAsFactors = TRUE)
str(dat)
summary(dat)
levels(dat$doi)

#dat = dat[!grepl("Kooi", dat$source), ]
# dat$tomex_binary = NA
# dat$tomex_binary[grepl("ToMEx", dat$source)] = "yes"
# dat$tomex_binary[!grepl("ToMEx", dat$source)] = "no"

# scale variables
colnames(dat)
vars = dat[,c(5:17,19)]
scaled = scale(vars)

# replace all NAs with 0 - set all unknown values to variable mean
scaled[is.na(scaled)] = 0

# NMDS
NMDS = metaMDS(scaled, distance="euclidean", k=2, autotransform=F, wascores=F)
NMDS$points
x <- NMDS$points[,1]      
y <- NMDS$points[,2] 

# get vectors for each variable
envfitall <- envfit(NMDS, scaled) 
vectors = as_tibble(envfitall$vectors$arrows, rownames = "properties")
vectors$vectornames = c("PE/PET/Polyester", "PP", "PS", "PA", "PU", "Acrylates", "PVC/PVA", "PTFE", "Other",
                        "Fragments", "Spheres", "Fibers", "Mean length", "Mean width")

adjust.x = c(0.4, 0.1, -0.1, 0, 0.1, 0.1, 0.1, 0, 0, -0.4, -0.3, 0, 0, -0.4)
adjust.y = c(-0.15, 0, 0, 0.1, 0, 0.1, 0.1, -0.1, 0.2, 0, 0.05, 0.1, 0.1, 0)

palette(c("lightskyblue", "royalblue3"))
op = par(bty = "l", las = 1)

plot(x,y, bg = dat$fresh_marine_binary, col = c(rep("red",2), rep("grey20",(nrow(dat)-2))), pch = c(rep(24,2), rep(21,(nrow(dat)-2))), #c(rep("red",2), rep("grey80",(nrow(dat)-2)))
     xlab = "NMDS1", ylab = "NMDS2", cex = 1.5)

for(i in 1:nrow(envfitall$vectors$arrows)){
arrows(0,0, envfitall$vectors$arrows[i,1]*2.5, envfitall$vectors$arrows[i,2]*2.5, 
       length = 0.1, col = "grey50")
text(envfitall$vectors$arrows[i,1]*2.5 + adjust.x[i], envfitall$vectors$arrows[i,2]*2.5 + adjust.y[i],
     vectors$vectornames[i], cex = 0.8, col = "grey50")
}          
text("ToMEx 2.0", x = -2.9, y = -1.05, col = "red", 
     cex = 0.8, pos = 4, font = 2)


### different symbols for different publication sources
plot(x,y, bg = c(rep("red",2), rep("grey80",(nrow(dat)-2))), col = dat$fresh_marine_binary, pch = c(21,23,22,24,15,25)[dat$doi],
     xlab = "NMDS1", ylab = "NMDS2", cex = 1.5)



