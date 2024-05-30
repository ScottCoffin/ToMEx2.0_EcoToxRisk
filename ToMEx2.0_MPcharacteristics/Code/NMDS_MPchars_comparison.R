# Load packages ---
library(vegan)

# Import data ---

dat = read.csv("Data/data_comp_to_env_compiled.csv", stringsAsFactors = TRUE)
str(dat)
summary(dat)
dat = dat[!grepl("Kooi", dat$source), ]
dat$tomex_binary = NA
dat$tomex_binary[grepl("ToMEx", dat$source)] = "yes"
dat$tomex_binary[!grepl("ToMEx", dat$source)] = "no"

# scale variables
colnames(dat)
vars = dat[,c(4:19,22,23)]
scaled = scale(vars)

# raplce all NAs with 0 - set all unknown values to variable mean
scaled[is.na(scaled)] = 0

# NMDS
NMDS = metaMDS(scaled, distance="euclidean", k=2, autotransform=F, wascores=F)
NMDS$points
x <- NMDS$points[,1]      
y <- NMDS$points[,2] 

# get vectors for each variable
envfitall <- envfit(NMDS, scaled) 

palette(c("lightskyblue", "royalblue3"))
op = par(bty = "l", las = 1)
plot(x,y, bg = dat$fresh_marine_binary, col = c(rep("darkred",2), rep("black",7)), pch = ifelse(dat$tomex_binary == "yes", 24, 21),
     xlab = "NMDS1", ylab = "NMDS2")

for(i in 1:nrow(envfitall$vectors$arrows)){
arrows(0,0, envfitall$vectors$arrows[i,1]*2, envfitall$vectors$arrows[i,2]*2, 
       length = 0.1, col = "grey50")
text(envfitall$vectors$arrows[i,1]*2+0.5, envfitall$vectors$arrows[i,2]*2+0.5,
     rownames(envfitall$vectors$arrows)[i], cex = 0.8, col = "grey50")
}          
text("ToMEx 2.0", x = -4.2, y = -0.8, col = "darkred", 
     cex = 0.8, pos = 4)

