
library(xlsx)
library(psych)
library(mc2d)
library(trapezoid)
# Code runs pSSD for number-based concentrations

#### Main parameters ####

# number of simulations for the triangular distributions of the endpoints and assessment factors
SIM <- 10

# coefficient of variation for the endpoint distributions
CV.DP <- 0.3

# coefficient of variation for the assessment factor distributions
CV.UF <- 0.5

## read functions
setwd("scripts")

source("PSSD/do.pssd.r")
source("PSSD/rmore.r")

#import dataset#

#### read input data ####

library(readxl)

getwd()

setwd("PSSD/updated")


# requires species values as endpoint value 

Datapoints <- as.matrix(read_excel("MPs_NoHONECs_mass_V.xlsx", 
                                   sheet = "pSSD_mass"))

#Datapoints_number <- as.matrix(read_excel("MPs_NoHONECs_mass_V.xlsx",
 #                                       sheet = "pSSD_number"))

# requires assessment factors per species endpoint value - AF for time
UFt <- as.matrix(read_excel("MPs_NoHONECs_mass_V.xlsx", 
                            sheet = "UFt"))

# AF for endpoint
UFd <- as.matrix(read_excel("MPs_NoHONECs_mass_V.xlsx", 
                            sheet = "UFd"))

# Polymers used for each endpoint
Polymer <- as.matrix(read_excel("MPs_NoHONECs_mass_V.xlsx", 
                                sheet = "Polymer"))

# Shape of microplastics for each endpoint
Shape <- as.matrix(read_excel("MPs_NoHONECs_mass_V.xlsx", 
                              sheet = "Shape"))

# Dose descriptors of toxicity
DDescriptor <- as.matrix(read_excel("MPs_NoHONECs_mass_V.xlsx", 
                                    sheet = "Dose descriptor"))

# Diameter of the nanoplastics
Diameter <- as.matrix(read_excel("MPs_NoHONECs_mass_V.xlsx", 
                                 sheet = "Diameter"))
Diameterlog <- 3+log10(Diameter)
########################

# run 10 separate pSSDs if necessary, why, this is probably a ram issue?
pSSD1 <- do.pSSD(DP = Datapoints,
                 UFt = UFt,
                 UFdd = UFd,
                 SIM/10, 
                 CV.DP,
                 CV.UF)
pSSD2 <- do.pSSD(DP = Datapoints,
                 UFt = UFt,
                 UFdd = UFd,
                 SIM/10,
                 CV.DP,
                 CV.UF)
pSSD3 <- do.pSSD(DP = Datapoints,
                 UFt = UFt,
                 UFdd = UFd,
                 SIM/10,
                 CV.DP,
                 CV.UF)
pSSD4 <- do.pSSD(DP = Datapoints,
                 UFt = UFt,
                 UFdd = UFd,
                 SIM/10,
                 CV.DP,
                 CV.UF)
pSSD5 <- do.pSSD(DP = Datapoints,
                 UFt = UFt,
                 UFdd = UFd,
                 SIM/10,
                 CV.DP,
                 CV.UF)
pSSD6 <- do.pSSD(DP = Datapoints,
                 UFt = UFt,
                 UFdd = UFd,
                 SIM/10,
                 CV.DP,
                 CV.UF)
pSSD7 <- do.pSSD(DP = Datapoints,
                 UFt = UFt,
                 UFdd = UFd,
                 SIM/10,
                 CV.DP,
                 CV.UF)
pSSD8 <- do.pSSD(DP = Datapoints,
                 UFt = UFt,
                 UFdd = UFd,
                 SIM/10,
                 CV.DP,
                 CV.UF)
pSSD9 <- do.pSSD(DP = Datapoints,
                 UFt = UFt,
                 UFdd = UFd,
                 SIM/10,
                 CV.DP,
                 CV.UF)
pSSD10 <- do.pSSD(DP = Datapoints,
                  UFt = UFt,
                  UFdd = UFd,
                  SIM/10,
                  CV.DP,
                  CV.UF)
pSSD <- cbind(pSSD1, pSSD2, pSSD3, pSSD4, pSSD5, pSSD6, pSSD7, pSSD8, pSSD9, pSSD10)

save(pSSD, file = "pSSD_mass.RData")

#corr.endpoints <- Datapoints/(UFd*UFt)

##### PLOT THE PSSDs ##############################################################################

library(psych) # To use the geometric.mean function

#open a pdf file to start plotting
{ pdf(file = "PSSD_NoHONECs_mass_V.pdf",
      height = 4,
      width  = 5,
      pointsize = 10)
  #some parameters for margins
  par(mfrow = c(1,1), mar = c(4,4.5,3,1), mgp = c(2.5,1,0), xpd = F)
  
  # empty plot with annotations
  plot(c(-4, 8), c(0,1),
       main="PSSD of microplastics in freshwater",
       xlab = "NOEC (ug/L)", ylab = "Cumulative probability",
       type = "n", axes = F)
  axis(1, at = c( -4, -2, 0, 2, 4, 6, 8), 
       labels = expression(10^-4, 10^-2, 1, 10^2, 10^4, 10^6, 10^8), lwd = 0.5)
  axis(2, lwd = 0.5)
  abline(h = seq(0, 1, 0.2), lty = 1, col = "gray90", lwd = 0.5)
  abline(v = seq(-4, 8, 1), lty = 1, col = "gray90", lwd = 0.5)
  box(lwd = 0.5)
  
  # plot some ecdfs
  iv <- seq(-3, 7, 0.001)
  ECDF.data <- matrix(NA, 10000, length(iv))
  
  for(i in 1:10000){
    # calculate the ecdf function
    the.ecdf.f <- ecdf(log(pSSD[,i], base = 10))
    ECDF.data[i,] <- the.ecdf.f(iv)
  }
  
  # to plot quantiles of ECDFs:
  polygon(c(iv,rev(iv)),
          c(apply(ECDF.data,2,function(x){ min(x) }),
            rev(apply(ECDF.data,2,function(x){ max(x) }))),
          border = NA, col = adjustcolor("coral", alpha.f = 0.2))
  
  polygon(c(iv,rev(iv)),
          c(apply(ECDF.data,2,function(x){ quantile(x,0.05) }),
            rev(apply(ECDF.data,2,function(x){ quantile(x,0.95) }))),
          border = NA, col = adjustcolor("coral", alpha.f = 0.5))
  
  polygon(c(iv,rev(iv)),
          c(apply(ECDF.data,2,function(x){ quantile(x,0.25) }),
            rev(apply(ECDF.data,2,function(x){ quantile(x,0.75) }))),
          border = NA, col = adjustcolor("coral"))
  
  
  # calculate and plot the median for each x value of every ecdf pSSD
  # lines(iv, apply(ECDF.data, 2, median), col = "firebrick4", lwd = 2)
  # 
  # # OR calculate and plot the mean for each x value of every ecdf pSSD
  lines(iv, apply(ECDF.data, 2, mean), col = "firebrick4", lwd = 2)
  # 
  # then calculate the geometric mean of the endpoints to show on curve
  
  # deterministic values of NOEC (for plotting)
  NOEC.det <- Datapoints /( UFd * UFt)
  
  # calculate the geometric mean
  NOEC.gmean <- apply(NOEC.det, 2, function(x) geometric.mean(x, na.rm =T ))
  
  # calculate where to put the endpoints on the y axis
  #prop <- 1:ncol(Datapoints)/(ncol(Datapoints)+1)
  prop <- rep(NA, length(NOEC.gmean))
  for (i in 1:length(NOEC.gmean)) {
    prop[i] <- i/(length(NOEC.gmean))-(1/(length(NOEC.gmean))/2)
  }
  
  ind <- rank(NOEC.gmean)
  ind[4] <-4
  ind[12]<-21
  ind[14]<-23
  
  # overlap for two species
  library(stringr)
  
  # whose with additives  
  for (sp in 1:ncol(Datapoints)){
    for (i in 1:nrow(Datapoints)){
      if (str_detect((Shape[i,sp]), "fiber") && !is.na(Shape[i,sp])){
        if (str_detect(Polymer[i,sp], "PET") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 4, cex = 1, col = "green4")
        } else if (str_detect(Polymer[i,sp], "PP") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 4, cex = 1, col = "blue")
        }
      }
      else if (str_detect((Shape[i,sp]), "sphere") && !is.na(Shape[i,sp])){
        if (str_detect(Polymer[i,sp], "PS") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 16, cex = 1, col = "gray18")
        } else if (str_detect(Polymer[i,sp], "PET") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 16, cex = 1, col = "green4") 
        } else if (str_detect(Polymer[i,sp], "PVC") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 16, cex = 1, col = "skyblue")
        } else if (str_detect(Polymer[i,sp], "PE") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 16, cex = 1, col = "darkgray")
        } else if (str_detect(Polymer[i,sp], "PLA") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 16, cex = 1, col = "lightgreen")  
        } else if (str_detect(Polymer[i,sp], "PA") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 16, cex = 1, col = "red")
        } else if (str_detect(Polymer[i,sp], "PUR") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 16, cex = 1, col = "yellow")
        } else if (str_detect(Polymer[i,sp], "PP") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 16, cex = 1, col = "blue")  
        } else if (str_detect(Polymer[i,sp], "AF") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 16, cex = 1, col = "mediumpurple")
        }
      }
      else if (str_detect((Shape[i,sp]), "fragment") && !is.na(Shape[i,sp])){
        if (str_detect(Polymer[i,sp], "PS") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 15, cex = 1, col = "gray18")
        } else if (str_detect(Polymer[i,sp], "PET") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 15, cex = 1, col = "green4") 
        } else if (str_detect(Polymer[i,sp], "PVC") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 15, cex = 1, col = "skyblue")
        } else if (str_detect(Polymer[i,sp], "PE") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 15, cex = 1, col = "darkgray")
        } else if (str_detect(Polymer[i,sp], "PLA") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 15, cex = 1, col = "lightgreen")  
        } else if (str_detect(Polymer[i,sp], "PA") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 15, cex = 1, col = "red")
        } else if (str_detect(Polymer[i,sp], "PUR") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 15, cex = 1, col = "yellow")  
        } else if (str_detect(Polymer[i,sp], "PP") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 15, cex = 1, col = "blue")  
        } else if (str_detect(Polymer[i,sp], "AF") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 15, cex = 1, col = "mediumpurple")
        }
      }
      else if (str_detect((Shape[i,sp]), "irregular") && !is.na(Shape[i,sp])){
        if (str_detect(Polymer[i,sp], "PS") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 17, cex = 1, col = "gray18")
        } else if (str_detect(Polymer[i,sp], "PET") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 17, cex = 1, col = "green4") 
        } else if (str_detect(Polymer[i,sp], "PVC") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 17, cex = 1, col = "skyblue")
        } else if (str_detect(Polymer[i,sp], "PE") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 17, cex = 1, col = "darkgray")
        } else if (str_detect(Polymer[i,sp], "PLA") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 17, cex = 1, col = "lightgreen")
        } else if (str_detect(Polymer[i,sp], "PA") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 17, cex = 1, col = "red")
        } else if (str_detect(Polymer[i,sp], "PUR") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 17, cex = 1, col = "yellow")
        } else if (str_detect(Polymer[i,sp], "PP") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 17, cex = 1, col = "blue")
        } else if (str_detect(Polymer[i,sp], "AF") == T){
          points(x = log10(NOEC.det[i,sp]),
                 y = prop[ind[sp]],
                 pch = 17, cex = 1, col = "mediumpurple")
        }
      }
      
      xoffset <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
      
      text(x = min(log10((NOEC.det[,sp])), na.rm = T) - xoffset[sp],
           adj = 1,
           y = prop[ind[sp]],
           labels = colnames(Datapoints)[sp],
           cex = 0.6, font = 3)
    } 
  }
  
  
  #legend("bottomright", c("Mean PSSD", "PMMA", "PHB", "PS", "PS-COOH", "PS-NH2","with Sodium Azide"),
  #pch = c(NA, 16, 16, 16, 16, 16, 2), col = c("firebrick4", "yellow", "forestgreen", "black","red", "blue","black"),
  #lty = c(1, NA, NA, NA, NA, NA, NA), cex = 0.7, bg = "white", box.lwd = 0.4)
  
  ### NEW PLOT FOR THE PNEC ON NEW PDF PAGE
  
  # calculate it
  # activate when running for the first time
  PNEC1 <- apply(pSSD, 2, function(x) quantile(x, probs = 0.05, type = 1))
  hist(PNEC1, freq = F, breaks = 100, col = rgb(255, 127, 80, max = 255, alpha = 50), border = "coral",
       xlab = "PNEC (ug/L)", ylab = "Probability density", main = "Probability density of the PNEC", lwd = 0.5)
  #lines(density(log(PNEC,10)), lwd = 2, col = "red")
  box(lwd = 0.5)
  
  
  #activate when running for the first time
    save(PNEC1, file = "PNEC_mass.RData")
  
  #close the graph
  dev.off()
}

# PNEC value
Mode_Y	<- function(x)
{
  dens	<- density(x)
  ind		<- which(dens$y==max(dens$y))
  dens$x[ind]
}

Stat_PNEC <- matrix(NA, 2, 8)

Stat_PNEC[1,] <- c("Min", "Q5", "Q25", "Mean", "Mode", "Q75", "Q95", "Max")

Stat_PNEC[2,1] <- min(PNEC1)
Stat_PNEC[2,2] <- quantile(PNEC1, 0.05)
Stat_PNEC[2,3] <- quantile(PNEC1, 0.25)
Stat_PNEC[2,4] <- mean(PNEC1)
Stat_PNEC[2,5] <- Mode_Y(PNEC1)
Stat_PNEC[2,6] <- quantile(PNEC1, 0.75)
Stat_PNEC[2,7] <- quantile(PNEC1, 0.95)
Stat_PNEC[2,8] <- max(PNEC1)

write.csv(Stat_PNEC, file = "Stats_PNEC_NoHONECs_V.csv")

Stat_PNEC
