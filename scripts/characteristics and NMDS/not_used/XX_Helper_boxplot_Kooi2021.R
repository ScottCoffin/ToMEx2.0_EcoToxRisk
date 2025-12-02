#################################
#### NOT USED IN MANUSCRIPT ####
#################################


## Aim: Plot data from Table S3 from Kooi et al. 2021 as boxplots  together with particle characteristics form ToMEx2.0

# tomex_values: measurements of particle characteristics from ToMEx2.0
# ref_dataset: dataset with references data from environmental samples given as: compartment, median, lower.quartile, upper.quartile, average, std; 
#               different groups in different rows
# color: vector of colors used for the boxes of the different groups 

boxplot_comp = function(tomex_values, ref_dataset, tomex_color, ref_color) {
  par(bty = "n")
  max_size = max(c(ref_dataset$average + ref_dataset$std, max(tomex_values, na.rm = T)))
  
  for(i in 1:nrow(ref_dataset)){
    
    label = ref_dataset[i,1]
    quarts = unlist(ref_dataset[i,2:6], use.names = FALSE)
    mean = quarts[4]
    sd_up = mean + quarts[5]
    sd_low = ifelse(mean - quarts[5] > 0, mean - quarts[5], 0)
    mean_tomex = mean(tomex_values, na.rm = T)
    sd_tomex_up = mean_tomex + sd(tomex_values, na.rm = T)
    sd_tomex_low = ifelse(mean_tomex - sd(tomex_values, na.rm = T) > 0, mean_tomex - sd(tomex_values, na.rm = T), 0)
    nplots = nrow(ref_dataset)+1
    bp <- boxplot(0, whisklty=0, staplelty=0, range=1.5, plot=FALSE)
    bp$stats[c(2:4), ] <- as.numeric(quarts[c(2, 1, 3)])
    max_size = max(c(ref_dataset$average + ref_dataset$std, sd_tomex_up))
    if(i == 1){
      par(mar = c(3,0,1,0), oma = c(0,6,0,0), xpd = TRUE)
      layout(matrix(c(1:nplots), nrow = 1, ncol = nplots, byrow = TRUE))
      boxplot(tomex_values, whisklty=0, staplelty=0, outline = FALSE, plot = TRUE, ylim = c(0,max_size + 10), boxfill = tomex_color)
      mtext("Particle length (µm)", side = 2, line = 3, cex = 0.8)
      mtext("ToMEx 2.0", side = 1, at = 1, cex = 0.8) 
      arrows(0.9, mean_tomex, 1.1, mean_tomex, length = 0, lwd = 3, col = "grey40")
      arrows(1, sd_tomex_up,1, sd_tomex_low, length= 0, lwd = 1, col = "grey40")
      abline(h = 0)
      par(mar = c(3,0,1,0))
      bxp(bp, whisklty=0, staplelty=0, boxfill = ref_color[i], outline = FALSE, ylim = c(0, max_size + 10),
          ylab = "", yaxt = "n")
      abline(h = 0)
    }
    if(i != 1) {
      par(mar = c(3,0,1,0))
      bxp(bp, whisklty=0, staplelty=0, boxfill = ref_color[i], outline = FALSE, ylim = c(0, max_size + 10),
          ylab = "", yaxt = "n")
      abline(h = 0)
    }
    mtext(label, side = 1, at = 1, cex = 0.8)
    arrows(0.9,mean,1.1, mean, length = 0, lwd = 3, col = "grey40")
    arrows(1, sd_up, 1, sd_low, length= 0, lwd = 1, col = "grey40")
  }
}



