
library(dplyr)
library(pals)

# import data
dat = readRDS("Data/prepared_data.RDS")

# prepare polymer groups: All polymer types with only one entry are considered as "others"
others.list = as.data.frame(table(dat$poly_f))$Var1[as.data.frame(table(dat$poly_f))$Freq == 1]
dat$poly_group = as.character(dat$poly_f)
dat$poly_group[dat$poly_f %in% others.list] = "Other"
dat$poly_group[dat$poly_group == "Mix - See Original Study"] = "Polymer mix"
dat$poly_group = as.factor(dat$poly_group)

# split freshwater from marine
fresh = dat[dat$environment == "Freshwater",]
marine = dat[dat$environment == "Marine",]


# Polymer-----

(polyfresh = fresh %>%  group_by(poly_group) %>% summarize(count = n()) %>% 
   mutate(freq = count/sum(count))) 
(polymarine = marine %>%  group_by(poly_group) %>% summarize(count = n()) %>% 
    mutate(freq = count/sum(count)))
polymerged = merge(polyfresh, polymarine, by = "poly_group", sort = FALSE, all.x = TRUE)
(polymerged = polymerged %>% arrange(desc(freq.x), .by_group = TRUE))

polybar = cbind(polymerged$freq.x, polymerged$freq.y)
polybar= as.matrix(t(polybar))
colnames(polybar) = polymerged$poly_group

op = par(mar = c(12,3,1,1))
barplot(polybar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,0.7))
legend("topright", legend = c("Freshwater", "Marine"), 
       fill = c("lightskyblue", "royalblue3"), bty = "n")
par(op)

# Shape ----
(shapefresh = fresh %>%  group_by(shape_f) %>% summarize(count = n()) %>% 
   mutate(freq = count/sum(count))) 
(shapemarine = marine %>%  group_by(shape_f) %>% summarize(count = n()) %>% 
    mutate(freq = count/sum(count)))
shapemerged = merge(shapefresh, shapemarine, by = "shape_f", sort = FALSE)
(shapemerged = shapemerged %>% arrange(desc(freq.x), .by_group = TRUE))

shapebar = cbind(shapemerged$freq.x, shapemerged$freq.y)
shapebar= as.matrix(t(shapebar))
colnames(shapebar) = shapemerged$shape_f

op = par(mar = c(12,3,1,1))
barplot(shapebar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,0.7))
legend("topright", legend = c("Freshwater", "Marine"), 
       fill = c("lightskyblue", "royalblue3"), bty = "n")
par(op)

# Particle length -----
op = par(bty = "l", mfrow = c(1,2))
boxplot(dat$size.length.um.used.for.conversions ~ dat$environment,
        bty = "l", las = 2, xlab = "", pch = "°",
        col = c("lightskyblue", "royalblue3")) #,border = "grey60"
summary(fresh$size.length.um.used.for.conversions)
sd(fresh$size.length.um.used.for.conversions, na.rm = TRUE)
summary(marine$size.length.um.used.for.conversions)
sd(marine$size.length.um.used.for.conversions, na.rm = TRUE)

# Particle width -----
boxplot(dat$size.width.um.used.for.conversions ~ dat$environment,
        bty = "l", las = 2, xlab = "", pch = "°",
        col = c("lightskyblue", "royalblue3"))
summary(fresh$size.width.um.used.for.conversions)
sd(fresh$size.width.um.used.for.conversions, na.rm = TRUE)
summary(marine$size.width.um.used.for.conversions)
sd(marine$size.width.um.used.for.conversions, na.rm = TRUE)

par(op)


# Sodium.azide ----

(sodfresh = fresh %>%  group_by(sodium.azide) %>% summarize(count = n()) %>% 
   mutate(freq = count/sum(count))) 
(sodmarine = marine %>%  group_by(sodium.azide) %>% summarize(count = n()) %>% 
    mutate(freq = count/sum(count)))
sodmerged = merge(sodfresh, sodmarine, by = "sodium.azide", sort = FALSE)
(sodmerged = sodmerged %>% arrange(desc(freq.x), .by_group = TRUE))

sodbar = cbind(sodmerged$freq.x, sodmerged$freq.y)
sodbar= as.matrix(t(sodbar))
colnames(sodbar) = sodmerged$sodium.azide

op = par(mar = c(12,3,1,1), mfrow = c(1,4))
barplot(sodbar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1))
legend("topright", legend = c("Freshwater", "Marine"), 
       fill = c("lightskyblue", "royalblue3"), bty = "n")
par(op)

# DOM_present ----

(domfresh = fresh %>%  group_by(DOM_present) %>% summarize(count = n()) %>% 
   mutate(freq = count/sum(count))) 
(dommarine = marine %>%  group_by(DOM_present) %>% summarize(count = n()) %>% 
    mutate(freq = count/sum(count)))
dommerged = merge(domfresh, dommarine, by = "DOM_present", sort = FALSE)
(dommerged = dommerged %>% arrange(desc(freq.x), .by_group = TRUE))

dombar = cbind(dommerged$freq.x, dommerged$freq.y)
dombar= as.matrix(t(dombar))
colnames(dombar) = dommerged$DOM_present

#op = par(mar = c(12,3,1,1))
barplot(dombar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1))
legend("topright", legend = c("Freshwater", "Marine"), 
       fill = c("lightskyblue", "royalblue3"), bty = "n")
#par(op)

# charge ----
(chargefresh = fresh %>%  group_by(charge) %>% summarize(count = n()) %>% 
   mutate(freq = count/sum(count))) 
(chargemarine = marine %>%  group_by(charge) %>% summarize(count = n()) %>% 
    mutate(freq = count/sum(count)))
chargemerged = merge(chargefresh, chargemarine, by = "charge", sort = FALSE)
(chargemerged = chargemerged %>% arrange(desc(freq.x), .by_group = TRUE))

chargebar = cbind(chargemerged$freq.x, chargemerged$freq.y)
chargebar= as.matrix(t(chargebar))
colnames(chargebar) = chargemerged$charge

#op = par(mar = c(12,3,1,1))
barplot(chargebar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1))
legend("topright", legend = c("Freshwater", "Marine"), 
       fill = c("lightskyblue", "royalblue3"), bty = "n")
#par(op)

# funcional group ----

(functfresh = fresh %>%  group_by(functional.group) %>% summarize(count = n()) %>% 
   mutate(freq = count/sum(count))) 
(functmarine = marine %>%  group_by(functional.group) %>% summarize(count = n()) %>% 
    mutate(freq = count/sum(count)))
functmerged = merge(functfresh, functmarine, by = "functional.group", sort = FALSE)
(functmerged = functmerged %>% arrange(desc(freq.x), .by_group = TRUE))

functbar = cbind(functmerged$freq.x, functmerged$freq.y)
functbar= as.matrix(t(functbar))
colnames(functbar) = functmerged$funct

#op = par(mar = c(12,3,1,1))
barplot(functbar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1))
legend("topright", legend = c("Freshwater", "Marine"), 
       fill = c("lightskyblue", "royalblue3"), bty = "n")
par(op)

## combined descriptive figure

layout(matrix(1:8, ncol = 4, byrow = TRUE))

op = par(mar = c(12,4,1,1), bty = "l")

# A - polymer
barplot(polybar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1), ylab = "Proportion")
legend("topright", legend = c("Freshwater", "Marine"), 
       fill = c("lightskyblue", "royalblue3"), bty = "n")
text("A", x = 2, y = 0.97, cex = 1.5, font = 1)

# B - shape
barplot(shapebar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1), ylab = "Proportion")
text("B", x = 2, y = 0.97, cex = 1.5, font = 1)

# C - length
boxplot(dat$size.length.um.used.for.conversions ~ dat$environment,
        bty = "l", las = 2, xlab = "", pch = "°",
        col = c("lightskyblue", "royalblue3"),
        ylab = "Particle length in µm") 
text("C", x = 0.6, y = 400, cex = 1.5, font = 1)


# D - width
boxplot(dat$size.width.um.used.for.conversions ~ dat$environment,
        bty = "l", las = 2, xlab = "", pch = "°",
        col = c("lightskyblue", "royalblue3"),
        ylab = "Particle width in µm")
text("D", x = 0.6, y = 350, cex = 1.5, font = 1)
