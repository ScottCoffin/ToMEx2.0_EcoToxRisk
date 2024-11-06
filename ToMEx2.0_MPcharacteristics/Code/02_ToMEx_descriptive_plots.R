
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


# We need summary stats and plot:

# - charge






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

# Particle width -----
boxplot(dat$size.width.um.used.for.conversions ~ dat$environment,
        bty = "l", las = 2, xlab = "", pch = "°",
        col = c("lightskyblue", "royalblue3"))

par(op)


# Sodium.azide ----



# DOM_present ----


# charge ----


# funcional group ----



