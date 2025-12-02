detach("package:plyr")
library(dplyr)
library(pals)
getwd()

# import data
dat = readRDS("ToMEx2.0_MPcharacteristics/Data/prepared_data.RDS")

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

# op = par(mar = c(12,3,1,1))
# barplot(polybar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
#         ylim = c(0,0.7))
# legend("topright", legend = c("Freshwater", "Marine"), 
#        fill = c("lightskyblue", "royalblue3"), bty = "n")
# par(op)

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

# op = par(mar = c(12,3,1,1))
# barplot(shapebar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
#         ylim = c(0,0.7))
# legend("topright", legend = c("Freshwater", "Marine"), 
#        fill = c("lightskyblue", "royalblue3"), bty = "n")
# par(op)

# Particle length -----
# op = par(bty = "l", mfrow = c(1,2))
# boxplot(dat$size.length.um.used.for.conversions ~ dat$environment,
#         bty = "l", las = 2, xlab = "", pch = "°",
#         col = c("lightskyblue", "royalblue3")) #,border = "grey60"
summary(fresh$size.length.um.used.for.conversions)
sd(fresh$size.length.um.used.for.conversions, na.rm = TRUE)
summary(marine$size.length.um.used.for.conversions)
sd(marine$size.length.um.used.for.conversions, na.rm = TRUE)

# Particle width -----
# boxplot(dat$size.width.um.used.for.conversions ~ dat$environment,
#         bty = "l", las = 2, xlab = "", pch = "°",
#         col = c("lightskyblue", "royalblue3"))
summary(fresh$size.width.um.used.for.conversions)
sd(fresh$size.width.um.used.for.conversions, na.rm = TRUE)
summary(marine$size.width.um.used.for.conversions)
sd(marine$size.width.um.used.for.conversions, na.rm = TRUE)

# Particle surface area -----
summary(fresh$particle.surface.area.um2)
summary(marine$particle.surface.area.um2)

# Particle density -----
summary(as.numeric(fresh$density.g.cm3))
summary(as.numeric(marine$density.g.cm3))


# Summary statistics entered manually to "Data/data_comp_to_env_compiled.csv" for NMDS comparison with environmental samples



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
# 
# op = par(mar = c(12,3,1,1), mfrow = c(1,4))
# barplot(sodbar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
#         ylim = c(0,1))
# legend("topright", legend = c("Freshwater", "Marine"), 
#        fill = c("lightskyblue", "royalblue3"), bty = "n")
# par(op)

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
# barplot(dombar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
#         ylim = c(0,1))
# legend("topright", legend = c("Freshwater", "Marine"), 
#        fill = c("lightskyblue", "royalblue3"), bty = "n")
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
# barplot(chargebar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
#         ylim = c(0,1))
# legend("topright", legend = c("Freshwater", "Marine"), 
#        fill = c("lightskyblue", "royalblue3"), bty = "n")
# #par(op)

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
# barplot(functbar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
#         ylim = c(0,1))
# legend("topright", legend = c("Freshwater", "Marine"), 
#        fill = c("lightskyblue", "royalblue3"), bty = "n")
# par(op)

## combined descriptive figure

png("ToMEx2.0_MPcharacteristics/Plots/descriptive_plots.png", width = 17, height = 15, units = "cm", res = 1000)

layout(matrix(1:8, ncol = 4, byrow = TRUE))

op = par(mar = c(10,4,2,1), bty = "l", xpd = TRUE)

# A - polymer
barplot(polybar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1), ylab = "", cex.names = 0.8, cex.axis = 0.8)
legend(x = 10, y = 0.6, legend = c("Freshwater", "Marine"), cex = 0.8, 
       fill = c("lightskyblue", "royalblue3"), bty = "n")
mtext("Proportion of particles", cex = 0.5, side = 2, line = 2)
text("A - Polymer type", x = 1, y = 1, cex = 1, font = 1, pos = 4)

# B - shape
barplot(shapebar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1), ylab = "", cex.names = 0.8, cex.axis = 0.8)
mtext("Proportion of particles", cex = 0.5, side = 2, line = 2)
text("B - Particle shape", x = 1, y = 1, cex = 1, font = 1, pos = 4)

# C - length
boxplot(dat$size.length.um.used.for.conversions ~ dat$environment,
        bty = "l", las = 2, xlab = "", pch = "°",
        col = c("lightskyblue", "royalblue3"),
        ylab = "µm") 
text("C - Particle length", x = 0.45, y = 400, cex = 1, font = 1, pos = 4)


# D - width
boxplot(dat$size.width.um.used.for.conversions ~ dat$environment,
        bty = "l", las = 2, xlab = "", pch = "°",
        col = c("lightskyblue", "royalblue3"),
        ylab = "µm")
text("D - Particle width", x = 0.45, y = 350, cex = 1, font = 1, pos = 4)

op = par(mar = c(6,4,1,1))

# E - sodium azide
barplot(sodbar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1), ylab = "Proportion of particles")
text("E - Sodium azide", x = 1, y = 0.96, cex = 1, font = 1, pos = 4)

# F - DOM
barplot(dombar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1), ylab = "Proportion of particles")
text("F - DOM present", x = 2, y = 0.96, cex = 1, font = 1, pos = 4)

# G - Surface charge
barplot(chargebar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1), ylab = "Proportion of particles")
text("G - Surface charge", x = 1, y = 0.96, cex = 1, font = 1, pos = 4)

# H - Surface functionalization
barplot(functbar, beside = TRUE, las = 2, col = c("lightskyblue", "royalblue3"),
        ylim = c(0,1), ylab = "Proportion of particles")
text("H - Functionalization", x = 1, y = 0.96, cex = 1, font = 1, pos = 4)

dev.off()


##### GGPLOT Approach ####
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(forcats) # For reordering factor levels

# Import data
dat <- readRDS("ToMEx2.0_MPcharacteristics/Data/prepared_data.RDS")

# Prepare polymer groups
others.list <- as.data.frame(table(dat$poly_f))$Var1[as.data.frame(table(dat$poly_f))$Freq == 1]
dat$poly_group <- as.character(dat$poly_f)
dat$poly_group[dat$poly_f %in% others.list] <- "Other"
dat$poly_group[dat$poly_group == "Mix - See Original Study"] <- "Polymer mix"
dat$poly_group <- as.factor(dat$poly_group)

# Split freshwater and marine data
fresh <- dat[dat$environment == "Freshwater", ]
marine <- dat[dat$environment == "Marine", ]


# Polymer ----
polyfresh <- fresh %>% group_by(poly_group) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
polymarine <- marine %>% group_by(poly_group) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
polymerged <- merge(polyfresh, polymarine, by = "poly_group", sort = FALSE, all.x = TRUE)
polymerged <- polymerged %>% arrange(desc(freq.x), .by_group = TRUE)

# Reshape to long format
polybar_long <- polymerged %>%
  pivot_longer(cols = c(freq.x, freq.y), names_to = "Environment", values_to = "Proportion") %>%
  mutate(Environment = ifelse(Environment == "freq.x", "Freshwater", "Marine"))

# Reorder x-axis levels by mean proportion
polybar_long <- polybar_long %>%
  group_by(poly_group) %>%
  mutate(mean_proportion = mean(Proportion, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(poly_group = fct_reorder(poly_group, mean_proportion))

# Create the polymer plot
plot_a <- ggplot(data = polybar_long, aes(x = poly_group, y = Proportion, fill = Environment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0,1),
                     n.breaks = 3) +
  scale_fill_manual(values = c("lightskyblue", "royalblue3")) +
  labs(y = "Proportion of particles", title = "A - Polymer type") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

# Particle Shape ----
shapefresh <- fresh %>% group_by(shape_f) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
shapemarine <- marine %>% group_by(shape_f) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
shapemerged <- merge(shapefresh, shapemarine, by = "shape_f", sort = FALSE)
shapemerged <- shapemerged %>% arrange(desc(freq.x), .by_group = TRUE)

# Reshape to long format
shapebar_long <- shapemerged %>%
  pivot_longer(cols = c(freq.x, freq.y), names_to = "Environment", values_to = "Proportion") %>%
  mutate(Environment = ifelse(Environment == "freq.x", "Freshwater", "Marine"))

# Reorder x-axis levels by mean proportion
shapebar_long <- shapebar_long %>%
  group_by(shape_f) %>%
  mutate(mean_proportion = mean(Proportion, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(shape_f = fct_reorder(shape_f, mean_proportion))

# Create the particle shape plot
plot_b <- ggplot(data = shapebar_long, aes(x = shape_f, y = Proportion, fill = Environment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightskyblue", "royalblue3")) +
  scale_y_continuous(limits = c(0,1),
                     n.breaks = 3) +
  labs(y = "Proportion of particles", title = "B - Particle shape") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# Particle length ----
plot_c <- ggplot(data = dat, aes(x = environment, y = size.length.um.used.for.conversions, fill = environment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("lightskyblue", "royalblue3")) +
  scale_y_log10(limits = c(0.001, 1000)) +
  labs(y = "µm", title = "C - Particle length") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title.x = element_blank(),
        legend.position = "none")

# Particle width ----
plot_d <- ggplot(data = dat, aes(x = environment, y = size.width.um.used.for.conversions, fill = environment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("lightskyblue", "royalblue3")) +
  scale_y_log10(limits = c(0.001, 1000)) +
  labs(y = "µm", title = "D - Particle width") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# Sodium azide ----
sodfresh <- fresh %>% group_by(sodium.azide) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
sodmarine <- marine %>% group_by(sodium.azide) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
sodmerged <- merge(sodfresh, sodmarine, by = "sodium.azide", sort = FALSE)

sodbar_long <- sodmerged %>%
  pivot_longer(cols = c(freq.x, freq.y), names_to = "Environment", values_to = "Proportion") %>%
  mutate(Environment = ifelse(Environment == "freq.x", "Freshwater", "Marine"))

plot_e <- ggplot(data = sodbar_long, aes(x = sodium.azide, y = Proportion, fill = Environment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightskyblue", "royalblue3")) +
  scale_y_continuous(limits = c(0,1),
                     n.breaks = 3
                      ) +
  labs(y = "Proportion of particles", title = "E - Sodium azide") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
        )

# DOM ----
domfresh <- fresh %>% group_by(DOM_present) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
dommarine <- marine %>% group_by(DOM_present) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
dommerged <- merge(domfresh, dommarine, by = "DOM_present", sort = FALSE)

dombar_long <- dommerged %>%
  pivot_longer(cols = c(freq.x, freq.y), names_to = "Environment", values_to = "Proportion") %>%
  mutate(Environment = ifelse(Environment == "freq.x", "Freshwater", "Marine"))

plot_f <- ggplot(data = dombar_long, aes(x = DOM_present, y = Proportion, fill = Environment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightskyblue", "royalblue3")) +
  scale_y_continuous(limits = c(0,1),
                     n.breaks = 3) +
  labs(y = "Proportion of particles", title = "F - DOM present") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# Surface charge ----
chargefresh <- fresh %>% group_by(charge) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
chargemarine <- marine %>% group_by(charge) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
chargemerged <- merge(chargefresh, chargemarine, by = "charge", sort = FALSE)

chargebar_long <- chargemerged %>%
  pivot_longer(cols = c(freq.x, freq.y), names_to = "Environment", values_to = "Proportion") %>%
  mutate(Environment = ifelse(Environment == "freq.x", "Freshwater", "Marine"))

plot_g <- ggplot(data = chargebar_long, aes(x = charge, y = Proportion, fill = Environment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightskyblue", "royalblue3")) +
  scale_y_continuous(limits = c(0,1),
                     n.breaks = 3) +
  labs(y = "Proportion of particles", title = "G - Surface charge") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# Functional group ----
functfresh <- fresh %>% group_by(functional.group) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
functmarine <- marine %>% group_by(functional.group) %>% summarize(count = n()) %>% mutate(freq = count / sum(count))
functmerged <- merge(functfresh, functmarine, by = "functional.group", sort = FALSE)

functbar_long <- functmerged %>%
  pivot_longer(cols = c(freq.x, freq.y), names_to = "Environment", values_to = "Proportion") %>%
  mutate(Environment = ifelse(Environment == "freq.x", "Freshwater", "Marine"))

plot_h <- ggplot(data = functbar_long, aes(x = functional.group, y = Proportion, fill = Environment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0,1),
                     n.breaks = 3) +
  scale_fill_manual(values = c("lightskyblue", "royalblue3")) +
  labs(y = "Proportion of particles", title = "H - Functionalization") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())


# Combine all plots into a grid layout with a shared legend
final_plot <- (plot_a | plot_b | plot_c | plot_d) /
  (plot_e | plot_f | plot_g | plot_h) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

final_plot
# Save the final plot
ggsave("ToMEx2.0_MPcharacteristics/Plots/descriptive_plots.png", plot = final_plot, width = 30, height = 20, units = "cm", dpi = 500)


