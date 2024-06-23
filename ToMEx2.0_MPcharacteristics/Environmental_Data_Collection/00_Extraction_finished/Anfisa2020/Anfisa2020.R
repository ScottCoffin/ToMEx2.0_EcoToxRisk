# load packages

library(dplyr)

# import data
anfisa.sub = read.csv("Environmental_Data_Collection/Anfisa2020/chem_sub.csv", stringsAsFactors = TRUE)
anfisa.surf = read.csv("Environmental_Data_Collection/Anfisa2020/chem.csv", stringsAsFactors =  TRUE)

str(anfisa.sub)
str(anfisa.surf)

# Repair numeric variables
anfisa.sub$Long = as.numeric(anfisa.sub$Long)
anfisa.sub$Short[anfisa.sub$Short == "121(30)"] = 121
anfisa.sub$Short = as.numeric(as.character(anfisa.sub$Short))

# remove particles longer than 5 mm from anfisa.sub
anfisa.sub = anfisa.sub[anfisa.sub$Long <= 5000,]

# assign polymer types to polymer groups
# PE, PET, PP, PS, PA, PUR, Acryl, Polyester, PVC, PTFE, Others
levels(anfisa.sub$Polymer)
anfisa.sub$Polymer = as.character(anfisa.sub$Polymer)
anfisa.sub$Polymer[anfisa.sub$Polymer %in% c("NBR","Phenoxy resin", "PMMA", "POM", "PPPO", "SAN")] = "Others" 
anfisa.sub$Polymer = as.factor(anfisa.sub$Polymer)

levels(anfisa.surf$Polymer)

# Overview over all samples - check with numbers reported in publication
mean(anfisa.sub$Square) # 0.04mm^2 (0.04 in publication)
mean(anfisa.sub$Long) #0.73 mm (0.7 in publication)

mean(anfisa.surf$Square) # 3.11 (3.2 in publication)
sd(anfisa.surf$Square) # 4.1 (4.1 in publication)

# calculate size measurements by habitat 
anfisa.sub %>% group_by(wm_new) %>% 
  summarize(mean_length = mean(Long, na.rm = TRUE), median_length = median(Long, na.rm = TRUE),
            mean_width = mean(Short, na.rm = TRUE), median_width = median(Short, na.rm = TRUE),
            mean_area = mean(Square, na.rm = TRUE), mean_density = mean(Dens))

anfisa.surf %>% group_by(wm_new) %>% 
  summarize(mean_area = mean(Square, na.rm = TRUE))

# get proportions of shapes per habitat

(shape.sub = anfisa.sub %>% group_by(wm_new) %>% count(Type))
(shape.sub.sum = shape.sub %>% group_by(wm_new) %>%  summarize(sum = sum(n)))
(shape.sub.merge = merge(shape.sub, shape.sub.sum, by = "wm_new", all.x = TRUE, all.y = FALSE))
shape.sub.merge$prop = shape.sub.merge$n/shape.sub.merge$sum
shape.sub.merge

(shape.surf = anfisa.surf %>% group_by(wm_new) %>% count(Type))
(shape.surf.sum = shape.surf %>% group_by(wm_new) %>%  summarize(sum = sum(n)))
(shape.surf.merge = merge(shape.surf, shape.surf.sum, by = "wm_new", all.x = TRUE, all.y = FALSE))
shape.surf.merge$prop = shape.surf.merge$n/shape.surf.merge$sum
shape.surf.merge

# get proportions of polymers per habitat

(poly.sub = anfisa.sub %>% group_by(wm_new) %>% count(Polymer))
(poly.sub.sum = poly.sub %>% group_by(wm_new) %>%  summarize(sum = sum(n)))
(poly.sub.merge = merge(poly.sub, poly.sub.sum, by = "wm_new", all.x = TRUE, all.y = FALSE))
poly.sub.merge$prop = round(poly.sub.merge$n/poly.sub.merge$sum, digits = 3)
poly.sub.merge


(poly.surf = anfisa.surf %>% group_by(wm_new) %>% count(Polymer))
(poly.surf.sum = poly.surf %>% group_by(wm_new) %>%  summarize(sum = sum(n)))
(poly.surf.merge = merge(poly.surf, poly.surf.sum, by = "wm_new", all.x = TRUE, all.y = FALSE))
poly.surf.merge$prop = round(poly.surf.merge$n/poly.surf.merge$sum, digits = 3)
poly.surf.merge



