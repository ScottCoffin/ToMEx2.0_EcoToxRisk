
dat = read.csv("Environmental_Data_Collection/00_Extraction_finished/SoininenEtAl2023/Lapland_RawData_restructured_MM.csv", stringsAsFactors = TRUE)
str(dat)

(lengths = dat %>%  group_by(Sampling_Site) %>% summarise(mean.length = mean(Major.dimension..µm.), 
                                                         median.length = median(Major.dimension..µm.),
                                                         mean.width = mean(Minor.dimension..µm.), 
                                                         median.width = median(Minor.dimension..µm.)))
polymers = dat %>%  group_by(Sampling_Site, Polymer.group) %>% summarise(n = n()) %>%  mutate(prop = n/sum(n))  %>% pivot_wider(names_from = Polymer.group, values_from = prop) %>% 
  summarise(across(everything(), function(x) sum(x, na.rm = TRUE)))

(final = merge(lengths, polymers, by = "Sampling_Site"))                  
final$other = final$abs + final$pmma + final$pan
final
