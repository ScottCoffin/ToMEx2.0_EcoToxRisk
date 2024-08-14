library(dplyr)

dat = read.csv("ManiEtAl2016_extracted_polymer_types.csv", header = TRUE, stringsAsFactors = TRUE)

dat %>% group_by(polymer.type) %>% summarise(sum = sum(number_of_particles)) %>% mutate(freq = sum/sum(sum))
