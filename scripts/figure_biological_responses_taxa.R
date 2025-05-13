###### Figure: Biological responses for taxa

###  Biological responses for taxa exposed to microplastic exposure ###

library(tidyverse)
library(forcats)
library(scales)

# read in ToMEx 2.o data
ToMExFinal <- readRDS("data/input/aoc_z_tomex2.RDS")


CleanFile <- ToMExFinal %>%
  drop_na(effect_f) %>% 
  filter(effect_f %in% c("No", "Yes")) %>%
  #filter for criteria used in thresholds
  filter(Group != "Bacterium",
           Group != "Plant",
          !env_f %in% c("Terrestrial", "Not Reported"),
         !is.na(effect.metric),
           exp_type_f == "Particle Only",
           effect.metric != "HONEC",
         size.length.um.used.for.conversions >= 1, #1 um and up used
           max.size.ingest.um >= size.length.um.used.for.conversions ,
           tier_zero_tech_f == "Red Criteria Passed",
           tier_zero_risk_f == "Red Criteria Passed", #All thresholds must pass technical and risk red criteria
           risk.13 != 0) %>%  #Drop studies that received a score of 0 for endpoints criteria (this also drops studies that have not yet been))
  mutate(bio_f = recode(bio_f,'subcell'='Subcell','cell' = 'Cell',
                        'tissue'='Tissue',
                        'Organism'='Individual',
                        'Population'='Population'))

group_counts <- CleanFile %>% 
  group_by(Group) %>% 
  summarize(counts = n())


CleanFile <- CleanFile %>% 
  left_join(group_counts) %>% 
  # Reorder Group levels by mean data availability
  mutate(Group = fct_reorder(Group, counts, .fun = mean, na.rm = TRUE))

CleanFile$bio_f <- factor(CleanFile$bio_f, levels=c('Subcell','Cell','Tissue','Individual','Population'))

# Define soft complementary colors
soft_red <- "#D95F02"    # Soft red
soft_green <- "#1B9E77"  # Soft green


Marine = "#1f78b4"   # Dark Blue
Freshwater = "#a6cee3"  # Light Blue

# Create the plot
Figure <- CleanFile %>%
  # Drop empty levels within each facet
  group_by(env_f, bio_f) %>%
  filter(n() > 0) %>%
  ungroup() %>%
  ggplot(aes(y = Group, fill = env_f)) + 
  geom_bar(stat = "count", #need to add 1 to shown 1's
           position = position_stack()) +
  facet_grid(env_f ~ bio_f, 
             labeller = labeller(bio_f = label_value),
             scales = "free_y") +  # Allow y-axis to vary by facet
  labs(x = "Data Points", 
       y = NULL, 
       title = NULL) +
    scale_x_log10(labels = comma_format()) +
   scale_fill_manual(values = c("Freshwater" = Freshwater, "Marine" = Marine)) +
  theme_bw(base_size = 15) + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.background.x = element_rect(fill = "white"),
    strip.background.y = element_rect(fill = "white"),
    strip.text.x = element_text(color = "black", face = "bold"),
    strip.text.y = element_text(color = "black", face = "bold", size = 18),
    text = element_text(color = "black", family = "serif"),
    axis.text = element_text(colour = "black"),
    legend.position = "none",
    panel.spacing = unit(1, "lines")  # Add some space between facets
  ) 

ggsave("output/Manuscript_Figs/bio_response_taxa.jpg",
       Figure,
       width = 7,
       height = 5,
       dpi = 400)

Figure
