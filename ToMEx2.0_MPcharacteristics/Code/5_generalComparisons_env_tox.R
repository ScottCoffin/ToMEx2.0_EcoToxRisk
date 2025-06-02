#### General comparisons (summary tables, etc.) between environmental and toxicity datasets ###

library(tidyverse)
# read in compiled environmental/toxicity summary data
dat = read.csv("ToMEx2.0_MPcharacteristics/Data/data_comp_to_env_compiled.csv", stringsAsFactors = TRUE) %>% 
  mutate(env_tox = case_when(
    doi == "ToMEx2.0" ~ "ToMEx",
    T ~ "Environment"
  ))

# compare relative proportions of polymers in ToMEx 2.0 vs. environment
poly_compare <- dat %>% 
  group_by(fresh_marine_binary, env_tox) %>% 
  summarize(PP = median(PP_prop, na.rm = T),
            PE_PET_Polyester = median(PE_PET_Polyester, na.rm = T),
            PS_min = min(PS_prop, na.rm = T),
            PS_max = max(PS_prop, na.rm = T),
            PA = median(PA_prop, na.rm = T),
            PUR = median(PUR_prop, na.rm = T),
            PVC = median(PVC, na.rm = T),
            PTFE = median(PTFE, na.rm = T),
            other = median(other, na.rm = T)) %>% 
  t() %>% 
  as.data.frame() #%>% 
 # slice(-1:3)

colnames(poly_compare) <- c("Freshwater - Env", "Freshwater - ToMEx",
                            "Marine - Env", "Marine - ToMEx")
poly_compare %>% 
  arrange(desc(`Marine - Env`))


#### Shape Compare
shape_compare <- dat %>% 
  group_by(fresh_marine_binary, env_tox) %>% 
  summarize(fiber_min = min(fibers_prop, na.rm = T),
            fiber_max = max(fibers_prop, na.rm = T),
            fiber_median = median(fibers_prop, na.rm = T))

# compare average particle lengths in ToMEx vs. environment grouped by environment
dat %>% 
  group_by(fresh_marine_binary, env_tox) %>% 
  summarize(min = min(length_median_um, na.rm = T),
            max = max(length_median_um, na.rm = T))
