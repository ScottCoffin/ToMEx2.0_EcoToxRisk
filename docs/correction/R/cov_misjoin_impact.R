# Impact of the mis-keyed merge of record-level Monte Carlo summaries
# (coefficients of variation, CoV) in the published analysis.
#
# EcoTox_MonteCarlo.Rmd numbered the MC results 1..n over the re-aligned
# records, while aoc_aligned_test keeps the input row number, and merged the
# two by unique_id. This script recomputes the CoVs from the published MC
# output exactly as the Rmd does, then joins them (a) as published and (b) to
# the correct record (the k-th MC record is the k-th row of aoc_aligned_test;
# verified by species), and compares the CoV statistics reported in the paper.
#
# Inputs: the published MC outputs (December 2025), archived at
#   https://doi.org/10.5281/zenodo.23582905 (folder monte_carlo_output_published/); place them in
#   scripts/monte carlo/output/pre_issue6/{results_df_sobol,aoc_MC_sobol}.rds
# Usage (repo root): Rscript docs/correction/R/cov_misjoin_impact.R [out.rds]

suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "docs/correction/results/cov_misjoin.rds"
dir_pub <- "scripts/monte carlo/output/pre_issue6"
stopifnot(file.exists(file.path(dir_pub, "results_df_sobol.rds")))

res <- readRDS(file.path(dir_pub, "results_df_sobol.rds"))
pub <- readRDS(file.path(dir_pub, "aoc_MC_sobol.rds"))
cov_cols <- c(grep("_ox_stress$|_food_dilution$", names(pub), value = TRUE))
cov_cols <- setdiff(cov_cols, c("particles_L_food_dilution_deterministic", "particles_L_ox_stress_deterministic"))

# CoV per MC record, as in EcoTox_MonteCarlo.Rmd ("Determine Coefficient of Variation")
cov_tt <- res %>%
  filter(translocatable %in% c("translocatable", "translocatable (all)", "translocatable (some)"),
         particles_L_ox_stress > 0) %>% tidyr::drop_na(particles_L_ox_stress) %>%
  group_by(unique_id) %>%
  summarise(sd_particles_L_ox_stress = sd(particles_L_ox_stress), mean_particles_L_ox_stress = mean(particles_L_ox_stress),
            median_particles_L_ox_stress = median(particles_L_ox_stress),
            CoV_ox_stress = sd_particles_L_ox_stress / mean_particles_L_ox_stress, .groups = "drop")
cov_fd <- res %>%
  filter(ingestible %in% c("ingestible", "ingestible (all)", "ingestible (some)"), group != "Algae",
         particles_L_food_dilution > 0) %>% tidyr::drop_na(particles_L_food_dilution) %>%
  group_by(unique_id) %>%
  summarise(sd_particles_L_food_dilution = sd(particles_L_food_dilution),
            mean_particles_L_food_dilution = mean(particles_L_food_dilution),
            median_particles_L_food_dilution = median(particles_L_food_dilution),
            CoV_food_dilution = sd_particles_L_food_dilution / mean_particles_L_food_dilution, .groups = "drop")
cov_rec <- full_join(cov_tt, cov_fd, by = "unique_id") %>% mutate(mc_row = as.integer(as.character(unique_id)))

# the records themselves (aoc_aligned_test), without the merged MC columns
records <- pub %>% select(-any_of(setdiff(names(cov_rec), c("unique_id", "mc_row")))) %>%
  mutate(unique_id = as.integer(as.character(unique_id))) %>% arrange(unique_id) %>%
  mutate(mc_row = row_number())

# verify the mapping: the species of MC record k is the species of the k-th record
sp_mc <- res %>% distinct(unique_id, species) %>% mutate(mc_row = as.integer(as.character(unique_id)))
chk <- inner_join(sp_mc, records %>% select(mc_row, Species), by = "mc_row")
mapping_ok <- all(as.character(chk$species) == as.character(chk$Species))
published_key_ok <- inner_join(sp_mc %>% transmute(unique_id = mc_row, species),
                               records %>% select(unique_id, Species), by = "unique_id")
stopifnot(mapping_ok)

joined <- list(
  published = records %>% left_join(cov_rec %>% select(-unique_id) %>% rename(unique_id = mc_row), by = "unique_id"),
  corrected = records %>% left_join(cov_rec %>% select(-unique_id), by = "mc_row")
)

summ <- function(aoc_MC, label) {
  fd <- aoc_MC %>% filter(ingestible %in% c("ingestible", "ingestible (all)", "ingestible (some)"), Group != "Algae",
                          particles_L_food_dilution_deterministic > 0) %>%
    group_by(environment) %>% summarise(ERM = "Food Dilution", n = sum(!is.na(CoV_food_dilution)),
      min = min(CoV_food_dilution, na.rm = TRUE), median = median(CoV_food_dilution, na.rm = TRUE),
      max = max(CoV_food_dilution, na.rm = TRUE), .groups = "drop")
  tt <- aoc_MC %>% filter(translocatable %in% c("translocatable", "translocatable (all)", "translocatable (some)"),
                          particles_L_ox_stress_deterministic > 0) %>%
    group_by(environment) %>% summarise(ERM = "Tissue Translocation", n = sum(!is.na(CoV_ox_stress)),
      min = min(CoV_ox_stress, na.rm = TRUE), median = median(CoV_ox_stress, na.rm = TRUE),
      max = max(CoV_ox_stress, na.rm = TRUE), .groups = "drop")
  sp <- function(y) { t <- suppressWarnings(cor.test(y, aoc_MC$max.size.ingest.um, method = "spearman"))
    tibble(rho = unname(t$estimate), p = t$p.value) }
  list(summary = bind_rows(fd, tt) %>% mutate(join = label),
       spearman = bind_rows(sp(aoc_MC$CoV_food_dilution) %>% mutate(ERM = "Food Dilution"),
                            sp(aoc_MC$CoV_ox_stress) %>% mutate(ERM = "Tissue Translocation")) %>% mutate(join = label),
       range_all = tibble(join = label, min = min(c(aoc_MC$CoV_food_dilution, aoc_MC$CoV_ox_stress), na.rm = TRUE),
                          max = max(c(aoc_MC$CoV_food_dilution, aoc_MC$CoV_ox_stress), na.rm = TRUE)))
}
s_pub <- summ(joined$published, "published (mis-keyed)")
s_cor <- summ(joined$corrected, "corrected")
# corrected simulations (v0.4.0 MC: E4/E5 fixed, merge keyed correctly), if run
cur_file <- "scripts/monte carlo/output/aoc_MC_sobol.rds"
s_new <- if (file.exists(cur_file)) summ(readRDS(cur_file), "corrected simulations (v0.4.0)") else NULL

result <- list(
  n_records = nrow(records), n_mc_records = nrow(cov_rec),
  published_join_matched = nrow(published_key_ok),
  published_join_wrong_species = sum(as.character(published_key_ok$species) != as.character(published_key_ok$Species)),
  published_join_unmatched = nrow(records) - nrow(published_key_ok),
  summary = bind_rows(s_pub$summary, s_cor$summary, s_new$summary),
  spearman = bind_rows(s_pub$spearman, s_cor$spearman, s_new$spearman),
  range_all = bind_rows(s_pub$range_all, s_cor$range_all, s_new$range_all)
)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, out)
print(result[c("n_records", "published_join_matched", "published_join_wrong_species", "published_join_unmatched")])
print(as.data.frame(result$range_all), digits = 3)
print(as.data.frame(result$spearman), digits = 3)
print(as.data.frame(result$summary), digits = 3)
