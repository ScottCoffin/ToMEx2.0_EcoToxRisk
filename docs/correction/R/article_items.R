# Published vs corrected values of the quantitative statements, tables and figures of
# Coffin et al. (2026), J. Hazard. Mater. 503:141021, for the item-by-item list of the
# correction (docs/correction/correction_PSSDpp.qmd, "Changes to the published article").
#
# Published values: Table 2 as printed (typed below), and the analysis outputs of the
# published analysis (git 7163736, 2 December 2025: output/data/supplemental_information.xlsx),
# which reproduce the printed SSD and MC + SSD values of Table 2 and Table S5. Corrected
# values: the current production outputs (output/data, scripts/monte carlo/output,
# data/output). Run from the repository root:
#   Rscript docs/correction/R/article_items.R
# Writes docs/correction/results/article_items.rds.

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
stopifnot(file.exists("ToMEx2.0_EcoToxRisk.Rproj"))
out <- "docs/correction/results/article_items.rds"

# ---- Table 2 as printed (medians and 95 % CIs; NA = not given) ----
pub <- tribble(
  ~Environment, ~ERM, ~Tier, ~Method, ~median, ~lower, ~upper,
  "Marine", "FD", 1, "SSD", 4.8e-5, NA, NA,      "Marine", "FD", 1, "MC", 1.0e-5, 1.0e-5, 1.0e-4,
  "Marine", "FD", 2, "SSD", 2.6e-3, 4.8e-5, 1.1, "Marine", "FD", 2, "MC", 2.6e-3, 2.3e-4, 2.2e-2,
  "Marine", "FD", 3, "SSD", 0.30, 0.03, 7.3,     "Marine", "FD", 3, "MC", 0.28, 0.029, 2.2,
  "Marine", "FD", 3, "PSSD", 0.140, 0.01, 1.2,
  "Marine", "FD", 4, "SSD", 1.5, 0.22, 44,       "Marine", "FD", 4, "MC", 1.3, 0.15, 9.5,
  "Marine", "FD", 4, "PSSD", 0.64, 0.081, 5.0,
  "Marine", "TT", 1, "SSD", 0.82, NA, NA,        "Marine", "TT", 1, "MC", 0.81, 0.21, 1.8,
  "Marine", "TT", 2, "SSD", 13, 0.82, 380,       "Marine", "TT", 2, "MC", 13, 6.4, 30,
  "Marine", "TT", 3, "SSD", 490, 140, 5300,      "Marine", "TT", 3, "MC", 480, 160, 1200,
  "Marine", "TT", 3, "PSSD", 320, 74, 1400,
  "Marine", "TT", 4, "SSD", 1600, 430, 16000,    "Marine", "TT", 4, "MC", 1550, 605, 3820,
  "Marine", "TT", 4, "PSSD", 710, 190, 2500,
  "Freshwater", "FD", 1, "SSD", 0.4, NA, NA,     "Freshwater", "FD", 1, "MC", 0.4, 0.03, 4,
  "Freshwater", "FD", 2, "SSD", 4.9, 0.4, 150,   "Freshwater", "FD", 2, "MC", 4.1, 0.34, 38,
  "Freshwater", "FD", 3, "SSD", 42, 28, 2400,    "Freshwater", "FD", 3, "MC", 36, 3.4, 250,
  "Freshwater", "FD", 3, "PSSD", 31, 0.14, 1500,
  "Freshwater", "FD", 4, "SSD", 140, 63, 6700,   "Freshwater", "FD", 4, "MC", 120, 14, 770,
  "Freshwater", "FD", 4, "PSSD", 140, 6.3, 1800,
  "Freshwater", "TT", 1, "SSD", 93, NA, NA,      "Freshwater", "TT", 1, "MC", 100, 66, 180,
  "Freshwater", "TT", 2, "SSD", 720, 93, 23000,  "Freshwater", "TT", 2, "MC", 710, 460, 1100,
  "Freshwater", "TT", 3, "SSD", 4800, 1900, 54000, "Freshwater", "TT", 3, "MC", 4600, 2000, 6700,
  "Freshwater", "TT", 3, "PSSD", 2150, 167, 38000,
  "Freshwater", "TT", 4, "SSD", 9500, 3200, 120000, "Freshwater", "TT", 4, "MC", 9500, 4100, 14000,
  "Freshwater", "TT", 4, "PSSD", 6400, 420, 56000)

# ---- corrected Table 2 from the production outputs ----
ssd <- read.csv("output/data/all_thresholds_combined_quality.csv", check.names = FALSE, row.names = 1) %>%
  filter(label == "risk and tech", Environment %in% c("Freshwater", "Marine")) %>%
  transmute(Environment, Tier = as.integer(sub("Tier", "", Tier)),
            FD_median = `Food Dilution (Default)`, FD_lower = `Food Dilution (5th %)`, FD_upper = `Food Dilution (95th %)`,
            TT_median = `Tissue Translocation (Default)`, TT_lower = `Tissue Translocation (5th %)`,
            TT_upper = `Tissue Translocation (95th %)`) %>%
  pivot_longer(-c(Environment, Tier), names_to = c("ERM", ".value"), names_sep = "_") %>%
  mutate(Method = "SSD", lower = ifelse(Tier == 1, NA, lower), upper = ifelse(Tier == 1, NA, upper))
mc <- readRDS("scripts/monte carlo/output/summary_stats_base_thresholds_sobol.rds") %>%
  filter(Metric %in% c("Food Dilution (Default)", "Tissue Translocation (Default)")) %>%
  transmute(Environment, Tier = as.integer(sub("Tier", "", Tier)), ERM = ifelse(grepl("Food", Metric), "FD", "TT"),
            Method = "MC", median = Median, lower = Percentile_5th, upper = Percentile_95th)
ps <- readRDS("data/output/PNEC_summary_table_wide.rds") %>%
  transmute(Environment, Tier = as.integer(sub("Tier", "", Tier)), ERM = ifelse(grepl("Food", ERM), "FD", "TT"),
            Method = "PSSD", median = Median, lower = Q5, upper = Q95)
corr <- bind_rows(ssd, mc, ps)
table2 <- full_join(pub, corr, by = c("Environment", "ERM", "Tier", "Method"), suffix = c(".pub", ".corr")) %>%
  mutate(rci.pub = (upper.pub - lower.pub) / median.pub, rci.corr = (upper.corr - lower.corr) / median.corr) %>%
  arrange(desc(Environment), ERM, Tier, factor(Method, c("SSD", "MC", "PSSD")))

# ---- statements derived from Table 2 (same rules for published and corrected) ----
statements <- function(t, v) {
  t <- t %>% mutate(m = .data[[paste0("median.", v)]], lo = .data[[paste0("lower.", v)]],
                    rci = .data[[paste0("rci.", v)]])
  w <- t %>% select(Environment, ERM, Tier, Method, m) %>% pivot_wider(names_from = Method, values_from = m)
  w34 <- w %>% filter(Tier >= 3)
  low <- t %>% filter(Tier >= 3) %>% group_by(Environment, ERM, Tier) %>%
    summarise(lowest_median = Method[which.min(m)], lowest_q5 = Method[which.min(lo)], .groups = "drop")
  env <- t %>% filter(Method == "MC") %>% select(Environment, ERM, Tier, m) %>%
    pivot_wider(names_from = Environment, values_from = m) %>% mutate(ratio = Freshwater / Marine)
  env_all <- t %>% select(Environment, ERM, Tier, Method, m) %>%
    pivot_wider(names_from = Environment, values_from = m) %>% mutate(log10_ratio = log10(Freshwater / Marine))
  rci <- t %>% filter(!is.na(rci))
  list(
    range_FD = range(t$m[t$ERM == "FD"]), range_TT = range(t$m[t$ERM == "TT"]),
    ssd_vs_mc = with(w, abs(MC / SSD - 1)), ssd_vs_mc_34 = with(w34, abs(MC / SSD - 1)),
    pssd_vs_min = with(w34, PSSD / pmin(SSD, MC)),
    n_pssd_lowest_median = sum(low$lowest_median == "PSSD"), n_pssd_lowest_q5 = sum(low$lowest_q5 == "PSSD"),
    lowest = low, env_ratio_MC = env, env_log10 = range(env_all$log10_ratio),
    rci_FD = range(rci$rci[rci$ERM == "FD"]), rci_TT = range(rci$rci[rci$ERM == "TT"]),
    rci_FW = range(rci$rci[rci$Environment == "Freshwater"]), rci_Mar = range(rci$rci[rci$Environment == "Marine"]),
    rci_pssd = rci %>% filter(Method == "PSSD") %>% select(Environment, ERM, Tier, rci),
    widest_rci = rci %>% filter(Tier >= 2) %>% group_by(Environment, ERM, Tier) %>%
      summarise(widest = Method[which.max(rci)], .groups = "drop"))
}
st <- list(pub = statements(table2, "pub"), corr = statements(table2, "corr"))

# ---- SI workbooks: published (Dec 2025) vs current ----
pub_xlsx <- tempfile(fileext = ".xlsx")
system2("git", c("show", "7163736:output/data/supplemental_information.xlsx"), stdout = pub_xlsx)
cur_xlsx <- "output/data/supplemental_information.xlsx"
rd <- function(f, s) openxlsx::read.xlsx(f, sheet = s, sep.names = " ")
rsd <- list(pub = rd(pub_xlsx, "Monte_Carlo_RSD"), corr = rd(cur_xlsx, "Monte_Carlo_RSD"))
dbsum <- list(pub = rd(pub_xlsx, "Database summary"), corr = rd(cur_xlsx, "Database summary"))
envsens <- list(pub = rd(pub_xlsx, "Fresh_Mar_Sp_Par_sensitivity"), corr = rd(cur_xlsx, "Fresh_Mar_Sp_Par_sensitivity"))

# ---- leave-one-out (Fig. 8) ----
loo_stats <- function(d) {
  d <- d %>% filter(!is.na(diff))
  top <- d %>% group_by(environment, ERM, tier) %>% slice_max(diff, n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(environment, ERM, tier, authoryear, diff)
  bottom <- d %>% group_by(environment, ERM, tier) %>% slice_min(diff, n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(environment, ERM, tier, authoryear, diff)
  list(n = nrow(d), minimal = mean(d$diff >= -0.5 & d$diff <= 2), more_than_double = mean(d$diff > 1),
       more_than_10fold = mean(d$diff > 10), less_than_half = mean(d$diff < -0.5),
       top = top, bottom = bottom, all = d,
       fw_tt_range = range(d$diff[d$environment == "Freshwater" & d$ERM == "Tissue Translocation"]))
}
loo <- list(pub = loo_stats(rd(pub_xlsx, "Leave_one_out_studies")), corr = loo_stats(rd(cur_xlsx, "Leave_one_out_studies")))

# ---- ingestible size range of the food dilution data (Section 4.3.1: "8-70,273 um") ----
fd <- readRDS("data/output/tier1_2_food.RDS")
ingest_range <- if ("max.size.ingest.um" %in% names(fd)) range(fd$max.size.ingest.um, na.rm = TRUE) else NA

saveRDS(list(table2 = table2, statements = st, rsd = rsd, dbsum = dbsum, envsens = envsens, loo = loo,
             ingest_range = ingest_range, date = Sys.Date()), out)
cat("written:", out, "\n")
