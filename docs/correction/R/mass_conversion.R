# Thresholds in mass units (ug/L) under alternative assumptions for the mean mass of an
# environmental particle in the default size range (1-5000 um), for the correction (E9, Table C12).
#
# A threshold in particles/L (aligned to 1-5000 um) is converted to ug/L by multiplying it with
# the mean particle mass of the environmental particle population, m = E[mass]. m cannot be
# measured over 1-5000 um (detection limits are 10-20 um), so it is extrapolated, and it depends
# strongly on the assumptions. Scenarios:
#   A  published conversion (dimensionally incorrect: the length limits 1 and 5000 used as masses)
#   B  power law of particle mass (alpha_M) between the masses of 1 and 5000 um spheres
#      (method of the Kooi et al. 2021 SI worked example; PSSDplusplus 0.4.1 convert_units_fxn)
#   C  as B, with the alignment's particle shape (ellipsoid: W = R x L, H = H:W x W)
#   D  power law of particle volume (alpha_V) between the volumes of 1 and 5000 um ellipsoids, x density
#   E  number-based power law of particle length (alpha_L; the distribution that rescales the counts),
#      ellipsoids: m = density x E[V(L)]
#   F  as E, particle thickness at most 20 um (fixed-thickness fragments, Chen et al. 2026)
#   G  as E, with the pooled, mesh-corrected surface-water exponent alpha_L = 1.53 (Zhao et al. 2026)
#   H  empirical: mean mass of the particles measured by Kooi et al. (2021, Table S3; detected
#      particles >= 20 um freshwater, >= 11 um marine) x the share of particles above that limit
#      (alpha_L); particles below the detection limit are taken to carry negligible mass
# Parameters: Kooi et al. (2021) Table S4 / Coffin et al. (2026) Table S4 (surface waters), as the
# defaults of PSSDplusplus::align_data().
# Run from the repository root: Rscript docs/correction/R/mass_conversion.R
# Writes docs/correction/results/mass_scenarios.rds.

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
stopifnot(file.exists("ToMEx2.0_EcoToxRisk.Rproj"))

env <- tibble(
  Environment = c("Freshwater", "Marine"),
  aL = c(2.64, 2.07), aV = c(1.68, 1.48), aM = c(1.65, 1.32),
  rho = c(1.04, 1.10), R = c(0.67, 0.77), HW = c(0.67, 0.77),
  kooi_mean_ug = c(8.38, 12.1), kooi_dl_um = c(20, 11))
L1 <- 1; L2 <- 5000

vol <- function(L, R, HW, hmax = Inf) { W <- R * L; pi / 6 * L * W * pmin(HW * W, hmax) } # um3
sph <- function(L) pi / 6 * L^3
mass <- function(v, rho) rho * v * 1e-6                                                    # ug
mux <- function(a, UL, LL) ((1 - a) / (2 - a)) * (UL^(2 - a) - LL^(2 - a)) / (UL^(1 - a) - LL^(1 - a))
# number-based length power law: E[f(L)] over [lo, hi], normalised over [L1, L2]
e_len <- function(f, a, lo = L1, hi = L2) {
  num <- integrate(function(u) f(exp(u)) * exp(u)^(1 - a), log(lo), log(hi), rel.tol = 1e-10)$value
  den <- integrate(function(u) exp(u)^(1 - a), log(L1), log(L2), rel.tol = 1e-10)$value
  num / den
}

m <- env %>% rowwise() %>% mutate(
  A = mux(aM, L2, L1),
  B = mux(aM, mass(sph(L2), rho), mass(sph(L1), rho)),
  C = mux(aM, mass(vol(L2, R, HW), rho), mass(vol(L1, R, HW), rho)),
  D = mass(mux(aV, vol(L2, R, HW), vol(L1, R, HW)), rho),
  E = e_len(function(L) mass(vol(L, R, HW), rho), aL),
  F = e_len(function(L) mass(vol(L, R, HW, hmax = 20), rho), aL),
  G = e_len(function(L) mass(vol(L, R, HW), rho), 1.53),
  H = kooi_mean_ug * e_len(function(L) 1, aL, lo = kooi_dl_um)
) %>% ungroup()

scenarios <- tribble(
  ~scenario, ~label, ~basis,
  "A", "Published conversion (incorrect)", "Length limits 1 and 5000 used as masses (the published analysis)",
  "B", "Mass power law, spheres", "alpha_M between the masses of 1 and 5000 um spheres (Kooi et al. 2021 SI method; v0.4.1 code)",
  "C", "Mass power law, ellipsoids", "As B with the alignment's particle shape (W = R x L, H = H:W x W)",
  "D", "Volume power law x density", "alpha_V between the volumes of 1 and 5000 um ellipsoids, x polymer density",
  "E", "Length power law, ellipsoids", "Number-based length distribution (alpha_L, as used to rescale the counts) with the alignment's shape",
  "F", "Length power law, thickness <= 20 um", "As E with particle thickness capped at 20 um (fixed-thickness fragments; Chen et al. 2026)",
  "G", "Pooled surface-water exponent", "As E with alpha_L = 1.53 (45 surface-water studies, mesh-corrected; Zhao et al. 2026)",
  "H", "Measured mean mass, detected share", "Mean mass of the particles measured by Kooi et al. (2021, Table S3) x share of particles above the detection limit")

mbar <- m %>% select(Environment, A:H) %>% pivot_longer(A:H, names_to = "scenario", values_to = "mean_mass_ug")

# particle thresholds (production run): medians of all methods
ssd <- read.csv("output/data/all_thresholds_combined_quality.csv", check.names = FALSE, row.names = 1) %>%
  filter(label == "risk and tech", Environment %in% c("Freshwater", "Marine")) %>%
  transmute(Environment, Tier, FD = `Food Dilution (Default)`, TT = `Tissue Translocation (Default)`) %>%
  pivot_longer(c(FD, TT), names_to = "ERM", values_to = "particles_L") %>% mutate(Method = "SSD")
mc <- readRDS("scripts/monte carlo/output/summary_stats_base_thresholds_sobol.rds") %>%
  filter(Metric %in% c("Food Dilution (Default)", "Tissue Translocation (Default)")) %>%
  transmute(Environment, Tier, ERM = ifelse(grepl("Food", Metric), "FD", "TT"), particles_L = Median, Method = "MC + SSD")
ps <- readRDS("data/output/PNEC_summary_table_wide.rds") %>%
  transmute(Environment, Tier, ERM = ifelse(grepl("Food", ERM), "FD", "TT"), particles_L = Median, Method = "PSSD++")
thr <- bind_rows(ssd, mc, ps)
mass_thr <- thr %>% inner_join(mbar, by = "Environment", relationship = "many-to-many") %>%
  mutate(ug_L = particles_L * mean_mass_ug)

# published MC + SSD thresholds in ug/L (outputs of the published analysis, git 7163736)
pub_xlsx <- tempfile(fileext = ".xlsx")
system2("git", c("show", "7163736:output/data/supplemental_information.xlsx"), stdout = pub_xlsx)
pub_mass <- openxlsx::read.xlsx(pub_xlsx, sheet = "Monte_Carlo_thresholds_mass", sep.names = " ")

write.csv(mass_thr %>% left_join(scenarios %>% select(scenario, label), by = "scenario") %>%
            select(Method, Environment, ERM, Tier, particles_L, scenario, assumption = label, mean_mass_ug, ug_L),
          "docs/correction/results/mass_thresholds.csv", row.names = FALSE)
saveRDS(list(scenarios = scenarios, parameters = env, mean_mass = mbar, thresholds = mass_thr,
             published_mc_mass = pub_mass, date = Sys.Date()), "docs/correction/results/mass_scenarios.rds")
print(m %>% select(Environment, A:H) %>% mutate(across(A:H, ~ signif(.x, 3))))
cat("written: docs/correction/results/mass_scenarios.rds\n")
