## data-raw/create_tomex2_mini.R
## Run once with `source("package/data-raw/create_tomex2_mini.R")` from the
## repo root to regenerate package/data/tomex2_mini.rda.

load("package/data/tomex2.rda")

# ---- Quality filters mirroring run_pssd_reprex() ----
passed <- tomex2[
  tomex2[["tier_zero_tech_f"]] == "Red Criteria Passed" &
  tomex2[["tier_zero_risk_f"]] == "Red Criteria Passed" &
  tomex2[["risk.13"]] != 0 &
  !is.na(tomex2[["effect.metric"]]) &
  tomex2[["effect.metric"]] != "HONEC" &
  !tomex2[["Group"]] %in% c("Bacterium", "Plant"),
]

# ---- Select 5 freshwater + 5 marine species ----
fw_species <- c(
  "Daphnia magna",
  "Hyalella azteca",
  "Ceriodaphnia dubia",
  "Danio rerio",
  "Oncorhynchus mykiss"
)
mar_species <- c(
  "Paracentrotus lividus",
  "Artemia franciscana",
  "Parvocalanus crassirostris",
  "Mytilus galloprovincialis",
  "Brachionus plicatilis"
)

fw  <- passed[passed[["env_f"]] == "Freshwater" & passed[["Species"]] %in% fw_species, ]
mar <- passed[passed[["env_f"]] == "Marine"     & passed[["Species"]] %in% mar_species, ]

# Supplement with any available species if chosen ones are missing
if (nrow(fw) == 0) {
  fw <- passed[passed[["env_f"]] == "Freshwater", ][
    passed[passed[["env_f"]] == "Freshwater", ][["Species"]] %in%
      head(unique(passed[passed[["env_f"]] == "Freshwater", ][["Species"]]), 5), ]
}
if (nrow(mar) == 0) {
  mar <- passed[passed[["env_f"]] == "Marine", ][
    passed[passed[["env_f"]] == "Marine", ][["Species"]] %in%
      head(unique(passed[passed[["env_f"]] == "Marine", ][["Species"]]), 5), ]
}

tomex2_mini <- rbind(fw, mar)

cat("tomex2_mini dimensions:", nrow(tomex2_mini), "x", ncol(tomex2_mini), "\n")
cat("Environments:", paste(unique(tomex2_mini$env_f), collapse = ", "), "\n")
cat("Species:", paste(unique(tomex2_mini$Species), collapse = ", "), "\n")

save(tomex2_mini, file = "package/data/tomex2_mini.rda", compress = "xz")
cat("Saved: package/data/tomex2_mini.rda\n")
