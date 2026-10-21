## data-raw/update_tomex2.R
## Refresh the bundled ToMEx 2.0 data from the ToMEx 2.0 database repository
## (SCCWRP/ToMEx_AquaticOrganisms). Run from the repo root:
##
##   TOMEX_DB <- "path/to/ToMEx_AquaticOrganisms"   # after running
##   # ToMEx2.0_Onboarding/ToMEx2.0_Data_Tidying.R there
##   source("package/data-raw/update_tomex2.R")
##
## then re-run package/data-raw/create_tomex2_mini.R.
##
## Writes data/input/aoc_z_tomex2.RDS and data/input/aoc_setup.RDS (manuscript
## pipeline), package/inst/extdata/aoc_z_tomex2.RDS and package/data/tomex2.rda.
##
## PSSDplusplus 0.4.0 data: ToMEx 2.0 v2.0.1 (2026-06-30, commit 137ceaf) rebuilt
## with the body-length and mouth-size logic of body_length_functions.R and the
## row-level corrections in body_length_measured.csv (GitHub issue #6;
## SCCWRP/ToMEx_AquaticOrganisms#5).

stopifnot(exists("TOMEX_DB"), dir.exists(TOMEX_DB))

aoc_z <- readRDS(file.path(TOMEX_DB, "aoc_z_tomex2.RDS"))
aoc_setup <- readRDS(file.path(TOMEX_DB, "aoc_setup.RDS"))

required <- c("body.length.cm", "body.length.basis", "max.size.ingest.mm",
              "max.size.ingest.um", "max.size.ingest.basis", "jams_domain")
missing <- setdiff(required, names(aoc_z))
if (length(missing) > 0) {
  stop("aoc_z_tomex2.RDS lacks ", paste(missing, collapse = ", "),
       "; rebuild it with the current ToMEx2.0_Data_Tidying.R.", call. = FALSE)
}
# never a default ingestible size
stopifnot(all(is.na(aoc_z$max.size.ingest.mm[aoc_z$max.size.ingest.basis == "unresolved"])))

saveRDS(aoc_z, "data/input/aoc_z_tomex2.RDS")
saveRDS(aoc_setup, "data/input/aoc_setup.RDS")
saveRDS(aoc_z, "package/inst/extdata/aoc_z_tomex2.RDS")
tomex2 <- aoc_z
save(tomex2, file = "package/data/tomex2.rda", compress = "xz")

cat("tomex2:", nrow(tomex2), "rows x", ncol(tomex2), "columns\n")
print(table(tomex2$max.size.ingest.basis, useNA = "ifany"))
cat("Saved: data/input/aoc_z_tomex2.RDS, data/input/aoc_setup.RDS,",
    "package/inst/extdata/aoc_z_tomex2.RDS, package/data/tomex2.rda\n")
