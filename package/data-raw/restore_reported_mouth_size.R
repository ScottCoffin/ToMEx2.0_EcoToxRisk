## data-raw/restore_reported_mouth_size.R
## Run once with `source("package/data-raw/restore_reported_mouth_size.R")`
## from the repo root, then re-run data-raw/create_tomex2_mini.R.
##
## GitHub issue #6: the ToMEx 2.0 data-prep pipeline replaced every maximum
## ingestible size (mouth size) with the body-length prediction (Jams et al.
## 2020), overwriting the sizes reported in ToMEx 1.0 (Koelmans et al. 2020).
## The method is: use a reported mouth size whenever one exists for the
## species / life stage; predict from body length only when none is reported.
## This restores the reported values in the ToMEx 2.0 dataset files:
##   data/input/aoc_z_tomex2.RDS
##   package/inst/extdata/aoc_z_tomex2.RDS
##   package/data/tomex2.rda

aoc_setup <- readRDS("data/input/aoc_setup.RDS") # ToMEx 1.0, tidied
aoc_z <- readRDS("data/input/aoc_z_tomex2.RDS")

reported <- aoc_setup[
  aoc_setup$max.size.ingest.reported.estimated %in% "reported" &
    !is.na(aoc_setup$max.size.ingest.mm),
]
reported <- unique(data.frame(
  Species = as.character(reported$species_f),
  life = tolower(as.character(reported$life_f)),
  max.size.ingest.mm = reported$max.size.ingest.mm
))
stopifnot(!anyDuplicated(reported[c("Species", "life")]))

key_z <- paste(as.character(aoc_z$Species), tolower(as.character(aoc_z$life_f)))
key_rep <- paste(reported$Species, reported$life)
idx <- match(key_z, key_rep)
hit <- !is.na(idx)

changed <- hit & (is.na(aoc_z$max.size.ingest.mm) |
  abs(aoc_z$max.size.ingest.mm - reported$max.size.ingest.mm[idx]) > 1e-12)
cat("Rows restored to a reported mouth size:", sum(changed), "\n")
print(table(droplevels(factor(key_z[changed]))))

aoc_z$max.size.ingest.mm[hit] <- reported$max.size.ingest.mm[idx[hit]]
aoc_z$max.size.ingest.um[hit] <- 1000 * aoc_z$max.size.ingest.mm[hit]

saveRDS(aoc_z, "data/input/aoc_z_tomex2.RDS")
saveRDS(aoc_z, "package/inst/extdata/aoc_z_tomex2.RDS")
tomex2 <- aoc_z
save(tomex2, file = "package/data/tomex2.rda", compress = "xz")
cat("Saved: data/input/aoc_z_tomex2.RDS, package/inst/extdata/aoc_z_tomex2.RDS,",
    "package/data/tomex2.rda\n")
