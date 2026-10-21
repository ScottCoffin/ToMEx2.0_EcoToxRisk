#### Body length and maximum ingestible size (mouth size) ####
# Shared by:
#   SCCWRP/ToMEx_AquaticOrganisms   ToMEx2.0_Onboarding/ToMEx2.0_Data_Tidying.R
#   ToMEx2.0_EcoToxRisk             scripts/utils/ToMEx2.0_Data_Tidying_functions.R
# Keep the two copies identical.
#
# Body length for each test row is resolved in this order (first available):
#   1. Row-level correction (body_length_measured.csv), keyed by DOI + species,
#      optionally life stage and exposure duration: a length measured in the
#      study, a stage-matched literature reference or a flagged imputation.
#   2. The database record's own length (ToMEx 1.0 row or ToMEx 2.0 template).
#   3. A species x life-stage reference (ToMEx 1.0 summary, or the adult
#      entries in gape_size.csv), joined by the tidying script.
#   4. Nothing: the length is "unresolved". It is never set to a default size.
# Maximum ingestible size is then, in order: a reported mouth size (the
# record's own, else one reported for the same species / life stage), an
# override from gape_size.csv, the Jams et al. (2020) prediction from body
# length, or NA (unresolved; excluded from food dilution downstream).
#
# Jams et al. (2020) model: log10(mouth, mm) = 0.9341 log10(body length, mm) - 1.12.
# Calibrated on animals >= 9 mm; less robust (tends to over-predict) < 30 mm.

library(dplyr)

# body.size.source values that are species-level (adult) size references
species_reference_pattern <- "sealifebase|fishbase|worms|marinespecies|invasions\\.si\\.edu"

#' Jams et al. (2020) calibration domain for a body length (cm)
jams_domain <- function(body.length.cm) {
  mm <- body.length.cm * 10
  dplyr::case_when(
    is.na(mm) ~ NA_character_,
    mm >= 30 ~ "within",
    mm >= 9 ~ "less_robust",
    TRUE ~ "outside_calibration"
  )
}

#' Jams et al. (2020) predicted maximum ingestible size (mm) from body length (cm)
predict_max_size_ingest_mm <- function(body.length.cm,
                                       beta_log10_body_length = 0.9341,
                                       body_length_intercept = 1.12) {
  ifelse(!is.na(body.length.cm) & body.length.cm > 0,
         10^(beta_log10_body_length * log10(body.length.cm * 10) - body_length_intercept),
         NA_real_)
}

normalize_doi <- function(x) {
  x <- tolower(trimws(as.character(x)))
  sub("^(https?://)?(dx\\.)?doi\\.org/", "", x)
}

#' Read the row-level body-length corrections (tier 1)
read_body_length_measured <- function(path) {
  readr::read_csv(path, show_col_types = FALSE,
                  col_types = readr::cols(.default = readr::col_character(),
                                          exposure.duration.d = readr::col_double(),
                                          body.length.cm = readr::col_double(),
                                          length_min_cm = readr::col_double(),
                                          length_max_cm = readr::col_double()))
}

#' Apply row-level body-length corrections.
#'
#' `measured` has columns doi (alternatives separated by "|"), species_f,
#' life_f (blank = any), exposure.duration.d (blank = any), body.length.cm,
#' body.length.basis, evidence.tier and body.size.source. Matching rows get
#' that length, basis, tier and source. `body.length.cm.database` keeps the
#' value the database had.
apply_body_length_measured <- function(df, measured) {
  df$body.length.cm.database <- df$body.length.cm
  if (!"body.length.basis" %in% names(df)) df$body.length.basis <- NA_character_
  if (!"body.length.tier" %in% names(df)) df$body.length.tier <- NA_character_
  doi <- normalize_doi(df$doi)
  species <- as.character(df$species_f)
  life <- tolower(as.character(df$life_f))
  for (k in seq_len(nrow(measured))) {
    m <- measured[k, ]
    hit <- doi %in% normalize_doi(strsplit(m$doi, "|", fixed = TRUE)[[1]]) &
      species == m$species_f
    if (!is.na(m$life_f) && nzchar(m$life_f)) hit <- hit & life == tolower(m$life_f)
    if (!is.na(m$exposure.duration.d)) hit <- hit & df$exposure.duration.d %in% m$exposure.duration.d
    hit[is.na(hit)] <- FALSE
    if (!any(hit)) {
      warning("body_length_measured.csv entry matched no rows: ", m$species_f, " / ", m$doi,
              call. = FALSE)
      next
    }
    df$body.length.cm[hit] <- m$body.length.cm
    df$body.length.basis[hit] <- m$body.length.basis
    df$body.length.tier[hit] <- m$evidence.tier
    df$body.size.source[hit] <- m$body.size.source
  }
  df
}

#' Classify where each row's body length came from (rows not set by a
#' correction) and flag non-adult rows that carry an adult size reference.
#'
#' `record_length` is TRUE where the database record itself carried the
#' length (before any species-level fill); `tomex1` is TRUE for ToMEx 1.0 rows.
classify_body_length_basis <- function(df, record_length, tomex1) {
  src <- tolower(iconv(as.character(df$body.size.source), "UTF-8", "ASCII", sub = ""))
  life <- tolower(as.character(df$life_f))
  basis <- dplyr::case_when(
    !is.na(df$body.length.basis) ~ df$body.length.basis,
    is.na(df$body.length.cm) ~ "unresolved",
    tomex1 & src %in% c("reported", "study") ~ "study-reported (ToMEx 1.0)",
    !tomex1 & record_length ~ "ToMEx 2.0 template",
    src == "koelmans 2020" ~ "literature (Koelmans et al. 2020)",
    grepl(species_reference_pattern, src) ~ "adult species reference",
    TRUE ~ "literature reference"
  )
  # adult length of each species, as carried by its adult rows
  adult_len <- tapply(df$body.length.cm[life == "adult"],
                      as.character(df$species_f)[life == "adult"],
                      max, na.rm = TRUE)
  ref <- unname(adult_len[as.character(df$species_f)])
  mismatch <- life %in% c("early", "juvenile") &
    basis %in% c("adult species reference", "literature (Koelmans et al. 2020)",
                 "literature reference") &
    (basis == "adult species reference" | (!is.na(ref) & df$body.length.cm >= ref))
  mismatch[is.na(mismatch)] <- FALSE
  basis[mismatch] <- "adult reference (stage mismatch)"
  df$body.length.basis <- basis
  df$jams_domain <- jams_domain(df$body.length.cm)
  df
}

#' Assign maximum ingestible size (mm and um) and its basis.
#'
#' df needs species_f, life_f, body.length.cm, max.size.ingest.mm and a logical
#' mouth.size.reported (TRUE where max.size.ingest.mm is a reported value).
#' `overrides` (optional) has species_f, life_f and max.size.ingest.override.um.
assign_max_size_ingest <- function(df, beta_log10_body_length = 0.9341,
                                   body_length_intercept = 1.12,
                                   overrides = NULL) {
  reported_gape <- df %>%
    filter(mouth.size.reported) %>%
    group_by(species_f, life_f) %>%
    summarise(max.size.ingest.mm.reported = first(max.size.ingest.mm), .groups = "drop")
  df <- df %>% left_join(reported_gape, by = c("species_f", "life_f"))
  if (!is.null(overrides) && nrow(overrides) > 0) {
    ov <- overrides %>%
      transmute(.sp = as.character(species_f), .life = tolower(as.character(life_f)),
                max.size.ingest.override.um = as.numeric(max.size.ingest.override.um)) %>%
      filter(!is.na(max.size.ingest.override.um)) %>%
      distinct(.sp, .life, .keep_all = TRUE)
    df <- df %>%
      mutate(.sp = as.character(species_f), .life = tolower(as.character(life_f))) %>%
      left_join(ov, by = c(".sp", ".life")) %>%
      select(-.sp, -.life)
  } else {
    df$max.size.ingest.override.um <- NA_real_
  }
  predicted <- predict_max_size_ingest_mm(df$body.length.cm, beta_log10_body_length,
                                          body_length_intercept)
  df %>%
    mutate(
      max.size.ingest.basis = case_when(
        mouth.size.reported ~ "reported (record)",
        !is.na(max.size.ingest.mm.reported) ~ "reported (species, life stage)",
        !is.na(max.size.ingest.override.um) ~ "override",
        !is.na(predicted) ~ "predicted (Jams et al. 2020)",
        TRUE ~ "unresolved"),
      max.size.ingest.mm = case_when(
        mouth.size.reported ~ max.size.ingest.mm,
        !is.na(max.size.ingest.mm.reported) ~ max.size.ingest.mm.reported,
        !is.na(max.size.ingest.override.um) ~ max.size.ingest.override.um / 1000,
        TRUE ~ predicted),
      max.size.ingest.reported.estimated = case_when(
        max.size.ingest.basis %in% c("reported (record)", "reported (species, life stage)") ~ "reported",
        max.size.ingest.basis == "override" ~ "override",
        max.size.ingest.basis == "predicted (Jams et al. 2020)" ~ "estimated"),
      max.size.ingest.um = 1000 * max.size.ingest.mm) %>%
    select(-c(max.size.ingest.mm.reported, mouth.size.reported, max.size.ingest.override.um))
}

#' Report test organisms whose maximum ingestible size could not be resolved.
#'
#' These rows are excluded from food dilution downstream (never set to a
#' default size). Writes `path` (CSV) when given, warns, and with
#' strict = TRUE stops instead.
report_unresolved_body_length <- function(df, path = NULL, strict = FALSE,
                                          not_applicable = c("Algae", "Bacterium", "Plant")) {
  group <- if ("org_f" %in% names(df)) df$org_f else df$Group
  u <- df[is.na(df$max.size.ingest.mm) & !as.character(group) %in% not_applicable, ]
  if (nrow(u) == 0) return(invisible(u))
  summary <- u %>%
    mutate(Group = if ("org_f" %in% names(u)) org_f else Group,
           reason = if ("vivo_f" %in% names(u)) {
             ifelse(as.character(vivo_f) %in% "In Vitro", "in vitro: not applicable", "unresolved")
           } else "unresolved") %>%
    count(doi, species_f, life_f, Group, reason, name = "n_rows")
  n_unresolved <- sum(summary$n_rows[summary$reason == "unresolved"])
  if (!is.null(path)) readr::write_csv(summary, path)
  if (n_unresolved == 0) return(invisible(summary))
  msg <- sprintf(paste0("%d test rows (%d species) have no resolved body length or mouth size; ",
                        "they are excluded from food dilution.%s"),
                 n_unresolved, dplyr::n_distinct(summary$species_f[summary$reason == "unresolved"]),
                 if (!is.null(path)) paste0(" See ", path, ".") else "")
  if (strict) stop(msg, call. = FALSE)
  warning(msg, call. = FALSE)
  invisible(summary)
}

#' Warn about mouth sizes predicted outside the robust range of the Jams model
warn_jams_domain <- function(df) {
  basis <- if ("max.size.ingest.basis" %in% names(df)) df$max.size.ingest.basis else
    rep("predicted (Jams et al. 2020)", nrow(df))
  flagged <- basis == "predicted (Jams et al. 2020)" & !is.na(df$jams_domain) &
    df$jams_domain != "within"
  flagged[is.na(flagged)] <- FALSE
  if (any(flagged)) {
    sp <- if ("species_f" %in% names(df)) df$species_f else df$Species
    warning(sprintf(paste0("%d rows (%d species) use the Jams et al. (2020) body-length model ",
                           "for animals < 30 mm (%d < 9 mm, outside its calibration); ",
                           "predicted ingestible sizes are less robust (see jams_domain)."),
                    sum(flagged), dplyr::n_distinct(sp[flagged]),
                    sum(flagged & df$jams_domain == "outside_calibration")),
            call. = FALSE)
  }
  invisible(flagged)
}
