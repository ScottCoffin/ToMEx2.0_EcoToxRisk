#### In vitro / in vivo status and SSD eligibility ####
# Shared by:
#   SCCWRP/ToMEx_AquaticOrganisms   (RDAmaker.R, ToMEx2.0_Onboarding/ToMEx2.0_Data_Tidying.R, app.R)
#   ToMEx2.0_EcoToxRisk             (scripts/utils/vivo_functions.R)
# Keep the two copies identical.
#
# Species sensitivity distributions (SSDs) use in vivo data only. In vitro
# records (e.g., cell cultures) have no organism-level exposure, body length or
# mouth size and must never reach an SSD.

library(dplyr)

#' Normalize in vitro / in vivo coding to "In Vivo" or "In Vitro".
#'
#' Accepts any capitalization, spacing or punctuation of "in vivo"/"invivo"
#' and "in vitro"/"invitro". Any other non-missing value is an error. Missing
#' values are an error with strict = TRUE; otherwise they stay NA (not
#' eligible for SSDs) with a warning listing the DOIs.
normalize_vivo <- function(x, doi = NULL, strict = FALSE) {
  key <- gsub("[^a-z]", "", tolower(as.character(x)))
  out <- rep(NA_character_, length(key))
  out[key == "invivo"] <- "In Vivo"
  out[key == "invitro"] <- "In Vitro"
  unknown <- !is.na(key) & nzchar(key) & is.na(out)
  if (any(unknown)) {
    stop("Unrecognized in vitro/in vivo value(s): ",
         paste(unique(as.character(x[unknown])), collapse = ", "), call. = FALSE)
  }
  missing <- is.na(out)
  if (any(missing)) {
    where <- if (!is.null(doi)) paste0(" (DOIs: ", paste(unique(doi[missing]), collapse = ", "), ")") else ""
    msg <- sprintf("%d record(s) have no in vitro/in vivo status%s; they are not eligible for SSDs.",
                   sum(missing), where)
    if (strict) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
  }
  out
}

#' Read record-level corrections of in vitro / in vivo status
read_vivo_corrections <- function(path) {
  readr::read_csv(path, show_col_types = FALSE,
                  col_types = readr::cols(.default = readr::col_character()))
}

#' Apply record-level corrections of in vitro / in vivo status.
#'
#' `corrections` has doi (alternatives separated by "|"), species_f, vivo_f
#' (the corrected value), source and note. Matching rows get the corrected
#' (normalized) vivo_f; `vivo_f.database` keeps the database value.
apply_vivo_corrections <- function(df, corrections) {
  df$vivo_f.database <- as.character(df$vivo_f)
  df$vivo_f <- as.character(df$vivo_f)
  doi <- tolower(trimws(as.character(df$doi)))
  for (k in seq_len(nrow(corrections))) {
    m <- corrections[k, ]
    hit <- doi %in% tolower(trimws(strsplit(m$doi, "|", fixed = TRUE)[[1]])) &
      as.character(df$species_f) == m$species_f
    hit[is.na(hit)] <- FALSE
    if (!any(hit)) {
      warning("vivo_corrections.csv entry matched no rows: ", m$species_f, " / ", m$doi, call. = FALSE)
      next
    }
    df$vivo_f[hit] <- normalize_vivo(m$vivo_f, strict = TRUE)
  }
  df$vivo_f <- factor(df$vivo_f, levels = c("In Vivo", "In Vitro"))
  df
}

#' Reason a record is excluded from SSDs (NA = eligible).
#'
#' "in vitro" and "vivo status missing" exclude a record from every SSD. An
#' in vivo test organism with no resolved maximum ingestible size is excluded
#' from food dilution SSDs only.
ssd_exclusion_reason <- function(df) {
  vivo <- as.character(df$vivo_f)
  group <- if ("org_f" %in% names(df)) as.character(df$org_f) else as.character(df$Group)
  no_mouth <- if ("max.size.ingest.mm" %in% names(df)) is.na(df$max.size.ingest.mm) else
    rep(FALSE, nrow(df))
  dplyr::case_when(
    is.na(vivo) ~ "vivo status missing",
    vivo == "In Vitro" ~ "in vitro",
    no_mouth & !group %in% c("Algae", "Bacterium", "Plant") ~
      "food dilution only: unresolved body length or mouth size",
    TRUE ~ NA_character_
  )
}

#' Keep only records eligible for SSDs (in vivo).
#'
#' Every code path that builds SSD input must call this. Rows that are in
#' vitro or have no in vitro/in vivo status are dropped; with strict = TRUE a
#' missing status is an error.
ssd_eligible <- function(df, strict = FALSE) {
  if (!"vivo_f" %in% names(df)) {
    stop("ssd_eligible(): data have no vivo_f column; cannot exclude in vitro records.", call. = FALSE)
  }
  vivo <- as.character(df$vivo_f)
  if (strict && any(is.na(vivo))) {
    stop("ssd_eligible(): ", sum(is.na(vivo)), " record(s) have no in vitro/in vivo status.", call. = FALSE)
  }
  df[!is.na(vivo) & vivo == "In Vivo", , drop = FALSE]
}

#' List in vitro (and status-missing) records, by DOI and species.
#' Writes `path` (CSV) when given and returns the summary invisibly.
report_excluded_in_vitro <- function(df, path = NULL) {
  vivo <- as.character(df$vivo_f)
  x <- df[is.na(vivo) | vivo != "In Vivo", , drop = FALSE]
  summary <- x %>%
    mutate(reason = ifelse(is.na(as.character(vivo_f)), "vivo status missing", "in vitro"),
           Group = if ("org_f" %in% names(x)) org_f else Group,
           species_f = if ("species_f" %in% names(x)) species_f else Species) %>%
    count(doi, species_f, Group, reason, name = "n_rows") %>%
    mutate(threshold_types = "all SSDs (food dilution, tissue translocation)")
  if (!is.null(path)) readr::write_csv(summary, path)
  invisible(summary)
}
