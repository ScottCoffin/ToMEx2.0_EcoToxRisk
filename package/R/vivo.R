#' Normalize in vitro / in vivo coding
#'
#' Maps any capitalization, spacing or punctuation of "in vivo"/"invivo" and
#' "in vitro"/"invitro" to `"In Vivo"` or `"In Vitro"`. Any other non-missing
#' value is an error. Missing values are an error with `strict = TRUE`;
#' otherwise they stay `NA` (not eligible for SSDs) with a warning.
#'
#' @param x Character vector of in vitro / in vivo codes.
#' @param doi Optional DOIs (same length as `x`), listed in the warning for
#'   missing values.
#' @param strict Logical; error on missing values.
#' @return Character vector of `"In Vivo"`, `"In Vitro"` or `NA`.
#' @examples
#' normalize_vivo(c("In Vivo", "in vitro", "invivo"))
#' @export
normalize_vivo <- function(x, doi = NULL, strict = FALSE) {
  key <- gsub("[^a-z]", "", tolower(as.character(x)))
  out <- rep(NA_character_, length(key))
  out[key == "invivo"] <- "In Vivo"
  out[key == "invitro"] <- "In Vitro"
  unknown <- !is.na(key) & nzchar(key) & is.na(out)
  if (any(unknown)) {
    stop(
      "Unrecognized in vitro/in vivo value(s): ",
      paste(unique(as.character(x[unknown])), collapse = ", "),
      call. = FALSE
    )
  }
  missing <- is.na(out)
  if (any(missing)) {
    where <- if (!is.null(doi)) {
      paste0(" (DOIs: ", paste(unique(doi[missing]), collapse = ", "), ")")
    } else {
      ""
    }
    msg <- sprintf(
      "%d record(s) have no in vitro/in vivo status%s; they are not eligible for SSDs.",
      sum(missing),
      where
    )
    if (strict) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
  }
  out
}

#' Keep only records eligible for species sensitivity distributions
#'
#' SSDs use in vivo data only: in vitro records (e.g. cell cultures) have no
#' organism-level exposure, body length or mouth size. Every code path that
#' builds SSD input in this package calls this function
#' (the internal `process_environment_data()`, [make_tiered_SSDs()], [run_pssd_reprex()]).
#'
#' @param df Data frame with a `vivo_f` column (`"In Vivo"` / `"In Vitro"`).
#' @param strict Logical; error if any record has no in vitro/in vivo status
#'   (otherwise such records are dropped).
#' @return `df` restricted to rows with `vivo_f == "In Vivo"`.
#' @examples
#' ssd_eligible(data.frame(vivo_f = c("In Vivo", "In Vitro", NA), x = 1:3))
#' @export
ssd_eligible <- function(df, strict = FALSE) {
  if (!"vivo_f" %in% names(df)) {
    stop(
      "ssd_eligible(): data have no vivo_f column; cannot exclude in vitro records.",
      call. = FALSE
    )
  }
  vivo <- as.character(df$vivo_f)
  if (strict && any(is.na(vivo))) {
    stop(
      "ssd_eligible(): ", sum(is.na(vivo)),
      " record(s) have no in vitro/in vivo status.",
      call. = FALSE
    )
  }
  df[!is.na(vivo) & vivo == "In Vivo", , drop = FALSE]
}
