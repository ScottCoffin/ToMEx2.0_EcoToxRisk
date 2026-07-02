#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = "") {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[[length(hit)]], fixed = TRUE)
}

fmt_int <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

clean_doi <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a")] <- NA_character_
  tolower(x)
}

doi_change_lines <- function(previous, current) {
  if (!("doi" %in% names(previous)) || !("doi" %in% names(current))) {
    return("- DOI check: `doi` column unavailable in one or both datasets")
  }

  previous_doi <- clean_doi(previous$doi)
  current_doi <- clean_doi(current$doi)
  previous_doi <- previous_doi[!is.na(previous_doi)]
  current_doi <- current_doi[!is.na(current_doi)]

  new_only <- sort(setdiff(unique(current_doi), unique(previous_doi)))

  previous_counts <- as.data.frame(table(previous_doi), stringsAsFactors = FALSE)
  names(previous_counts) <- c("doi", "previous_rows")
  current_counts <- as.data.frame(table(current_doi), stringsAsFactors = FALSE)
  names(current_counts) <- c("doi", "current_rows")
  counts <- merge(current_counts, previous_counts, by = "doi", all.x = TRUE)
  counts$previous_rows[is.na(counts$previous_rows)] <- 0L
  counts$delta <- counts$current_rows - counts$previous_rows
  counts <- counts[counts$delta > 0L, c("doi", "previous_rows", "current_rows", "delta")]
  counts <- counts[order(counts$doi), , drop = FALSE]

  lines <- c(
    paste0("- New unique DOI values: ", fmt_int(length(new_only))),
    paste0("- DOI values with increased row counts: ", fmt_int(nrow(counts)))
  )

  if (length(new_only) > 0L) {
    lines <- c(lines, "New DOI values:")
    lines <- c(lines, paste0("- ", new_only))
  }

  if (nrow(counts) > 0L) {
    lines <- c(lines, "DOI values with increased rows:")
    lines <- c(
      lines,
      paste0(
        "- ", counts$doi, ": ", fmt_int(counts$previous_rows), " -> ",
        fmt_int(counts$current_rows), " (", ifelse(counts$delta >= 0L, "+", ""),
        fmt_int(counts$delta), ")"
      )
    )
  }

  lines
}

current_path <- arg_value("current", "data/input/aoc_z_tomex2.RDS")
previous_path <- arg_value("previous", "")
output_path <- arg_value("output", "tomex2-size-report.md")

if (!file.exists(current_path)) {
  stop("Current dataset file not found: ", current_path)
}

current <- readRDS(current_path)
if (!is.data.frame(current)) {
  stop("Current dataset is not a data frame: ", paste(class(current), collapse = ", "))
}

current_rows <- nrow(current)
current_cols <- ncol(current)

if (nzchar(previous_path) && file.exists(previous_path)) {
  previous <- readRDS(previous_path)
  if (!is.data.frame(previous)) {
    stop("Previous dataset is not a data frame: ", paste(class(previous), collapse = ", "))
  }

  previous_rows <- nrow(previous)
  previous_cols <- ncol(previous)
  row_delta <- current_rows - previous_rows
  col_delta <- current_cols - previous_cols
  rows_changed <- row_delta != 0L
  cols_changed <- col_delta != 0L

  lines <- c(
    paste0("- Previous size: ", fmt_int(previous_rows), " rows x ", fmt_int(previous_cols), " columns"),
    paste0("- Current size: ", fmt_int(current_rows), " rows x ", fmt_int(current_cols), " columns"),
    paste0(
      "- Change: ", ifelse(row_delta >= 0L, "+", ""), fmt_int(row_delta),
      " rows; ", ifelse(col_delta >= 0L, "+", ""), fmt_int(col_delta), " columns"
    ),
    paste0(
      "- Flag: ",
      if (rows_changed || cols_changed) {
        paste(
          c(
            if (rows_changed) "row count changed" else NULL,
            if (cols_changed) "column count changed" else NULL
          ),
          collapse = "; "
        )
      } else {
        "row and column counts unchanged"
      }
    )
  )

  lines <- c(lines, "", "DOI check:", doi_change_lines(previous, current))
} else {
  previous_rows <- NA_integer_
  previous_cols <- NA_integer_
  row_delta <- NA_integer_
  col_delta <- NA_integer_
  rows_changed <- NA
  cols_changed <- NA
  lines <- c(
    "- Previous size: unavailable",
    paste0("- Current size: ", fmt_int(current_rows), " rows x ", fmt_int(current_cols), " columns"),
    "- Flag: previous dataset size unavailable for comparison"
  )
}

writeLines(lines, output_path, useBytes = TRUE)

github_output <- Sys.getenv("GITHUB_OUTPUT", "")
if (nzchar(github_output)) {
  cat(
    paste0("current_rows=", current_rows, "\n"),
    paste0("current_cols=", current_cols, "\n"),
    paste0("previous_rows=", ifelse(is.na(previous_rows), "", previous_rows), "\n"),
    paste0("previous_cols=", ifelse(is.na(previous_cols), "", previous_cols), "\n"),
    paste0("row_delta=", ifelse(is.na(row_delta), "", row_delta), "\n"),
    paste0("col_delta=", ifelse(is.na(col_delta), "", col_delta), "\n"),
    paste0("rows_changed=", tolower(as.character(rows_changed)), "\n"),
    paste0("cols_changed=", tolower(as.character(cols_changed)), "\n"),
    file = github_output,
    append = TRUE,
    sep = ""
  )
}
