#!/usr/bin/env Rscript

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[[length(hit)]], fixed = TRUE)
}

arg_flag <- function(name) {
  paste0("--", name) %in% args
}

read_text_file <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return("")
  }
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

one_line <- function(x) {
  x <- gsub("\r\n|\r|\n", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

bump_patch_version <- function(version) {
  parts <- strsplit(version, ".", fixed = TRUE)[[1]]
  if (length(parts) != 3 || any(!grepl("^[0-9]+$", parts))) {
    stop("Version must use numeric major.minor.patch form; got: ", version)
  }
  parts[[3]] <- as.character(as.integer(parts[[3]]) + 1L)
  paste(parts, collapse = ".")
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
description_path <- file.path(repo_root, "package", "DESCRIPTION")
if (!file.exists(description_path)) {
  stop("Run this script from the repository root; package/DESCRIPTION was not found.")
}

source_url <- arg_value(
  "source-url",
  Sys.getenv(
    "TOMEX2_SOURCE_URL",
    "https://raw.githubusercontent.com/SCCWRP/ToMEx_AquaticOrganisms/master/aoc_z_tomex2.RDS"
  )
)
source_repo <- arg_value("source-repo", Sys.getenv("TOMEX2_SOURCE_REPO", "SCCWRP/ToMEx_AquaticOrganisms"))
source_path <- arg_value("source-path", Sys.getenv("TOMEX2_SOURCE_PATH", "aoc_z_tomex2.RDS"))
source_commit <- arg_value("source-commit", Sys.getenv("TOMEX2_SOURCE_COMMIT", ""))
source_commit_url <- arg_value("source-commit-url", Sys.getenv("TOMEX2_SOURCE_COMMIT_URL", ""))
source_commit_date <- arg_value("source-commit-date", Sys.getenv("TOMEX2_SOURCE_COMMIT_DATE", ""))
source_message <- arg_value("source-message", Sys.getenv("TOMEX2_SOURCE_MESSAGE", ""))
source_message <- paste(c(source_message, read_text_file(arg_value("source-message-file", ""))), collapse = "\n")
source_message <- trimws(source_message)
source_release_notes <- read_text_file(arg_value("source-release-notes-file", ""))
bump_version <- !arg_flag("no-bump-version") &&
  tolower(Sys.getenv("TOMEX2_BUMP_VERSION", "true")) %in% c("1", "true", "yes")

tmp <- tempfile(fileext = ".RDS")
on.exit(unlink(tmp), add = TRUE)

message("Downloading ", source_url)
utils::download.file(source_url, tmp, mode = "wb", quiet = TRUE)

new_data <- readRDS(tmp)
if (!is.data.frame(new_data)) {
  stop("Downloaded object is not a data frame: ", paste(class(new_data), collapse = ", "))
}
if (nrow(new_data) == 0L || ncol(new_data) == 0L) {
  stop("Downloaded data frame is empty.")
}

root_rds <- file.path(repo_root, "data", "input", "aoc_z_tomex2.RDS")
package_rds <- file.path(repo_root, "package", "inst", "extdata", "aoc_z_tomex2.RDS")
package_rda <- file.path(repo_root, "package", "data", "tomex2.rda")
news_path <- file.path(repo_root, "package", "NEWS.md")

load_tomex2 <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("tomex2", envir = env, inherits = FALSE)) {
    stop(path, " does not contain an object named tomex2.")
  }
  get("tomex2", envir = env, inherits = FALSE)
}

rds_matches_download <- function(path) {
  file.exists(path) && identical(unname(tools::md5sum(path)), unname(tools::md5sum(tmp)))
}

rda_matches_download <- function(path) {
  file.exists(path) && identical(load_tomex2(path), new_data)
}

old_data <- if (file.exists(package_rda)) {
  tryCatch(load_tomex2(package_rda), error = function(e) NULL)
} else {
  NULL
}
old_dim <- if (is.data.frame(old_data)) dim(old_data) else c(NA_integer_, NA_integer_)
new_dim <- dim(new_data)

needs_update <- !rds_matches_download(root_rds) ||
  !rds_matches_download(package_rds) ||
  !rda_matches_download(package_rda)

if (needs_update) {
  dir.create(dirname(root_rds), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(package_rds), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(package_rda), recursive = TRUE, showWarnings = FALSE)

  file.copy(tmp, root_rds, overwrite = TRUE)
  file.copy(tmp, package_rds, overwrite = TRUE)
  tomex2 <- new_data
  save(tomex2, file = package_rda, compress = "xz")
  message("Updated tomex2 data files.")
} else {
  message("No tomex2 data file changes detected.")
}

root_data <- readRDS(root_rds)
package_rds_data <- readRDS(package_rds)
package_data <- load_tomex2(package_rda)
if (!identical(root_data, package_rds_data)) {
  stop("Root RDS and package inst/extdata RDS are not identical as R objects.")
}
if (!identical(root_data, package_data)) {
  stop("package/data/tomex2.rda object is not identical to the source RDS object.")
}

if (needs_update && bump_version) {
  desc <- readLines(description_path, warn = FALSE, encoding = "UTF-8")
  version_line <- grep("^Version:[[:space:]]*", desc)
  if (length(version_line) != 1L) {
    stop("Expected exactly one Version field in package/DESCRIPTION.")
  }
  old_version <- sub("^Version:[[:space:]]*", "", desc[[version_line]])
  new_version <- bump_patch_version(old_version)
  desc[[version_line]] <- paste("Version:", new_version)
  writeLines(desc, description_path, useBytes = TRUE)

  source_ref <- if (nzchar(source_commit)) {
    paste0(source_repo, "@", substr(source_commit, 1L, 7L))
  } else {
    source_repo
  }
  source_link <- if (nzchar(source_commit_url)) source_commit_url else source_url
  source_date <- if (nzchar(source_commit_date)) paste0(" (", source_commit_date, ")") else ""
  commit_message <- one_line(source_message)

  heading <- paste("PSSDplusplus", new_version)
  entry <- c(
    heading,
    paste(rep("=", nchar(heading)), collapse = ""),
    "",
    "Data updates",
    "------------",
    paste0(
      "- Synced bundled `tomex2`, `inst/extdata/aoc_z_tomex2.RDS`, and ",
      "`data/input/aoc_z_tomex2.RDS` with upstream `", source_path, "` from ",
      source_ref, source_date, "."
    ),
    paste0("- Upstream source: ", source_link)
  )
  if (!anyNA(old_dim)) {
    entry <- c(
      entry,
      paste0(
        "- Dataset size changed from ", format(old_dim[[1]], big.mark = ","),
        " rows x ", format(old_dim[[2]], big.mark = ","), " columns to ",
        format(new_dim[[1]], big.mark = ","), " rows x ",
        format(new_dim[[2]], big.mark = ","), " columns."
      )
    )
  }
  if (nzchar(commit_message)) {
    entry <- c(entry, paste0("- Upstream commit message: ", commit_message))
  }
  if (nzchar(source_release_notes)) {
    entry <- c(entry, "- Upstream release notes were included in the GitHub release notes.")
  }
  entry <- c(entry, "")

  old_news <- if (file.exists(news_path)) {
    readLines(news_path, warn = FALSE, encoding = "UTF-8")
  } else {
    character()
  }
  writeLines(c(entry, old_news), news_path, useBytes = TRUE)
  message("Bumped package version from ", old_version, " to ", new_version, ".")
}

message(
  "Validated tomex2 synchronization: ",
  nrow(root_data), " rows x ", ncol(root_data), " columns."
)
