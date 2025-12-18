#' Validate an installed package version
#'
#' Lightweight helper to enforce a specific dependency version without
#' attempting to install anything (CRAN-friendly).
#'
#' @param package Character scalar, package name.
#' @param required_version Character scalar specifying the exact version.
#' @return Invisibly returns `TRUE` when the required version is present;
#'   otherwise throws an informative error.
#' @export
check_package_version <- function(package, required_version) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' is not installed; please install version %s.", package, required_version),
      call. = FALSE
    )
  }

  installed_version <- as.character(utils::packageVersion(package))
  if (installed_version != required_version) {
    stop(
      sprintf(
        "Package '%s' must be at version %s (found %s). Please install the required version.",
        package,
        required_version,
        installed_version
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Ensure `ssdtools` is at the expected version (0.3.7)
#'
#' This package relies on `ssdtools` 0.3.7 to reproduce published PSSD fits.
#'
#' @param required_version Character version string. Defaults to `"0.3.7"`.
#' @export
check_ssdtools_version <- function(required_version = "0.3.7") {
  check_package_version("ssdtools", required_version)
}
