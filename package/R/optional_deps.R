# Graceful fallbacks for optional (Suggests-only) dependencies.
#
# These packages are cosmetic or performance conveniences, not required for
# correct output: crayon (colored console text), tictoc (timing messages),
# doSNOW/foreach (multi-core alignment), future/future.apply (multi-core
# PSSD fitting), and progressr (progress bars). When a package is not
# installed, functions in this file degrade to a plain-text/sequential/no-op
# equivalent instead of raising an error, so core workflows keep working with
# a smaller mandatory dependency footprint.

.pssd_has_ns <- function(pkg) {
  isTRUE(requireNamespace(pkg, quietly = TRUE))
}

.pssd_color <- function(text, color) {
  if (.pssd_has_ns("crayon")) {
    return(getExportedValue("crayon", color)(text))
  }
  text
}

.pssd_yellow <- function(text) .pssd_color(text, "yellow")
.pssd_green <- function(text) .pssd_color(text, "green")
.pssd_red <- function(text) .pssd_color(text, "red")
.pssd_blue <- function(text) .pssd_color(text, "blue")

# tic()/toc() replacements. When tictoc is unavailable, tic() is a no-op and
# toc() returns a list shaped like tictoc's return value (`tic`/`toc` in
# seconds via proc.time()) so `elapsed$toc - elapsed$tic` keeps working.
.pssd_tic_stack <- new.env(parent = emptyenv())
.pssd_tic_stack$times <- numeric(0)

.pssd_tic <- function(msg = NULL) {
  if (.pssd_has_ns("tictoc")) {
    if (is.null(msg)) tictoc::tic() else tictoc::tic(msg)
    return(invisible(NULL))
  }
  .pssd_tic_stack$times <- c(.pssd_tic_stack$times, proc.time()[["elapsed"]])
  invisible(NULL)
}

.pssd_toc <- function() {
  if (.pssd_has_ns("tictoc")) {
    return(tictoc::toc())
  }
  n <- length(.pssd_tic_stack$times)
  start <- if (n > 0) .pssd_tic_stack$times[n] else proc.time()[["elapsed"]]
  if (n > 0) {
    .pssd_tic_stack$times <- .pssd_tic_stack$times[-n]
  }
  list(tic = start, toc = proc.time()[["elapsed"]])
}
