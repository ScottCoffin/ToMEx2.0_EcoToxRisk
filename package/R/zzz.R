.onLoad <- function(libname, pkgname) {
  # Suppress noisy "replacing previous import" conflict warnings on attach/load.
  # Save old option so it can be restored on unload.
  op <- getOption("PSSDplusplus.warn.conflicts.old")
  if (is.null(op)) {
    options(PSSDplusplus.warn.conflicts.old = getOption("warn.conflicts"))
  }
  options(warn.conflicts = FALSE)
}

.onUnload <- function(libpath) {
  old <- getOption("PSSDplusplus.warn.conflicts.old")
  if (!is.null(old)) {
    options(warn.conflicts = old)
    options(PSSDplusplus.warn.conflicts.old = NULL)
  }
}
