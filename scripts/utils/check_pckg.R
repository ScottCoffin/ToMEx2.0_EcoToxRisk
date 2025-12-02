#### CHECK PACKAGE ####
# Function to check and install a specific version of an R package
check_and_install_version <- function(package, version) {
  
  # Attempt to unload if already loaded
  if (paste0("package:", package) %in% search()) {
    message("Detaching loaded package...")
    detach(paste0("package:", package), unload = TRUE, character.only = TRUE, force = TRUE)
  }
  
  # Attempt to unload DLLs
  dll <- paste0(package, ".dll")
  if (dll %in% names(getLoadedDLLs())) {
    message("Unloading DLL...")
    try(unloadDLL(dll), silent = TRUE)
  }
  
  # Check versions
  if (requireNamespace(package, quietly = TRUE)) {
    installed_version <- as.character(packageVersion(package))
    
    if (installed_version == version) {
      message(paste("Package", package, "is already at required version:", version))
    } else {
      message(paste("Installing required version:", version))
      remotes::install_version(package, version = version, upgrade = "never")
    }
  } else {
    message(paste("Package", package, "not installed. Installing version", version))
    remotes::install_version(package, version = version, upgrade = "never")
  }
  
  # Verify
  if (requireNamespace(package, quietly = TRUE)) {
    message(paste("Package", package, "is now at version:",
                  as.character(packageVersion(package))))
  } else {
    stop(paste("Failed to install package", package, "version", version))
  }
}


# Example usage
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Check and install a specific version of ggplot2
#check_and_install_version("ggplot2", "3.3.5")