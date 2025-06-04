#### CHECK PACKAGE ####
# Function to check and install a specific version of an R package
check_and_install_version <- function(package, version) {
  # Check if the package is installed
  if (requireNamespace(package, quietly = TRUE)) {
    # Get the installed version
    installed_version <- as.character(packageVersion(package))
    
    # Compare the installed version with the specified version
    if (installed_version == version) {
      message(paste("Package", package, "is already at the required version:", version))
    } else {
      message(paste("Package", package, "is installed but at version", installed_version, 
                    "instead of", version, ". Installing the specified version..."))
      remotes::install_version(package, version = version)
    }
  } else {
    # If the package is not installed, install the specified version
    message(paste("Package", package, "is not installed. Installing version", version, "..."))
    remotes::install_version(package, version = version)
  }
  
  # Verify the installation
  if (requireNamespace(package, quietly = TRUE)) {
    message(paste("Package", package, "is now at version:", as.character(packageVersion(package))))
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