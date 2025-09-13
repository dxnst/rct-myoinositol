# RCT Myo-inositol Analysis Scripts

## setup_environment.R
# Script to install and load required R packages for the analysis

# Function to install and load packages
install_if_missing <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
}

# Required packages for meta-analysis and forest plots
required_packages <- c(
  "meta",       # Meta-analysis package
  "metafor",    # Meta-analysis framework
  "forestplot", # Forest plot creation
  "dplyr",      # Data manipulation
  "readr",      # Reading CSV files
  "ggplot2",    # Plotting
  "grid",       # Grid graphics
  "gridExtra"   # Extra grid functionality
)

# Install and load packages
install_if_missing(required_packages)

# Display session info
cat("R packages loaded successfully!\n")
sessionInfo()