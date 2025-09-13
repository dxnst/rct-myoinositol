# Install Required R Packages for RCT Myo-inositol Analysis

# List of required packages
required_packages <- c(
  "meta",          # Meta-analysis functions
  "metafor",       # Advanced meta-analysis methods
  "forestplot",    # Forest plot creation
  "dplyr",         # Data manipulation
  "ggplot2",       # Data visualization
  "readr",         # Reading CSV files
  "readxl",        # Reading Excel files
  "tidyr",         # Data tidying
  "stringr",       # String manipulation
  "knitr",         # Report generation
  "rmarkdown"      # R Markdown documents
)

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  if(length(new_packages) > 0) {
    cat("Installing packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, dependencies = TRUE)
  } else {
    cat("All required packages are already installed.\n")
  }
}

# Install missing packages
install_if_missing(required_packages)

# Load packages to verify installation
cat("\nLoading packages to verify installation:\n")
for(pkg in required_packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat("✓", pkg, "loaded successfully\n")
  }, error = function(e) {
    cat("✗", pkg, "failed to load:", e$message, "\n")
  })
}

cat("\nPackage installation complete!\n")
cat("You can now run the analysis scripts in the analysis/ directory.\n")