# RCT Myo-inositol Project Configuration

# Project paths
data_path <- "data/"
raw_data_path <- "data/raw/"
processed_data_path <- "data/processed/"
analysis_path <- "analysis/"
results_path <- "results/"
figures_path <- "results/figures/"
tables_path <- "results/tables/"
reports_path <- "results/reports/"
docs_path <- "docs/"

# Analysis parameters
alpha_level <- 0.05
confidence_level <- 0.95
random_effects_method <- "DL"  # DerSimonian-Laird
heterogeneity_threshold <- 0.4  # I² threshold for substantial heterogeneity

# Data validation parameters
required_columns <- c(
  "nombre_estudio",
  "grupo", 
  "formulacion",
  "parametro_metabolico",
  "n",
  "media",
  "desviacion_estandar",
  "edad_media",
  "obesidad"
)

# Study groups
control_groups <- c("control", "placebo", "standard_care")
intervention_groups <- c("intervention", "treatment", "myo-inositol")

# Obesity classification
obesity_cutoff_bmi <- 30  # BMI >= 30 for obese classification

# Output file settings
figure_width <- 12
figure_height <- 8
figure_dpi <- 300
figure_format <- "png"

# Load this configuration in analysis scripts with:
# source("config.R")