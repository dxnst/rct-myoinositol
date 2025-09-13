# Analysis Directory

This directory contains R scripts for data analysis and meta-analysis.

## Scripts Overview

- `01_data_preprocessing.R`: Data cleaning and validation
- `02_descriptive_analysis.R`: Descriptive statistics and exploratory data analysis
- `03_meta_analysis.R`: Main meta-analysis functions
- `04_forest_plots.R`: Forest plot generation
- `05_subgroup_analysis.R`: Analysis by obesity status and other subgroups
- `06_sensitivity_analysis.R`: Sensitivity analysis and bias assessment

## Dependencies

The analysis requires the following R packages:
- `meta`: For meta-analysis
- `metafor`: Advanced meta-analysis methods
- `forestplot`: For creating forest plots
- `dplyr`: Data manipulation
- `ggplot2`: Data visualization
- `readr`: Data import/export

## Usage

Scripts should be run in numerical order as each may depend on outputs from previous analyses.