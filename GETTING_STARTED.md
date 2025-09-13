# Quick Start Guide

## Prerequisites

- R version 4.0 or higher
- RStudio (recommended but not required)

## Getting Started

1. **Install R dependencies**:
   ```r
   source("install_packages.R")
   ```

2. **Load project configuration**:
   ```r
   source("config.R")
   ```

3. **Prepare your data**:
   - Place raw data files in `data/raw/`
   - Use `data/sample_data_template.csv` as a reference for data format
   - Ensure your CSV files follow the required column structure

4. **Run analysis scripts** (in order):
   ```r
   # Data preprocessing
   source("analysis/01_data_preprocessing.R")
   
   # Descriptive analysis
   source("analysis/02_descriptive_analysis.R")
   
   # Meta-analysis
   source("analysis/03_meta_analysis.R")
   
   # Generate forest plots
   source("analysis/04_forest_plots.R")
   
   # Subgroup analysis
   source("analysis/05_subgroup_analysis.R")
   
   # Sensitivity analysis
   source("analysis/06_sensitivity_analysis.R")
   ```

5. **Check results**:
   - Figures will be saved in `results/figures/`
   - Tables will be saved in `results/tables/`
   - Reports will be saved in `results/reports/`

## Data Format Requirements

Your CSV files must include these columns:
- `nombre_estudio`: Study identifier
- `grupo`: "control" or "intervention"
- `formulacion`: Treatment details
- `parametro_metabolico`: Outcome measure
- `n`: Sample size
- `media`: Mean value
- `desviacion_estandar`: Standard deviation
- `edad_media`: Mean age
- `obesidad`: "obese" or "non-obese"

## Troubleshooting

- Check `docs/methodology.md` for detailed analysis methodology
- Ensure all required packages are installed
- Verify data format matches the template
- Check for missing values in required columns

## Support

For questions about the methodology, see `docs/` directory.
For technical issues, check package documentation and R help files.