# RCT Myo-inositol Analysis Workflow

## Overview

This repository contains a comprehensive workflow for analyzing Randomized Controlled Trial (RCT) data to conduct a systematic review of myo-inositol effects on metabolic parameters in obese and non-obese women with polycystic ovary syndrome (PCOS).

## Data Structure

The analysis expects CSV data with the following header structure:
```
nombre_estudio,grupo,formulacion,parametro_metabolico,n,media,desviacion_estandar,edad_media,obesidad
```

### Field Descriptions:
- **nombre_estudio**: Study name/identifier
- **grupo**: Group type (control/intervencion)
- **formulacion**: Treatment formulation (e.g., placebo, mio-inositol_2g, d-chiro-inositol_1g)
- **parametro_metabolico**: Metabolic parameter measured (e.g., IMC, glucosa_ayunas, insulina_basal)
- **n**: Sample size
- **media**: Mean value
- **desviacion_estandar**: Standard deviation
- **edad_media**: Mean age
- **obesidad**: Obesity status (si/no)

## Directory Structure

```
rct-myoinositol/
├── data/
│   └── example_data.csv          # Example dataset
├── scripts/
│   ├── setup_environment.R       # R package setup
│   ├── data_analysis.R          # Main analysis script
│   └── forest_plots.R           # Forest plot generation
├── output/
│   ├── forest_plots/            # Generated forest plots
│   └── meta_analysis_summary.csv # Summary results
├── docs/
│   └── workflow_guide.md        # This documentation
└── README.md
```

## Usage

### Prerequisites

The analysis requires R with the following packages:
- meta
- metafor
- forestplot
- dplyr
- readr
- ggplot2
- grid
- gridExtra

### Running the Analysis

1. **Setup Environment** (first time only):
   ```r
   source("scripts/setup_environment.R")
   ```

2. **Basic Analysis**:
   ```r
   source("scripts/data_analysis.R")
   results <- main()
   ```

3. **Comprehensive Analysis with Forest Plots**:
   ```r
   source("scripts/forest_plots.R")
   summary_table <- run_comprehensive_analysis()
   ```

### Custom Data Analysis

To analyze your own data:

1. Place your CSV file in the `data/` directory
2. Ensure it follows the required header structure
3. Run the analysis:
   ```r
   source("scripts/forest_plots.R")
   results <- run_comprehensive_analysis("data/your_data.csv")
   ```

## Output

The analysis generates:

1. **Forest Plots**:
   - Standard forest plots for each metabolic parameter
   - Subgroup analysis by obesity status
   - Saved as PNG files in `output/forest_plots/`

2. **Summary Table**:
   - Effect sizes and confidence intervals
   - Statistical significance tests
   - Heterogeneity measures (I² and τ²)
   - Saved as `output/meta_analysis_summary.csv`

3. **Individual Parameter Data**:
   - Meta-analysis data for each parameter
   - Saved as separate CSV files in `output/`

## Metabolic Parameters Analyzed

The workflow can analyze various metabolic parameters including:
- BMI (IMC)
- Fasting glucose (glucosa_ayunas)
- Basal insulin (insulina_basal)
- HOMA-IR
- Free testosterone (testosterona_libre)
- And any other parameters present in the data

## Statistical Methods

- **Effect Size**: Mean Difference (MD)
- **Meta-analysis Model**: Random Effects (REML)
- **Heterogeneity Assessment**: I² statistic and τ²
- **Subgroup Analysis**: By obesity status (obese vs non-obese)

## Interpretation

- **Negative effect sizes**: Favor intervention (reduction in parameter)
- **Positive effect sizes**: Favor control (increase in parameter)
- **Confidence intervals**: 95% CI provided for all estimates
- **Statistical significance**: P < 0.05 considered significant

## Example Results

The example dataset includes 5 studies examining various metabolic parameters with different inositol formulations in both obese and non-obese women with PCOS.

## Support

For questions or issues with the analysis workflow, please refer to the R documentation for the meta-analysis packages or consult the Cochrane Handbook for Systematic Reviews of Interventions.