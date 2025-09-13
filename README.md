# RCT Myo-inositol Systematic Review

A comprehensive workflow for analyzing Randomized Controlled Trial (RCT) data to conduct systematic reviews of myo-inositol effects on metabolic parameters in women with polycystic ovary syndrome (PCOS).

## Quick Start

1. Place your CSV data in the `data/` directory with the required format
2. Run the comprehensive analysis:
   ```r
   source("scripts/forest_plots.R")
   results <- run_comprehensive_analysis("data/your_data.csv")
   ```
3. Check the `output/` directory for forest plots and summary results

## Features

- **Meta-analysis**: Random effects models for pooled effect estimates
- **Forest Plots**: Automated generation with subgroup analysis
- **Subgroup Analysis**: By obesity status (obese vs non-obese)
- **Multiple Parameters**: Supports various metabolic outcomes
- **Export Functions**: CSV summaries and high-quality PNG plots

## Data Format

CSV with headers: `nombre_estudio,grupo,formulacion,parametro_metabolico,n,media,desviacion_estandar,edad_media,obesidad`

## Documentation

See [docs/workflow_guide.md](docs/workflow_guide.md) for detailed usage instructions.

## License

BSD 3-Clause License - see [LICENSE](LICENSE) for details.