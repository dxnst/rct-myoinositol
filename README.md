# RCT Myo-inositol Systematic Review

This repository contains data analysis workflows for systematic review of randomized controlled trials (RCTs) studying the effects of myo-inositol on metabolic parameters in women with polycystic ovary syndrome (PCOS).

## Project Overview

The project aims to analyze RCT data to perform a systematic review and meta-analysis of myo-inositol supplementation effects on metabolic parameters, comparing outcomes between obese and non-obese women with PCOS.

## Research Focus

- **Intervention**: Myo-inositol supplementation
- **Population**: Women with polycystic ovary syndrome (PCOS)
- **Comparison**: Obese vs. non-obese patients
- **Outcomes**: Metabolic parameters (glucose, insulin, lipid profiles, etc.)
- **Study Design**: Randomized controlled trials (RCTs)

## Repository Structure

```
├── data/           # Raw and processed data files
├── analysis/       # R scripts for data analysis
├── results/        # Generated outputs, plots, and reports
├── docs/          # Documentation and methodology
└── README.md      # This file
```

## Getting Started

This project uses R for statistical analysis and meta-analysis. The analysis workflow includes:

1. Data preprocessing and validation
2. Descriptive statistics by study characteristics
3. Meta-analysis with forest plot generation
4. Subgroup analysis by obesity status
5. Sensitivity and heterogeneity assessment

## Data Format

The expected CSV data format includes the following columns:
- `nombre_estudio`: Study name/identifier
- `grupo`: Treatment group (control/intervention)
- `formulacion`: Myo-inositol formulation details
- `parametro_metabolico`: Metabolic parameter measured
- `n`: Sample size
- `media`: Mean value
- `desviacion_estandar`: Standard deviation
- `edad_media`: Mean age
- `obesidad`: Obesity status (obese/non-obese)

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.