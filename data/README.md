# Data Directory

This directory contains the datasets used for the systematic review and meta-analysis.

## Files Structure

- `raw/`: Original, unprocessed data files
- `processed/`: Cleaned and processed datasets ready for analysis
- `README.md`: This file, describing the data structure and sources

## Expected Data Format

CSV files should include the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| `nombre_estudio` | Study name/identifier | "Smith2020", "Jones2019" |
| `grupo` | Treatment group | "control", "intervention" |
| `formulacion` | Myo-inositol formulation | "2g twice daily", "4g once daily" |
| `parametro_metabolico` | Metabolic parameter | "glucose", "insulin", "HOMA-IR" |
| `n` | Sample size | 25, 30, 45 |
| `media` | Mean value | 95.5, 12.3 |
| `desviacion_estandar` | Standard deviation | 8.2, 3.1 |
| `edad_media` | Mean age | 28.5, 32.1 |
| `obesidad` | Obesity status | "obese", "non-obese" |

## Data Sources

Data collected from systematic literature search of randomized controlled trials studying myo-inositol supplementation in women with PCOS.