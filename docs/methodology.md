# Methodology for Myo-inositol PCOS Systematic Review

## Objective

To systematically review and meta-analyze randomized controlled trials examining the effects of myo-inositol supplementation on metabolic parameters in women with polycystic ovary syndrome (PCOS), with particular attention to differences between obese and non-obese populations.

## Research Questions

1. **Primary**: What is the overall effect of myo-inositol supplementation on metabolic parameters in women with PCOS?
2. **Secondary**: Do the effects of myo-inositol differ between obese and non-obese women with PCOS?

## Study Selection Criteria

### Inclusion Criteria
- Randomized controlled trials (RCTs)
- Participants: Women diagnosed with PCOS
- Intervention: Myo-inositol supplementation (any dose, duration, formulation)
- Control: Placebo, no treatment, or standard care
- Outcomes: Metabolic parameters (glucose, insulin, HOMA-IR, lipid profiles)
- Available quantitative data for meta-analysis

### Exclusion Criteria
- Non-randomized studies
- Case reports, case series
- Studies without control groups
- Insufficient data for effect size calculation
- Studies focusing exclusively on fertility outcomes without metabolic data

## Data Extraction

### Study Characteristics
- Author, year, country
- Study design and duration
- Participant characteristics (age, BMI, PCOS diagnostic criteria)
- Intervention details (dose, frequency, duration, formulation)
- Control group description

### Outcome Measures
- Sample sizes by group
- Mean and standard deviation for continuous outcomes
- Number of events for dichotomous outcomes
- Time points of measurement

### Subgroup Classifications
- **Obesity status**: Based on BMI (≥30 kg/m² for obese, <30 kg/m² for non-obese)
- **Intervention duration**: Short-term (<3 months), long-term (≥3 months)
- **Myo-inositol dose**: Low, medium, high (to be defined based on data)

## Statistical Analysis Plan

### Meta-Analysis Methods
- Random-effects model using DerSimonian-Laird method
- Standardized mean difference (SMD) for continuous outcomes
- Risk ratio (RR) for dichotomous outcomes
- 95% confidence intervals for all estimates

### Heterogeneity Assessment
- I² statistic to quantify heterogeneity
- Cochran's Q test for heterogeneity significance
- Tau² for between-study variance

### Subgroup Analysis
- Pre-specified subgroup analysis by obesity status
- Meta-regression for continuous moderators (age, BMI, duration)
- Sensitivity analysis excluding studies with high risk of bias

### Publication Bias
- Funnel plot visual inspection
- Egger's test for publication bias
- Contour-enhanced funnel plots

## Software and Packages

### R Environment
- R version 4.0 or higher
- Packages: `meta`, `metafor`, `forestplot`, `dplyr`, `ggplot2`

### Reporting Standards
- PRISMA 2020 guidelines for systematic reviews
- Forest plots for effect size visualization
- Summary of findings tables

## Quality Assessment

Risk of bias assessment using Cochrane RoB 2 tool:
- Randomization process
- Deviations from intended interventions
- Missing outcome data
- Measurement of the outcome
- Selection of the reported result

## Data Management

- Data extraction performed in duplicate
- Disagreements resolved by consensus or third reviewer
- All data stored in structured CSV format
- Version control using Git for reproducibility