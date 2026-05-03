# Modelling Online Discussion Supplement

This repository contains code to reproduce the results presented in "Modelling online discussion under circadian rhythms and superspreading dynamics" by Joe Meagher and Nial Friel.

All data is distributed with the companion R package [onlineMessageboardActivity](https://github.com/jpmeagher/onlineMessageboardActivity). The repository does not include any data or fitted model objects — scripts must be run before rendering the figures.

## Repository Structure

- **scripts/** — R scripts for model fitting and analysis (run in sequence)
- **figures/** — Quarto documents that produce all manuscript figures and tables
- **code_modules/** — Shared utility functions (plotting theme, summary helpers)
- **models/**, **evidence/**, **predictions/** — Populated by running the scripts

## Reproducing the Analysis

### 1. Fit models and compute results

Run the R scripts in `scripts/` in the following order. Timings are based on a 2.3 GHz Quad-Core Intel Core i5 MacBook Pro using parallel processing.

| Script | Description | Time |
|--------|-------------|------|
| `model_fit_and_evidence.R` | Fits branching process models M1–M5 and computes model evidence via bridge sampling | ~30 min |
| `assess_predictive_performance.R` | Computes ELPD and CRPS metrics for each model | ~100 min |
| `goodness_of_fit.R` | Simulates discussion size distributions to assess goodness-of-fit | — |
| `appendix_frequency_basis.R` | Compares sinusoidal basis functions (K=1,2,3,4) for the appendix | — |
| `appendix_periodic_poisson_point_process.R` | Fits a periodic Poisson point process to immigrant arrivals | — |

### 2. Generate figures and tables

Once the scripts have been run, render each `.qmd` file in the `figures/` folder. These documents are structured to follow the manuscript and will produce all reported figures and tables.
