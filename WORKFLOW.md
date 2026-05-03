# Analysis Workflow

This document describes the complete workflow for reproducing the analysis in "Modelling online discussion under circadian rhythms and superspreading dynamics."

## Prerequisites

### Software Requirements
- R (version 4.0+)
- RStudio or similar IDE
- C++ compiler (Rtools on Windows, Xcode on macOS)

### Package Installation

```r
# Install companion package from GitHub
devtools::install_github("jpmeagher/onlineMessageboardActivity")

# Install CRAN dependencies
install.packages(c(
  "rstan", "bridgesampling",
  "dplyr", "magrittr", "lubridate", "zoo",
  "ggplot2", "ggpubr", "latex2exp",
  "HDInterval", "verification",
  "xtable", "reshape2"
))
```

### Stan Configuration

```r
# Set Stan options for optimal performance
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

## Full Analysis Pipeline

### Step 1: Fit Models and Estimate Evidence
**Script:** `scripts/model_fit_and_evidence.R`  
**Duration:** ~30 minutes  
**Output:** 
- `models/fit_M1.RDS` through `models/fit_M5.RDS`
- `evidence/evidence_M1.RDS` through `evidence/evidence_M5.RDS`

**What it does:**
- Fits five candidate Bayesian models (M1-M5) using Stan
- Computes model evidence via bridge sampling
- Models differ in:
  - Homogeneous vs. separate immigrant/offspring parameters
  - Presence of circadian rhythm (sinusoidal component)
  - Heterogeneity in reproduction numbers (dispersion parameters)

**Key models:**
- M1: Simplest (homogeneous, no rhythm, no heterogeneity)
- M4: Full model (separate parameters, rhythm, full heterogeneity)
- M5: Preferred model (separate parameters, rhythm, immigrant heterogeneity only)

### Step 2: Assess Predictive Performance
**Script:** `scripts/assess_predictive_performance.R`  
**Duration:** ~100 minutes  
**Output:** `predictions/prediction_summary.RDS`

**What it does:**
- Evaluates out-of-sample predictive performance on test data
- Computes ELPD (Expected Log Predictive Density)
- Computes CRPS (Continuous Ranked Probability Score)
- Evaluates predictions at multiple observation intervals (0, 1, 2, 4, 8 hours)
- Skill scores computed relative to empirical baseline

**Metrics:**
- Higher ELPD = better predictive performance
- Lower CRPS = better calibration
- Skill score = 1 - (CRPS_model / CRPS_baseline)

### Step 3: Goodness-of-Fit Assessment
**Script:** `scripts/goodness_of_fit.R`  
**Duration:** Variable (depends on simulation parameters)  
**Output:** `predictions/discussion_size.RDS`

**What it does:**
- Simulates discussion sizes from each fitted model
- Compares simulated distributions to observed data
- Performs Kolmogorov-Smirnov tests
- Analyzes time-of-day effects on discussion size

**Analyses:**
- Bootstrap estimates of KS test statistic
- Mean discussion size by hour of day
- Quantile-quantile plots

### Step 4: Appendix - Frequency Basis Selection
**Script:** `scripts/appendix_frequency_basis.R`  
**Duration:** Variable  
**Output:**
- `models/fit_freq_1.RDS` through `models/fit_freq_4.RDS`
- `evidence/evidence_freq_1.RDS` through `evidence/evidence_freq_4.RDS`
- `models/fit_pppm.RDS` (periodic Poisson process model)

**What it does:**
- Compares models with K=1, 2, 3, 4 sinusoidal basis functions
- Determines appropriate number of basis functions for circadian rhythm
- Spectral analysis of hourly activity counts
- Compares immigrant vs. offspring circadian patterns

**Finding:** K=2 provides best balance (evidence strongly favors K>1)

### Step 5: Generate Figures and Tables
**Document:** `figures_and_tables.Rmd`  
**Duration:** ~5 minutes (assuming all data generated)  
**Output:** 
- `figures_and_tables.pdf` - Manuscript supplement
- `figures_and_tables.tex` - LaTeX source
- `figures_and_tables_files/` - Individual figure files

**Contents:**
- **Background section:** Exploratory data analysis (Figures 1-3)
- **Methods section:** Model illustration (Figure 4)
- **Results section:** 
  - Parameter estimates (Table 2)
  - Circadian rhythm inference (Figure 5)
  - Model comparison (Table 3)
  - Predictive performance (Table 4, Figure 6)
  - Goodness-of-fit (Figure 7)
- **Appendices:**
  - Frequency basis selection (Table 5, Figures 8-9)
  - Immigrant arrivals analysis (Figure 10)

## Quick Start: Running Everything

```r
# Set working directory
setwd("path/to/modelling_online_discussion_supplement")

# Run analysis scripts in sequence
source("scripts/model_fit_and_evidence.R")           # ~30 min
source("scripts/assess_predictive_performance.R")    # ~100 min
source("scripts/goodness_of_fit.R")                  # variable
source("scripts/appendix_frequency_basis.R")         # variable

# Generate manuscript outputs
rmarkdown::render("figures_and_tables.Rmd")          # ~5 min
```

**Total time:** ~2-3 hours for complete pipeline

## Partial Workflows

### Only Main Results (Skip Appendix)

```r
source("scripts/model_fit_and_evidence.R")
source("scripts/assess_predictive_performance.R")
source("scripts/goodness_of_fit.R")

# Knit will work but appendix sections will have errors
# Consider commenting out appendix sections in Rmd
```

### Only Model Fitting (No Predictions)

```r
source("scripts/model_fit_and_evidence.R")

# Can still view parameter estimates and model evidence
fit_M4 <- readRDS("models/fit_M4.RDS")
print(fit_M4)
```

### Update Figures Only

If you've already run all scripts and only need to regenerate figures:

```r
# Ensure all .RDS files exist in models/, evidence/, predictions/
rmarkdown::render("figures_and_tables.Rmd")
```

## Data Flow Diagram

```
onlineMessageboardActivity package
  ├── train_df (2017 discussions, 48 hours)
  └── test_df (held-out data)
           ↓
model_fit_and_evidence.R
  ├── models/fit_M1.RDS ... fit_M5.RDS
  └── evidence/evidence_M1.RDS ... evidence_M5.RDS
           ↓
assess_predictive_performance.R
  └── predictions/prediction_summary.RDS
           ↓
goodness_of_fit.R
  └── predictions/discussion_size.RDS
           ↓
appendix_frequency_basis.R
  ├── models/fit_freq_*.RDS
  ├── evidence/evidence_freq_*.RDS
  └── models/fit_pppm.RDS
           ↓
figures_and_tables.Rmd
  ├── figures_and_tables.pdf
  └── figures_and_tables.tex
```

## Verification and Diagnostics

### Check Stan Convergence

```r
fit <- readRDS("models/fit_M4.RDS")
print(fit)

# Look for:
# - Rhat ≈ 1.00 (values > 1.1 indicate poor convergence)
# - n_eff > 100 (effective sample size)
# - No divergent transitions
```

### Check Bridge Sampling Accuracy

```r
evidence <- readRDS("evidence/evidence_M4.RDS")
bridgesampling::error_measures(evidence)

# Coefficient of variation (cv) should be < 0.01
# All models in paper have cv < 0.005
```

### Verify Output Files

```bash
# Check all expected outputs exist
ls models/fit_M*.RDS | wc -l     # Should be 5
ls evidence/evidence_M*.RDS | wc -l  # Should be 5
test -f predictions/prediction_summary.RDS && echo "Predictions exist"
test -f predictions/discussion_size.RDS && echo "GoF data exists"
```

## Common Issues and Solutions

### Issue: Stan compilation fails

**Solution:**
```r
# Reinstall package with verbose output
devtools::install_github("jpmeagher/onlineMessageboardActivity", force = TRUE)

# Check for compiler
pkgbuild::has_compiler()  # Should return TRUE
```

### Issue: Out of memory during Stan sampling

**Solution:**
```r
# Reduce parallel cores
options(mc.cores = 2)

# Or sample fewer iterations (though not recommended)
fit <- fit_branching_point_process(..., iter = 1000, chains = 2)
```

### Issue: Bridge sampling fails or gives unreliable estimates

**Solution:**
```r
# Check Stan fit convergence first
print(fit)

# Increase bridge sampling iterations
evidence <- bridgesampling::bridge_sampler(fit, repetitions = 10)
```

### Issue: Predictions script takes too long

**Solution:**
The script performs many posterior predictive calculations. You can:
- Reduce `n_sims` in the script
- Use fewer posterior samples
- Run overnight or use HPC resources

### Issue: Figures don't render correctly

**Solution:**
```r
# Check all dependencies installed
install.packages(c("ggplot2", "ggpubr", "latex2exp", "HDInterval"))

# Ensure all data files exist
file.exists("models/fit_M4.RDS")
file.exists("predictions/prediction_summary.RDS")

# Try rendering specific chunks interactively
rmarkdown::render("figures_and_tables.Rmd", output_format = "pdf_document")
```

## Customization

### Changing Model Hyperparameters

In `model_fit_and_evidence.R`:

```r
# Modify observation interval (hours)
a <- 48  # Default: 48 hours

# Change circadian period (hours)
period <- 24  # Default: 24 hours

# Adjust number of sinusoidal basis functions
K <- 2  # Default: 2 (justified in appendix)
omega <- 2 * pi * (1:K) / period
```

### Changing Stan Sampling Parameters

```r
fit <- fit_branching_point_process(
  ...,
  iter = 2000,    # Default: 2000 (1000 warmup + 1000 sampling)
  chains = 4,     # Default: 4
  seed = 104      # For reproducibility
)
```

### Modifying Plots

In `figures_and_tables.Rmd`, adjust the `my_theme` object:

```r
my_theme <- theme_classic() +
  theme(
    axis.text = element_text(size = 10),    # Increase font size
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14)
  )
```

## Performance Optimization

### Parallel Processing
Already optimized - uses all available cores by default:
```r
options(mc.cores = parallel::detectCores())
```

### Reduce Computation Time
- Use fewer posterior samples (not recommended for final results)
- Reduce number of simulation replications in goodness-of-fit
- Use faster machine or HPC cluster

### Memory Management
Stan models keep posterior samples in memory. For very large models:
```r
# Extract and save only needed parameters
posterior <- rstan::extract(fit, pars = c("mu", "eta", "psi"))
saveRDS(posterior, "posterior_samples.RDS")
rm(fit)  # Free memory
```

## Reproducibility Notes

- All scripts use `set.seed()` for random number generation
- Stan models use explicit seeds in fitting functions
- Original analysis performed in April 2019 - December 2022
- R version: 4.0+
- Stan version: 2.21+
- Platform: macOS (but should work on Windows/Linux)

## Getting Help

If you encounter issues:
1. Check this workflow document
2. Review `.github/copilot-instructions.md` for technical details
3. Consult `onlineMessageboardActivity` package documentation
4. Review manuscript for methodological questions
5. Open issue on GitHub repository
