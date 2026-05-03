# GitHub Copilot Instructions for Modelling Online Discussion Supplement

## Project Overview
This repository contains supplementary materials for the manuscript "Modelling online discussion under circadian rhythms and superspreading dynamics" by Joe Meagher and Nial Friel. It provides reproducible analysis using the companion R package [onlineMessageboardActivity](https://github.com/jpmeagher/onlineMessageboardActivity).

## Repository Structure

### Core Analysis Files
- **scripts/** - R scripts for model fitting and analysis (run in sequence)
- **models/** - Fitted Stan models saved as .RDS files
- **evidence/** - Model evidence estimates from bridge sampling
- **predictions/** - Predictive performance metrics and simulations
- **figures_and_tables.Rmd** - Main document generating all manuscript outputs

### Data
All datasets are included in the `onlineMessageboardActivity` package:
- `messageboard_df` - Raw r/ireland subreddit data (April 2019)
- `train_df` - Training data (first 48 hours, 2017 discussions)
- `test_df` - Test data (subsequent observations)

## Model Nomenclature

Five candidate models are compared (M1 through M5):

| Model | Immigrant/Offspring | Circadian Rhythm | Heterogeneity |
|-------|-------------------|------------------|---------------|
| **M1** | Homogeneous | None | None |
| **M2** | Separate μ, η | None | None |
| **M3** | Separate μ, η | Sinusoidal α(t) | None |
| **M4** | Separate μ, η | Sinusoidal α(t) | Both ψ₁, ψ₂ |
| **M5** | Separate μ, η | Sinusoidal α(t) | Immigrant ψ₁ only |

**Key Parameters:**
- **μ** - Expected reproduction number (offspring per point)
- **η** - Generation interval decay rate (exponential distribution)
- **ψ** - Dispersion parameter (negative binomial heterogeneity)
- **α(t)** - Circadian rhythm function (sinusoidal basis with K=2)

## Workflow and Dependencies

### Script Execution Order

Scripts must be run in sequence as later scripts depend on outputs from earlier ones:

1. **model_fit_and_evidence.R** (~30 min)
   - Fits models M1-M5 using Stan
   - Computes model evidence via bridge sampling
   - Outputs: `models/fit_M*.RDS`, `evidence/evidence_M*.RDS`

2. **assess_predictive_performance.R** (~100 min)
   - Computes ELPD and CRPS for each model
   - Evaluates predictions at multiple time horizons
   - Output: `predictions/prediction_summary.RDS`

3. **goodness_of_fit.R** (~time varies)
   - Simulates discussion sizes from fitted models
   - Performs Kolmogorov-Smirnov tests
   - Output: `predictions/discussion_size.RDS`

4. **appendix_frequency_basis.R** (~time varies)
   - Compares models with K=1,2,3,4 basis functions
   - Outputs: `models/fit_freq_*.RDS`, `evidence/evidence_freq_*.RDS`

5. **figures_and_tables.Rmd**
   - Knits manuscript figures and tables
   - Requires all prior scripts completed
   - Output: `figures_and_tables.pdf` and `figures_and_tables.tex`

### Computation Notes
- Original timings based on 2.3 GHz Quad-Core Intel Core i5 MacBook Pro
- Stan models use parallel processing: `options(mc.cores = parallel::detectCores())`
- Bridge sampling typically completes in <1 minute per model

## Code Style & Standards

### R Style
- Follow tidyverse conventions
- Use `|>` pipe operator (not `%>%`)
- Use `dplyr` verbs for data manipulation
- Package functions called with `::` notation in analysis scripts

## Key Statistical Concepts

### Poisson Cluster Process
- **Immigrants**: Points with `parent_id == 0` (new discussion posts)
- **Offspring**: Points with `parent_id > 0` (replies/comments)
- **Cluster**: Tree of points descended from single immigrant
- **Generation interval**: Time between parent and child point

### Circadian Rhythm Modeling
- 24-hour periodicity: `period = 24`
- Sinusoidal basis: `omega = 2 * pi * (1:K) / 24`
- Activity function: `α(t) = 1 + Σ[α_k * sin(ω_k * t) + α_{k+K} * cos(ω_k * t)]`

### Superspreading
- Heterogeneous reproduction numbers via negative binomial
- Dispersion parameter ψ controls variance (smaller = more heterogeneous)
- `ψ → ∞` approaches Poisson (no heterogeneity)

### Model Evaluation Metrics
- **ELPD** (Expected Log Predictive Density) - Higher is better
- **CRPS** (Continuous Ranked Probability Score) - Lower is better
- **Skill Score** - Relative to empirical baseline: `1 - (CRPS_model / CRPS_baseline)`
- **Bayes Factors** - Compare model evidence (log scale)

## Figure and Table Generation

### Plot Settings
Consistent theme applied via `my_theme` (customized `theme_classic()`):
- Font sizes: axis text (8pt), axis titles (8pt), plot title (12pt)
- Colors: viridis discrete palette for model comparisons
- Log scales for generation intervals and discussion sizes

## Package Integration

### Required Package Functions
From `onlineMessageboardActivity`:
- `refactor_branching_structure()` - Convert to within-discussion IDs
- `fit_branching_point_process()` - Fit BPP models
- `fit_periodic_point_process()` - Fit periodic Poisson process
- `simulate_gpm_cluster_process()` - Simulate clusters
- `compute_pwm_crps()` - Compute CRPS metric
- `sinusoidal_function()` - Evaluate sinusoidal basis

### Stan Backend
- Models compiled via `rstantools` during package installation
- Use `bridgesampling::bridge_sampler()` for evidence estimation
- Check convergence: `Rhat ≈ 1`, adequate `n_eff`

## Working with Results

### Loading Fitted Models
```r
fit_M4 <- readRDS("models/fit_M4.RDS")
evidence_M4 <- readRDS("evidence/evidence_M4.RDS")
predictions <- readRDS("predictions/prediction_summary.RDS")
```

### Extracting Posterior Samples
```r
# Get posterior samples as data frame
posterior_df <- fit_M4 |>
  as.data.frame() |>
  select(starts_with("mu"), starts_with("eta"))

# Compute posterior summaries
posterior_df |>
  summarise(across(everything(), list(mean = mean, sd = sd)))
```

### Computing Bayes Factors
```r
bf <- bridgesampling::bayes_factor(evidence_M5, evidence_M4, log = TRUE)
# Positive log BF indicates evidence for M5 over M4
```

## Common Tasks

### Re-fitting a Single Model
If you need to re-fit one model without running the full pipeline:

```r
library(onlineMessageboardActivity)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Prepare data
train_df <- train_df |>
  arrange(discussion) |>
  group_by(discussion) |>
  mutate(
    within_discussion_parent_id = refactor_branching_structure(
      id = id, parent_id = parent_id,
      is_immigrant = parent_id == 0, S = length(id)
    )
  ) |>
  ungroup()

# Fit model (example: M4)
fit_M4 <- fit_branching_point_process(
  t = train_df$t,
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = 48,
  point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = 2,
  omega = 2 * pi * (1:2) / 24,
  is_hetero = c(TRUE, TRUE),
  seed = 104
)

# Compute evidence
evidence_M4 <- bridgesampling::bridge_sampler(fit_M4)

# Save results
saveRDS(fit_M4, "models/fit_M4.RDS")
saveRDS(evidence_M4, "evidence/evidence_M4.RDS")
```

### Regenerating Specific Figures
You can extract and run specific chunks from `figures_and_tables.Rmd` interactively. Most figures depend on fitted models being available in `models/` directory.

## Troubleshooting

### Stan Compilation Issues
- Ensure `onlineMessageboardActivity` package is properly installed
- Try `pkgbuild::compile_dll()` in package directory if needed
- Check R and Rtools versions are compatible

### Memory Issues
- Stan models are memory-intensive with large datasets
- Consider reducing `mc.cores` if encountering memory errors
- Monitor with `options(mc.cores = 2)` if needed

### Bridge Sampling Failures
- Check convergence of Stan fit first (`print(fit)`)
- Ensure adequate posterior samples (default: 4000 post-warmup)
- Check coefficient of variation: `bridgesampling::error_measures(evidence)$cv`

### Missing Dependencies
Key packages required:
- `rstan`, `bridgesampling` - Bayesian inference
- `dplyr`, `magrittr` - Data manipulation
- `ggplot2`, `ggpubr` - Visualization
- `HDInterval`, `verification` - Statistical metrics
- `lubridate`, `zoo` - Time series handling

## Research Context

### Dataset
- **Source**: r/ireland subreddit (Reddit)
- **Period**: April 2019 (3 weeks)
- **Structure**: Discussion threads with parent-child relationships
- **Training**: First 48 hours (2017 discussions)
- **Testing**: Remaining observations

### Scientific Questions
1. Do circadian rhythms affect online discussion activity?
2. Is there evidence of superspreading (heterogeneous reproduction)?
3. Do immigrants and offspring have different offspring processes?
4. Can we predict discussion size from early observations?

### Main Findings
- Strong evidence for 24-hour circadian rhythm
- Decisive support for heterogeneous immigrant reproduction
- Offspring processes relatively homogeneous
- Model M5 provides best predictive performance
