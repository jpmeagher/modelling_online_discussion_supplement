
# packages ----------------------------------------------------------------

# devtools::install_github("jpmeagher/onlineMessageboardActivity")
library(onlineMessageboardActivity)
library(dplyr)
# library(xtable)
# library(bridgesampling)
# library(rstan)

# stan options ------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# data --------------------------------------------------------------------

train_df <- train_df %>% 
  filter(parent_id == 0)

# specify hyper-parameters ------------------------------------------------

days <- 21 # number of days in training data
period <- 24 # 24 hour period
K <- 2 # max frequency
omega <- structure(2 * pi * (1:K / period), dim = K)


# fit periodic poisson process model --------------------------------------

fit_pppm <- fit_periodic_point_process(
  t = train_df$t, immigrant_observation_interval = days * period,
  K = K, omega = omega,
  seed = 105
)


# save results ------------------------------------------------------------

## models
saveRDS(fit_pppm, "models/fit_pppm.RDS")
