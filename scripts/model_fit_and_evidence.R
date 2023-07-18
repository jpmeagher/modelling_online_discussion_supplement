# Packages required -------------------------------------------------------

# devtools::install_github("jpmeagher/onlineMessageboardActivity")
# uncomment line above if onlineMessageboardActivity is not yet installed
library(onlineMessageboardActivity)
library(dplyr)
library(xtable)
library(bridgesampling)
library(rstan)

# Stan options ------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Prepare data ------------------------------------------------------------

train_df <- train_df %>% 
  arrange(discussion) %>% 
  group_by(discussion) %>% # group discussions 
  mutate(
    within_discussion_parent_id = refactor_branching_structure(
      id = id, 
      parent_id = parent_id, 
      is_immigrant = parent_id == 0, 
      S = length(id)
    )
  ) %>% # create discussion branching structures
  mutate(
    within_discussion_id = seq_along(id)
    ) %>% 
  ungroup()
a <- 48

# Specify hyperparameters -------------------------------------------------

period <- 24 # 24 hour period
K <- 2 # Number of sinusoidal basis
omega <- structure(2 * pi * (1:K / period), dim = K) # basis frequencies

# Fit and assess M_1 ------------------------------------------------------

fit_M1 <- fit_branching_point_process(
  t = train_df$t, 
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a,
  seed = 101
) # approx 2 min
evidence_M1 <- bridge_sampler(fit_M1) # < 1 min

# Fit and assess M_2 ------------------------------------------------------

fit_M2 <- fit_branching_point_process(
  t = train_df$t, 
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a, 
  point_type = (train_df$within_discussion_parent_id > 0) + 1,
  seed = 102
) # approx 2 min
evidence_M2 <- bridge_sampler(fit_M2) # < 1 min

# Fit and assess M_3 ------------------------------------------------------

fit_M3 <- fit_branching_point_process(
  t = train_df$t, 
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a, 
  point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = K, omega = omega,
  seed = 103
) # approx 10 min
evidence_M3 <- bridge_sampler(fit_M3) # < 1 min

# Fit and assess M_4 ------------------------------------------------------

fit_M4 <- fit_branching_point_process(
  t = train_df$t, 
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a, 
  point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = K, omega = omega,
  is_hetero = c(T, T),
  seed = 104
) # approx 10 minutes
evidence_M4 <- bridge_sampler(fit_M4) # < 1 min

# Fit and assess M_5 ------------------------------------------------------

fit_M5 <- fit_branching_point_process(
  t = train_df$t, 
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a, 
  point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = K, omega = omega,
  is_hetero = c(T, F),
  seed = 105
) # approx 10 minutes
evidence_M5 <- bridge_sampler(fit_M5) # < 1 min

# Save results for analysis -----------------------------------------------

# Models
saveRDS(fit_M1, "models/fit_M1.RDS")
saveRDS(fit_M2, "models/fit_M2.RDS")
saveRDS(fit_M3, "models/fit_M3.RDS")
saveRDS(fit_M4, "models/fit_M4.RDS")
saveRDS(fit_M5, "models/fit_M5.RDS")

# Evidence
saveRDS(evidence_M1, "evidence/evidence_M1.RDS")
saveRDS(evidence_M2, "evidence/evidence_M2.RDS")
saveRDS(evidence_M3, "evidence/evidence_M3.RDS")
saveRDS(evidence_M4, "evidence/evidence_M4.RDS")
saveRDS(evidence_M5, "evidence/evidence_M5.RDS")

