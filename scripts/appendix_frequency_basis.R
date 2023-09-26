
# packages ----------------------------------------------------------------

# devtools::install_github("jpmeagher/onlineMessageboardActivity")
library(onlineMessageboardActivity)
library(dplyr)
library(xtable)
library(bridgesampling)
library(rstan)

# stan options ------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# data --------------------------------------------------------------------

train_df <- train_df %>% 
  arrange(discussion) %>% 
  group_by(discussion) %>% 
  mutate(
    within_discussion_parent_id = refactor_branching_structure(
      id = id, parent_id = parent_id, is_immigrant = parent_id == 0, S = length(id)
    )
  ) %>% 
  mutate(within_discussion_id = seq_along(id)) %>% 
  ungroup()

# specify hyper-parameters ------------------------------------------------

a <- 48
period <- 24 # 24 hour period
K <- 4 # max frequency
omega <- structure(2 * pi * (1:K / period), dim = K)


# fit one frequency model -------------------------------------------------

fit_freq_1 <- fit_branching_point_process(
  t = train_df$t, branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a, point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = 1, omega = omega[1, drop=F],
  is_hetero = c(T, T),
  seed = 105
)
evidence_freq_1 <- bridge_sampler(fit_freq_1)


# fit two frequency model -------------------------------------------------

fit_freq_2 <- fit_branching_point_process(
  t = train_df$t, branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a, point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = 2, omega = omega[1:2, drop=F],
  is_hetero = c(T, T),
  seed = 105
)
evidence_freq_2 <- bridge_sampler(fit_freq_2)

# fit three frequency model -----------------------------------------------

fit_freq_3 <- fit_branching_point_process(
  t = train_df$t, branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a, point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = 3, omega = omega[1:3, drop=F],
  is_hetero = c(T, T),
  seed = 105
)
evidence_freq_3 <- bridge_sampler(fit_freq_3)

# fit four frequency model ------------------------------------------------

fit_freq_4 <- fit_branching_point_process(
  t = train_df$t, branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a, point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = K, omega = omega,
  is_hetero = c(T, T),
  seed = 105
)
evidence_freq_4 <- bridge_sampler(fit_freq_4)

# save results ------------------------------------------------------------

## models
saveRDS(fit_freq_1, "models/fit_freq_1.RDS")
saveRDS(fit_freq_2, "models/fit_freq_2.RDS")
saveRDS(fit_freq_3, "models/fit_freq_3.RDS")
saveRDS(fit_freq_3, "models/fit_freq_3.RDS")
saveRDS(fit_freq_4, "models/fit_freq_4.RDS")

## evidence
saveRDS(evidence_freq_1, "evidence/evidence_freq_1.RDS")
saveRDS(evidence_freq_2, "evidence/evidence_freq_2.RDS")
saveRDS(evidence_freq_3, "evidence/evidence_freq_3.RDS")
saveRDS(evidence_freq_3, "evidence/evidence_freq_3.RDS")
saveRDS(evidence_freq_4, "evidence/evidence_freq_4.RDS")

