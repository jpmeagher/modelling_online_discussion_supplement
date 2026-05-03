# Packages required -------------------------------------------------------

# devtools::install_github("jpmeagher/onlineMessageboardActivity")
# uncomment line above if onlineMessageboardActivity is not yet installed
library(onlineMessageboardActivity)
library(dplyr)
library(xtable)
library(bridgesampling)
library(rstan)
library(tictoc)

# Stan options ------------------------------------------------------------

options(mc.cores = 4)

# Prepare data ------------------------------------------------------------

# Note that it takes a couple of mintes to prepare the data
# the refactor_branching_structure function is the bottleneck
train_df <- onlineMessageboardActivity::train_df |>
  arrange(discussion) |>
  group_by(discussion) |>
  mutate(
    within_discussion_parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent_id,
      is_immigrant = parent_id == 0,
      S = length(id)
    )
  ) |>
  mutate(
    within_discussion_id = seq_along(id)
  ) |>
  ungroup()
a <- 48

# Specify hyperparameters -------------------------------------------------

period <- 24 # 24 hour period
K <- 2 # Number of sinusoidal basis
omega <- structure(2 * pi * (1:K / period), dim = K) # basis frequencies

# Fit and assess M_1 ------------------------------------------------------

tic("M1 fit")
fit_M1 <- fit_branching_point_process(
  t = train_df$t,
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a,
  seed = 101
)
toc(log = TRUE)

tic("M1 evidence")
evidence_M1 <- bridge_sampler(fit_M1)
toc(log = TRUE)

# Fit and assess M_2 ------------------------------------------------------

tic("M2 fit")
fit_M2 <- fit_branching_point_process(
  t = train_df$t,
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a,
  point_type = (train_df$within_discussion_parent_id > 0) + 1,
  seed = 102
)
toc(log = TRUE)

tic("M2 evidence")
evidence_M2 <- bridge_sampler(fit_M2)
toc(log = TRUE)

# Fit and assess M_3 ------------------------------------------------------

tic("M3 fit")
fit_M3 <- fit_branching_point_process(
  t = train_df$t,
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a,
  point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = K,
  omega = omega,
  seed = 103
)
toc(log = TRUE)

tic("M3 evidence")
evidence_M3 <- bridge_sampler(fit_M3)
toc(log = TRUE)

# Fit and assess M_4 ------------------------------------------------------

tic("M4 fit")
fit_M4 <- fit_branching_point_process(
  t = train_df$t,
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a,
  point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = K,
  omega = omega,
  is_hetero = c(TRUE, TRUE),
  seed = 104
)
toc(log = TRUE)

tic("M4 evidence")
evidence_M4 <- bridge_sampler(fit_M4)
toc(log = TRUE)

# Fit and assess M_5 ------------------------------------------------------

tic("M5 fit")
fit_M5 <- fit_branching_point_process(
  t = train_df$t,
  branching_structure = train_df$within_discussion_parent_id,
  observation_interval = a,
  point_type = (train_df$within_discussion_parent_id > 0) + 1,
  K = K,
  omega = omega,
  is_hetero = c(TRUE, FALSE),
  seed = 105
)
toc(log = TRUE)

tic("M5 evidence")
evidence_M5 <- bridge_sampler(fit_M5)
toc(log = TRUE)

# Collect timing results --------------------------------------------------

timing_log <- tic.log(format = FALSE)

timing_df <- tibble(
  task = purrr::map_chr(timing_log, "msg"),
  elapsed_sec = purrr::map_dbl(timing_log, ~ .x$toc - .x$tic)
) |>
  tidyr::separate_wider_delim(
    task,
    delim = " ",
    names = c("model", "stage")
  ) |>
  tidyr::pivot_wider(
    names_from = stage,
    values_from = elapsed_sec,
    names_prefix = "elapsed_sec_"
  ) |>
  mutate(
    elapsed_min_fit = elapsed_sec_fit / 60,
    elapsed_min_evidence = elapsed_sec_evidence / 60
  )

tic.clearlog()

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

# Timing
saveRDS(timing_df, "models/timing_df.RDS")
