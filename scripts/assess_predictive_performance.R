# script description ------------------------------------------------------

# Assesses the predictive performance of each candidate model.
# Computes ELPD, CRPS, and stores simulated discussion sizes for
# empirical coverage probability and interval score analysis.
# Expected to run in about 3 hours with 8 parallel workers.

# packages ----------------------------------------------------------------

# devtools::install_github("jpmeagher/onlineMessageboardActivity")
library(onlineMessageboardActivity)
library(dplyr)
library(lubridate)
library(magrittr)
library(future.apply)
library(matrixStats)
library(tictoc)
library(cli)

tic("Total assessment")
cli_h1("Model predictive performance assessment")
cli_alert_info("Starting at {Sys.time()}")

# specify parallel futures ------------------------------------------------

plan("multisession", workers = 8)

# data --------------------------------------------------------------------

# specify the branching structure within each discussion
df <- onlineMessageboardActivity::test_df |>
  arrange(discussion) |>
  group_by(discussion) |>
  mutate(tau = t - min(t)) |>
  mutate(
    within_discussion_parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent_id,
      is_immigrant = parent_id == 0,
      S = length(id)
    )
  ) |>
  mutate(within_discussion_id = seq_along(id)) |>
  ungroup()

# set up data frame to store assessment statistics
prediction_summary <- df |>
  filter(parent_id == 0) |>
  select(discussion, t) |>
  left_join(
    df |>
      group_by(discussion) |>
      summarise(size = length(discussion)) |>
      ungroup(),
    join_by(discussion)
  )

# discussion observation period
a <- 48
cli_alert_success(
  "Data preparation complete ({nrow(df)} observations, {nrow(prediction_summary)} discussions)"
)

# hyperparameters ---------------------------------------------------------

# sinusoidal basis
K <- 2
day <- 24
omega <- 2 * pi * (1:K / day)

# observation intervals
s <- c(0, 2^(0:5))

# number of posterior samples to use for assessment
n_sim <- 100

# random seed
set.seed(123)

# fitted models -----------------------------------------------------------

fit_M1 <- readRDS("models/fit_M1.RDS") |> as.data.frame()
fit_M2 <- readRDS("models/fit_M2.RDS") |> as.data.frame()
fit_M3 <- readRDS("models/fit_M3.RDS") |> as.data.frame()
fit_M4 <- readRDS("models/fit_M4.RDS") |> as.data.frame()
fit_M5 <- readRDS("models/fit_M5.RDS") |> as.data.frame()

# helper: compute CRPS from stored simulated sizes -----------------------

compute_crps_from_sims <- function(sim_sizes, observed_sizes) {
  crps_mat <- matrix(nrow = length(observed_sizes), ncol = length(sim_sizes))
  for (h in seq_along(sim_sizes)) {
    crps_mat[, h] <- sapply(seq_along(observed_sizes), function(i) {
      compute_pwm_crps(
        observed_sizes[i],
        sim_sizes[[h]][, i],
        perform_checks = FALSE
      )
    })
  }
  crps_mat
}

# assess M_1 --------------------------------------------------------------

tmp_M <- fit_M1
n_samples <- nrow(tmp_M)

# assess M_1 elpd ---------------------------------------------------------

cli_h2("Assessing M1")
tic("M1 elpd")
prediction_summary$elpd_M1 <-
  future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |> filter(discussion == i)
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_z <- sapply(1:tmp_N, function(i) sum(tmp_beta == i))
      tmp_a <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        tmp_mu <- rep(tmp_M$`mu[1]`[j], tmp_N)
        tmp_psi <- rep(tmp_M$`psi[1]`[j], tmp_N)
        tmp_eta <- rep(tmp_M$`eta[1]`[j], tmp_N)
        heterogeneous_branching_point_process_marginal_likelihood(
          t = tmp_t,
          beta = tmp_beta,
          a = tmp_a,
          mu_vec = tmp_mu,
          psi_vec = tmp_psi,
          eta_vec = tmp_eta,
          z = tmp_z,
          perform_checks = FALSE
        )
      }) |>
        logSumExp() |>
        subtract(log(n_sim))
    },
    future.seed = TRUE
  )
toc(log = TRUE)
cli_alert_success("M1 ELPD complete")


# assess M_1 simulated sizes and crps ------------------------------------

tic("M1 crps")

sim_sizes_M1 <- vector("list", length(s))
for (h in seq_along(s)) {
  sim_sizes_M1[[h]] <- future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |>
        filter(discussion == i) |>
        filter(tau <= s[h])
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_a <- tmp_t[1] + s[h]
      tmp_b <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        simulate_gpm_cluster_process(
          t_seed = tmp_t,
          branching_structure_seed = tmp_beta,
          observation_horizon = tmp_a,
          simulation_horizon = tmp_b,
          immigrant_reproduction_number = tmp_M$`mu[1]`[j],
          immigrant_dispersion_parameter = tmp_M$`psi[1]`[j],
          immigrant_gi_exp_decay_rate = tmp_M$`eta[1]`[j],
          offspring_reproduction_number = tmp_M$`mu[1]`[j],
          offspring_dispersion_parameter = tmp_M$`psi[1]`[j],
          offspring_gi_exp_decay_rate = tmp_M$`eta[1]`[j],
          refactor_output = FALSE,
          perform_checks = FALSE
        ) |>
          nrow()
      })
    },
    future.seed = TRUE
  )
}
names(sim_sizes_M1) <- paste0("s", s)
crps_M1 <- compute_crps_from_sims(sim_sizes_M1, prediction_summary$size)
prediction_summary <- cbind(prediction_summary, crps_M1 = crps_M1)
rm(crps_M1, tmp_M)
toc(log = TRUE)
cli_alert_success("M1 CRPS complete")

# assess M_2 --------------------------------------------------------------

tmp_M <- fit_M2
n_samples <- nrow(tmp_M)


# assess M_2 elpd ---------------------------------------------------------

cli_h2("Assessing M2")
tic("M2 elpd")
prediction_summary$elpd_M2 <-
  future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |> filter(discussion == i)
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_z <- sapply(1:tmp_N, function(i) sum(tmp_beta == i))
      tmp_a <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        tmp_mu <- c(tmp_M$`mu[1]`[j], rep(tmp_M$`mu[2]`[j], tmp_N - 1))
        tmp_psi <- c(tmp_M$`psi[1]`[j], rep(tmp_M$`psi[2]`[j], tmp_N - 1))
        tmp_eta <- c(tmp_M$`eta[1]`[j], rep(tmp_M$`eta[2]`[j], tmp_N - 1))
        heterogeneous_branching_point_process_marginal_likelihood(
          t = tmp_t,
          beta = tmp_beta,
          a = tmp_a,
          mu_vec = tmp_mu,
          psi_vec = tmp_psi,
          eta_vec = tmp_eta,
          z = tmp_z,
          perform_checks = FALSE
        )
      }) |>
        logSumExp() |>
        subtract(log(n_sim))
    },
    future.seed = TRUE
  )
toc(log = TRUE)
cli_alert_success("M2 ELPD complete")

# assess M_2 simulated sizes and crps ------------------------------------

tic("M2 crps")

sim_sizes_M2 <- vector("list", length(s))
for (h in seq_along(s)) {
  sim_sizes_M2[[h]] <- future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |>
        filter(discussion == i) |>
        filter(tau <= s[h])
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_a <- tmp_t[1] + s[h]
      tmp_b <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        simulate_gpm_cluster_process(
          t_seed = tmp_t,
          branching_structure_seed = tmp_beta,
          observation_horizon = tmp_a,
          simulation_horizon = tmp_b,
          immigrant_reproduction_number = tmp_M$`mu[1]`[j],
          immigrant_dispersion_parameter = tmp_M$`psi[1]`[j],
          immigrant_gi_exp_decay_rate = tmp_M$`eta[1]`[j],
          offspring_reproduction_number = tmp_M$`mu[2]`[j],
          offspring_dispersion_parameter = tmp_M$`psi[2]`[j],
          offspring_gi_exp_decay_rate = tmp_M$`eta[2]`[j],
          refactor_output = FALSE,
          perform_checks = FALSE
        ) |>
          nrow()
      })
    },
    future.seed = TRUE
  )
}
names(sim_sizes_M2) <- paste0("s", s)
crps_M2 <- compute_crps_from_sims(sim_sizes_M2, prediction_summary$size)
prediction_summary <- cbind(prediction_summary, crps_M2 = crps_M2)
rm(crps_M2, tmp_M)
toc(log = TRUE)
cli_alert_success("M2 CRPS complete")

# assess M_3 --------------------------------------------------------------

tmp_M <- fit_M3
n_samples <- nrow(tmp_M)

alpha_ind <- grepl("alpha", colnames(tmp_M))

# assess M_3 elpd ---------------------------------------------------------

cli_h2("Assessing M3")
tic("M3 elpd")
prediction_summary$elpd_M3 <-
  future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |> filter(discussion == i)
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_a <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        tmp_mu <- c(tmp_M$`mu[1]`[j], rep(tmp_M$`mu[2]`[j], tmp_N - 1))
        tmp_psi <- c(tmp_M$`psi[1]`[j], rep(tmp_M$`psi[2]`[j], tmp_N - 1))
        tmp_eta <- c(tmp_M$`eta[1]`[j], rep(tmp_M$`eta[2]`[j], tmp_N - 1))
        tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
        heterogeneous_branching_point_process_marginal_likelihood(
          t = tmp_t,
          beta = tmp_beta,
          a = tmp_a,
          mu_vec = tmp_mu,
          psi_vec = tmp_psi,
          eta_vec = tmp_eta,
          omega = omega,
          alpha = tmp_alpha,
          perform_checks = FALSE
        )
      }) |>
        logSumExp() |>
        subtract(log(n_sim))
    },
    future.seed = TRUE
  )
toc(log = TRUE)
cli_alert_success("M3 ELPD complete")

# assess M_3 simulated sizes and crps ------------------------------------

tic("M3 crps")

sim_sizes_M3 <- vector("list", length(s))
for (h in seq_along(s)) {
  sim_sizes_M3[[h]] <- future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |>
        filter(discussion == i) |>
        filter(tau <= s[h])
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_a <- tmp_t[1] + s[h]
      tmp_b <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
        simulate_gpm_cluster_process(
          t_seed = tmp_t,
          branching_structure_seed = tmp_beta,
          observation_horizon = tmp_a,
          simulation_horizon = tmp_b,
          immigrant_reproduction_number = tmp_M$`mu[1]`[j],
          immigrant_dispersion_parameter = tmp_M$`psi[1]`[j],
          immigrant_gi_exp_decay_rate = tmp_M$`eta[1]`[j],
          offspring_reproduction_number = tmp_M$`mu[2]`[j],
          offspring_dispersion_parameter = tmp_M$`psi[2]`[j],
          offspring_gi_exp_decay_rate = tmp_M$`eta[2]`[j],
          dominating_scalar = 2.5,
          sinusoid_coefficients = tmp_alpha,
          sinusoid_frequencies = omega,
          refactor_output = FALSE,
          perform_checks = FALSE
        ) |>
          nrow()
      })
    },
    future.seed = TRUE
  )
}
names(sim_sizes_M3) <- paste0("s", s)
crps_M3 <- compute_crps_from_sims(sim_sizes_M3, prediction_summary$size)
prediction_summary <- cbind(prediction_summary, crps_M3 = crps_M3)
rm(crps_M3, tmp_M)
toc(log = TRUE)
cli_alert_success("M3 CRPS complete")

# assess M_4 --------------------------------------------------------------

tmp_M <- fit_M4
n_samples <- nrow(tmp_M)

alpha_ind <- grepl("alpha", colnames(tmp_M))

# assess M_4 elpd ---------------------------------------------------------

cli_h2("Assessing M4")
tic("M4 elpd")
prediction_summary$elpd_M4 <-
  future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |> filter(discussion == i)
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_z <- sapply(1:tmp_N, function(i) sum(tmp_beta == i))
      tmp_a <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        tmp_mu <- c(tmp_M$`mu[1]`[j], rep(tmp_M$`mu[2]`[j], tmp_N - 1))
        tmp_psi <- c(tmp_M$`psi[1]`[j], rep(tmp_M$`psi[2]`[j], tmp_N - 1))
        tmp_eta <- c(tmp_M$`eta[1]`[j], rep(tmp_M$`eta[2]`[j], tmp_N - 1))
        tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
        heterogeneous_branching_point_process_marginal_likelihood(
          t = tmp_t,
          beta = tmp_beta,
          a = tmp_a,
          mu_vec = tmp_mu,
          psi_vec = tmp_psi,
          eta_vec = tmp_eta,
          omega = omega,
          alpha = tmp_alpha,
          z = tmp_z,
          perform_checks = FALSE
        )
      }) |>
        logSumExp() |>
        subtract(log(n_sim))
    },
    future.seed = TRUE
  )
toc(log = TRUE)
cli_alert_success("M4 ELPD complete")

# assess M_4 simulated sizes and crps ------------------------------------

tic("M4 crps")

sim_sizes_M4 <- vector("list", length(s))
for (h in seq_along(s)) {
  sim_sizes_M4[[h]] <- future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |>
        filter(discussion == i) |>
        filter(tau <= s[h])
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_a <- tmp_t[1] + s[h]
      tmp_b <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
        simulate_gpm_cluster_process(
          t_seed = tmp_t,
          branching_structure_seed = tmp_beta,
          observation_horizon = tmp_a,
          simulation_horizon = tmp_b,
          immigrant_reproduction_number = tmp_M$`mu[1]`[j],
          immigrant_dispersion_parameter = tmp_M$`psi[1]`[j],
          immigrant_gi_exp_decay_rate = tmp_M$`eta[1]`[j],
          offspring_reproduction_number = tmp_M$`mu[2]`[j],
          offspring_dispersion_parameter = tmp_M$`psi[2]`[j],
          offspring_gi_exp_decay_rate = tmp_M$`eta[2]`[j],
          dominating_scalar = 2.5,
          sinusoid_coefficients = tmp_alpha,
          sinusoid_frequencies = omega,
          refactor_output = FALSE,
          perform_checks = FALSE
        ) |>
          nrow()
      })
    },
    future.seed = TRUE
  )
}
names(sim_sizes_M4) <- paste0("s", s)
crps_M4 <- compute_crps_from_sims(sim_sizes_M4, prediction_summary$size)
prediction_summary <- cbind(prediction_summary, crps_M4 = crps_M4)
rm(crps_M4, tmp_M)
toc(log = TRUE)
cli_alert_success("M4 CRPS complete")

# assess M_5 --------------------------------------------------------------

tmp_M <- fit_M5
n_samples <- nrow(tmp_M)
alpha_ind <- grepl("alpha", colnames(tmp_M))

# assess M_5 elpd ---------------------------------------------------------

cli_h2("Assessing M5")
tic("M5 elpd")
prediction_summary$elpd_M5 <-
  future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |> filter(discussion == i)
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_a <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        tmp_mu <- c(tmp_M$`mu[1]`[j], rep(tmp_M$`mu[2]`[j], tmp_N - 1))
        tmp_psi <- c(tmp_M$`psi[1]`[j], rep(tmp_M$`psi[2]`[j], tmp_N - 1))
        tmp_eta <- c(tmp_M$`eta[1]`[j], rep(tmp_M$`eta[2]`[j], tmp_N - 1))
        tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
        heterogeneous_branching_point_process_marginal_likelihood(
          t = tmp_t,
          beta = tmp_beta,
          a = tmp_a,
          mu_vec = tmp_mu,
          psi_vec = tmp_psi,
          eta_vec = tmp_eta,
          omega = omega,
          alpha = tmp_alpha,
          perform_checks = FALSE
        )
      }) |>
        logSumExp() |>
        subtract(log(n_sim))
    },
    future.seed = TRUE
  )
toc(log = TRUE)
cli_alert_success("M5 ELPD complete")

# assess M_5 simulated sizes and crps ------------------------------------

tic("M5 crps")

sim_sizes_M5 <- vector("list", length(s))
for (h in seq_along(s)) {
  sim_sizes_M5[[h]] <- future_sapply(
    as.numeric(prediction_summary$discussion),
    function(i) {
      tmp_df <- df |>
        filter(discussion == i) |>
        filter(tau <= s[h])
      tmp_N <- nrow(tmp_df)
      tmp_t <- tmp_df$t
      tmp_beta <- tmp_df$within_discussion_parent_id
      tmp_a <- tmp_t[1] + s[h]
      tmp_b <- tmp_t[1] + a
      sapply(seq(1, n_samples, length.out = n_sim), function(j) {
        tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
        simulate_gpm_cluster_process(
          t_seed = tmp_t,
          branching_structure_seed = tmp_beta,
          observation_horizon = tmp_a,
          simulation_horizon = tmp_b,
          immigrant_reproduction_number = tmp_M$`mu[1]`[j],
          immigrant_dispersion_parameter = tmp_M$`psi[1]`[j],
          immigrant_gi_exp_decay_rate = tmp_M$`eta[1]`[j],
          offspring_reproduction_number = tmp_M$`mu[2]`[j],
          offspring_dispersion_parameter = tmp_M$`psi[2]`[j],
          offspring_gi_exp_decay_rate = tmp_M$`eta[2]`[j],
          dominating_scalar = 2.5,
          sinusoid_coefficients = tmp_alpha,
          sinusoid_frequencies = omega,
          refactor_output = FALSE,
          perform_checks = FALSE
        ) |>
          nrow()
      })
    },
    future.seed = TRUE
  )
}
names(sim_sizes_M5) <- paste0("s", s)
crps_M5 <- compute_crps_from_sims(sim_sizes_M5, prediction_summary$size)
prediction_summary <- cbind(prediction_summary, crps_M5 = crps_M5)
rm(crps_M5, tmp_M)
toc(log = TRUE)
cli_alert_success("M5 CRPS complete")

# save model assessment ---------------------------------------------------

cli_h2("Saving results")
saveRDS(prediction_summary, "predictions/prediction_summary.RDS")

# save simulated sizes for coverage and interval score analysis
saveRDS(
  list(
    M1 = sim_sizes_M1,
    M2 = sim_sizes_M2,
    M3 = sim_sizes_M3,
    M4 = sim_sizes_M4,
    M5 = sim_sizes_M5
  ),
  "predictions/simulated_sizes.RDS"
)

toc(log = TRUE)
cli_alert_success("Assessment complete at {Sys.time()}")

# print timing summary
timing_log <- tic.log(format = TRUE)
cli_h3("Timing summary")
purrr::walk(timing_log, cli_alert_info)
tic.clearlog()

# ── Timing summary
# ℹ M1 elpd: 5.107 sec elapsed
# ℹ M1 crps: 988.868 sec elapsed
# ℹ M2 elpd: 4.729 sec elapsed
# ℹ M2 crps: 1036.95 sec elapsed
# ℹ M3 elpd: 21.924 sec elapsed
# ℹ M3 crps: 1115.423 sec elapsed
# ℹ M4 elpd: 19.608 sec elapsed
# ℹ M4 crps: 1174.839 sec elapsed
# ℹ M5 elpd: 21.665 sec elapsed
# ℹ M5 crps: 1139.523 sec elapsed
# ℹ Total assessment: 5885.456 sec elapsed
