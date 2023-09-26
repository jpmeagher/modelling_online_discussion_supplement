# script description ------------------------------------------------------

# assesses the predictive performance of each candidate model expected to run in
# about 24 hours of compute time -  easy to run in parallel.
# we reduce this to 3 hours by paralellising the analysis
timestamp()
print('running model assessment')

# packages ----------------------------------------------------------------

# devtools::install_github("jpmeagher/onlineMessageboardActivity")
library(onlineMessageboardActivity)
library(dplyr)
library(lubridate)
library(magrittr)
library(future.apply)
library(matrixStats)

# specify parallel futures ------------------------------------------------

plan("multisession", workers = 8)

# data --------------------------------------------------------------------

# specify the branching structure within each discussion
df <- test_df %>%
  arrange(discussion) %>%
  group_by(discussion) %>%
  mutate(tau = t - min(t)) %>% 
  mutate(
    within_discussion_parent_id = refactor_branching_structure(
      id = id,
      parent_id = parent_id,
      is_immigrant = parent_id == 0,
      S = length(id)
    )
  ) %>%
  mutate(within_discussion_id = seq_along(id)) %>%
  ungroup()
# expect to take about 5 minutes with 8 workers

# set up data frame to store assessment statistics 
prediction_summary <- df %>% 
  filter(parent_id == 0) %>% 
  select(discussion, t) %>% 
  left_join(
    df %>% 
      group_by(discussion) %>% 
      summarise(size = length(discussion)) %>% 
      ungroup()
  )

# discussion observation period
a <- 48
timestamp()
print('data preparation complete')

# hyperparameters ---------------------------------------------------------

# sinusoidal basis
K <- 2
day <- 24
omega <- 2 * pi * (1:K / day)

# observation intervals
s <- c(0, 2^(0:5))

# fitted models -----------------------------------------------------------

fit_M1 <- readRDS("models/fit_M1.RDS") %>% as.data.frame()
fit_M2 <- readRDS("models/fit_M2.RDS") %>% as.data.frame()
fit_M3 <- readRDS("models/fit_M3.RDS") %>% as.data.frame()
fit_M4 <- readRDS("models/fit_M4.RDS") %>% as.data.frame()
fit_M5 <- readRDS("models/fit_M5.RDS") %>% as.data.frame()

# assess M_1 --------------------------------------------------------------

tmp_M <- fit_M1
n_samples <- nrow(tmp_M)
n_sim <- 100

# assess M_1 elpd ---------------------------------------------------------

prediction_summary$elpd_M1 <-
  future_sapply(as.numeric(prediction_summary$discussion), function(i) {
    tmp_df <- df %>% 
      filter(discussion == i)
    tmp_N <- nrow(tmp_df)
    tmp_t <- tmp_df$t
    tmp_beta <- tmp_df$within_discussion_parent_id
    tmp_z <- sapply(1:tmp_N, function(i) sum(tmp_beta == i))
    tmp_a <- tmp_t[1] + a
    sapply(seq(1, n_samples, length.out = n_sim),
           function(j) {
             tmp_mu <- rep(tmp_M$`mu[1]`[j], tmp_N)
             tmp_psi <- rep(tmp_M$`psi[1]`[j], tmp_N)
             tmp_eta <- rep(tmp_M$`eta[1]`[j], tmp_N)
             heterogeneous_branching_point_process_marginal_likelihood(
               t = tmp_t, beta = tmp_beta, a = tmp_a,
               mu_vec = tmp_mu, psi_vec = tmp_psi, eta_vec = tmp_eta, 
               z = tmp_z,
               perform_checks = F
             )
           }) %>%
      logSumExp() %>% 
      subtract(log(n_sim))
  }, future.seed=TRUE)
timestamp()
print('M_1 elpd complete')


# assess M_1 crps ---------------------------------------------------------

crps_M1 <- matrix(nrow = nrow(prediction_summary), ncol = length(s))
for (h in seq_along(s)) {
  crps_M1[, h] <- future_sapply(as.numeric(prediction_summary$discussion),
                                function(i) {
                                  tmp_df <- df %>%
                                    filter(discussion == i) %>%
                                    filter(tau <= s[h])
                                  tmp_size <- prediction_summary %>%
                                    filter(discussion == i) %>%
                                    use_series(size)
                                  tmp_N <- nrow(tmp_df)
                                  tmp_t <- tmp_df$t
                                  tmp_beta <- tmp_df$within_discussion_parent_id
                                  tmp_a <- tmp_t[1] + s[h]
                                  tmp_b <- tmp_t[1] + a
                                  sapply(seq(1, n_samples, length.out = n_sim),
                                         function(j) {
                                           tmp_mu_imm <- tmp_M$`mu[1]`[j]
                                           tmp_psi_imm <- tmp_M$`psi[1]`[j]
                                           tmp_eta_imm <- tmp_M$`eta[1]`[j]
                                           tmp_mu_off <- tmp_M$`mu[1]`[j]
                                           tmp_psi_off <- tmp_M$`psi[1]`[j]
                                           tmp_eta_off <- tmp_M$`eta[1]`[j]
                                           simulate_gpm_cluster_process(
                                             t_seed = tmp_t, branching_structure_seed = tmp_beta,
                                             observation_horizon = tmp_a,
                                             simulation_horizon = tmp_b,
                                             immigrant_reproduction_number = tmp_mu_imm,
                                             immigrant_dispersion_parameter = tmp_psi_imm,
                                             immigrant_gi_exp_decay_rate = tmp_eta_imm,
                                             offspring_reproduction_number = tmp_mu_off,
                                             offspring_dispersion_parameter = tmp_psi_off,
                                             offspring_gi_exp_decay_rate = tmp_eta_off,
                                             refactor_output = FALSE,
                                             perform_checks = FALSE
                                           ) %>%
                                             nrow()
                                         }) %>%
                                    compute_pwm_crps(tmp_size, ., perform_checks = FALSE)
                                }, future.seed=TRUE)
}
prediction_summary <- cbind(prediction_summary, crps_M1 = crps_M1)
# housekeeping
rm(crps_M1, tmp_M)
timestamp()
print('M_1 crps complete')

# assess M_2 --------------------------------------------------------------

tmp_M <- fit_M2
n_samples <- nrow(tmp_M) 
n_sim <- 100

# assess M_2 elpd ---------------------------------------------------------

prediction_summary$elpd_M2 <-
  future_sapply(as.numeric(prediction_summary$discussion), function(i) {
    tmp_df <- df %>% 
      filter(discussion == i)
    tmp_N <- nrow(tmp_df)
    tmp_t <- tmp_df$t
    tmp_beta <- tmp_df$within_discussion_parent_id
    tmp_z <- sapply(1:tmp_N, function(i) sum(tmp_beta == i))
    tmp_a <- tmp_t[1] + a
    sapply(seq(1, n_samples, length.out = n_sim),
           function(j) {
             tmp_mu <- c(tmp_M$`mu[1]`[j], rep(tmp_M$`mu[2]`[j], tmp_N-1))
             tmp_psi <- c(tmp_M$`psi[1]`[j], rep(tmp_M$`psi[2]`[j], tmp_N-1))
             tmp_eta <- c(tmp_M$`eta[1]`[j], rep(tmp_M$`eta[2]`[j], tmp_N-1))
             heterogeneous_branching_point_process_marginal_likelihood(
               t = tmp_t, beta = tmp_beta, a = tmp_a,
               mu_vec = tmp_mu, psi_vec = tmp_psi, eta_vec = tmp_eta, 
               z = tmp_z,
               perform_checks = F
             )
           }) %>%
      logSumExp() %>% 
      subtract(log(n_sim))
  }, future.seed=TRUE)
timestamp()
print('M_2 elpd complete')

# assess M_2 crps ---------------------------------------------------------

crps_M2 <- matrix(nrow = nrow(prediction_summary), ncol = length(s))
for (h in seq_along(s)) {
  crps_M2[, h] <- future_sapply(as.numeric(prediction_summary$discussion),
                                function(i) {
                                  tmp_df <- df %>%
                                    filter(discussion == i) %>%
                                    filter(tau <= s[h])
                                  tmp_size <- prediction_summary %>%
                                    filter(discussion == i) %>%
                                    use_series(size)
                                  tmp_N <- nrow(tmp_df)
                                  tmp_t <- tmp_df$t
                                  tmp_beta <- tmp_df$within_discussion_parent_id
                                  tmp_a <- tmp_t[1] + s[h]
                                  tmp_b <- tmp_t[1] + a
                                  sapply(seq(1, n_samples, length.out = n_sim),
                                         function(j) {
                                           tmp_mu_imm <- tmp_M$`mu[1]`[j]
                                           tmp_psi_imm <- tmp_M$`psi[1]`[j]
                                           tmp_eta_imm <- tmp_M$`eta[1]`[j]
                                           tmp_mu_off <- tmp_M$`mu[2]`[j]
                                           tmp_psi_off <- tmp_M$`psi[2]`[j]
                                           tmp_eta_off <- tmp_M$`eta[2]`[j]
                                           simulate_gpm_cluster_process(
                                             t_seed = tmp_t, branching_structure_seed = tmp_beta,
                                             observation_horizon = tmp_a,
                                             simulation_horizon = tmp_b,
                                             immigrant_reproduction_number = tmp_mu_imm,
                                             immigrant_dispersion_parameter = tmp_psi_imm,
                                             immigrant_gi_exp_decay_rate = tmp_eta_imm,
                                             offspring_reproduction_number = tmp_mu_off,
                                             offspring_dispersion_parameter = tmp_psi_off,
                                             offspring_gi_exp_decay_rate = tmp_eta_off,
                                             refactor_output = FALSE,
                                             perform_checks = FALSE
                                           ) %>%
                                             nrow()
                                         }) %>%
                                    compute_pwm_crps(tmp_size, ., perform_checks = FALSE)
                                }, future.seed=TRUE)
}
prediction_summary <- cbind(prediction_summary, crps_M2 = crps_M2)
rm(crps_M2, tmp_M)
timestamp()
print('M_2 crps complete')

# assess M_3 --------------------------------------------------------------

tmp_M <- fit_M3
n_samples <- nrow(tmp_M) 
n_sim <- 100
alpha_ind <- grepl("alpha", colnames(tmp_M))

# assess M_3 elpd ---------------------------------------------------------

prediction_summary$elpd_M3 <-
  future_sapply(as.numeric(prediction_summary$discussion), function(i) {
    tmp_df <- df %>% 
      filter(discussion == i)
    tmp_N <- nrow(tmp_df)
    tmp_t <- tmp_df$t
    tmp_beta <- tmp_df$within_discussion_parent_id
    tmp_a <- tmp_t[1] + a
    sapply(seq(1, n_samples, length.out = n_sim),
           function(j) {
             tmp_mu <- c(tmp_M$`mu[1]`[j], rep(tmp_M$`mu[2]`[j], tmp_N-1))
             tmp_psi <- c(tmp_M$`psi[1]`[j], rep(tmp_M$`psi[2]`[j], tmp_N-1))
             tmp_eta <- c(tmp_M$`eta[1]`[j], rep(tmp_M$`eta[2]`[j], tmp_N-1))
             tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
             heterogeneous_branching_point_process_marginal_likelihood(
               t = tmp_t, beta = tmp_beta, a = tmp_a,
               mu_vec = tmp_mu, psi_vec = tmp_psi, eta_vec = tmp_eta, 
               omega = omega, alpha = tmp_alpha,
               perform_checks = F
             )
           }) %>%
      logSumExp() %>% 
      subtract(log(n_sim))
  }, future.seed=TRUE)
timestamp()
print('M_3 elpd complete')

# assess M_3 crps ---------------------------------------------------------

# predict cluster sizes
crps_M3 <- matrix(nrow = nrow(prediction_summary), ncol = length(s))
for (h in seq_along(s)) {
  crps_M3[, h] <- future_sapply(as.numeric(prediction_summary$discussion),
                                function(i) {
                                  tmp_df <- df %>%
                                    filter(discussion == i) %>%
                                    filter(tau <= s[h])
                                  tmp_size <- prediction_summary %>%
                                    filter(discussion == i) %>%
                                    use_series(size)
                                  tmp_N <- nrow(tmp_df)
                                  tmp_t <- tmp_df$t
                                  tmp_beta <- tmp_df$within_discussion_parent_id
                                  tmp_a <- tmp_t[1] + s[h]
                                  tmp_b <- tmp_t[1] + a
                                  sapply(seq(1, n_samples, length.out = n_sim),
                                         function(j) {
                                           tmp_mu_imm <- tmp_M$`mu[1]`[j]
                                           tmp_psi_imm <- tmp_M$`psi[1]`[j]
                                           tmp_eta_imm <- tmp_M$`eta[1]`[j]
                                           tmp_mu_off <- tmp_M$`mu[2]`[j]
                                           tmp_psi_off <- tmp_M$`psi[2]`[j]
                                           tmp_eta_off <- tmp_M$`eta[2]`[j]
                                           tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
                                           simulate_gpm_cluster_process(
                                             t_seed = tmp_t, branching_structure_seed = tmp_beta,
                                             observation_horizon = tmp_a,
                                             simulation_horizon = tmp_b,
                                             immigrant_reproduction_number = tmp_mu_imm,
                                             immigrant_dispersion_parameter = tmp_psi_imm,
                                             immigrant_gi_exp_decay_rate = tmp_eta_imm,
                                             offspring_reproduction_number = tmp_mu_off,
                                             offspring_dispersion_parameter = tmp_psi_off,
                                             offspring_gi_exp_decay_rate = tmp_eta_off,
                                             dominating_scalar = 2.5,
                                             sinusoid_coefficients = tmp_alpha,
                                             sinusoid_frequencies = omega,
                                             refactor_output = FALSE,
                                             perform_checks = FALSE
                                           ) %>%
                                             nrow()
                                         }) %>%
                                    compute_pwm_crps(tmp_size, ., perform_checks = FALSE)
                                }, future.seed=TRUE)
}

prediction_summary <- cbind(prediction_summary, crps_M3 = crps_M3)
rm(crps_M3, tmp_M)
timestamp()
print('M_3 crps complete')

# assess M_4 --------------------------------------------------------------

tmp_M <- fit_M4
n_samples <- nrow(tmp_M) 
n_sim <- 100
alpha_ind <- grepl("alpha", colnames(tmp_M))

# assess M_4 elpd ---------------------------------------------------------

prediction_summary$elpd_M4 <-
  future_sapply(as.numeric(prediction_summary$discussion), function(i) {
    tmp_df <- df %>% 
      filter(discussion == i)
    tmp_N <- nrow(tmp_df)
    tmp_t <- tmp_df$t
    tmp_beta <- tmp_df$within_discussion_parent_id
    tmp_z <- sapply(1:tmp_N, function(i) sum(tmp_beta == i))
    tmp_a <- tmp_t[1] + a
    sapply(seq(1, n_samples, length.out = n_sim),
           function(j) {
             tmp_mu <- c(tmp_M$`mu[1]`[j], rep(tmp_M$`mu[2]`[j], tmp_N-1))
             tmp_psi <- c(tmp_M$`psi[1]`[j], rep(tmp_M$`psi[2]`[j], tmp_N-1))
             tmp_eta <- c(tmp_M$`eta[1]`[j], rep(tmp_M$`eta[2]`[j], tmp_N-1))
             tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
             heterogeneous_branching_point_process_marginal_likelihood(
               t = tmp_t, beta = tmp_beta, a = tmp_a,
               mu_vec = tmp_mu, psi_vec = tmp_psi, eta_vec = tmp_eta, 
               omega = omega, alpha = tmp_alpha,
               z = tmp_z,
               perform_checks = FALSE
             )
           }) %>%
      logSumExp() %>% 
      subtract(log(n_sim))
  }, future.seed = T)
timestamp()
print('M_4 elpd complete')

# assess M_4 crps ---------------------------------------------------------

crps_M4 <- matrix(nrow = nrow(prediction_summary), ncol = length(s))
for (h in seq_along(s)) {
  crps_M4[, h] <- future_sapply(as.numeric(prediction_summary$discussion),
                                function(i) {
                                  tmp_df <- df %>%
                                    filter(discussion == i) %>%
                                    filter(tau <= s[h])
                                  tmp_size <- prediction_summary %>%
                                    filter(discussion == i) %>%
                                    use_series(size)
                                  tmp_N <- nrow(tmp_df)
                                  tmp_t <- tmp_df$t
                                  tmp_beta <- tmp_df$within_discussion_parent_id
                                  tmp_a <- tmp_t[1] + s[h]
                                  tmp_b <- tmp_t[1] + a
                                  sapply(seq(1, n_samples, length.out = n_sim),
                                         function(j) {
                                           tmp_mu_imm <- tmp_M$`mu[1]`[j]
                                           tmp_psi_imm <- tmp_M$`psi[1]`[j]
                                           tmp_eta_imm <- tmp_M$`eta[1]`[j]
                                           tmp_mu_off <- tmp_M$`mu[2]`[j]
                                           tmp_psi_off <- tmp_M$`psi[2]`[j]
                                           tmp_eta_off <- tmp_M$`eta[2]`[j]
                                           tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
                                           simulate_gpm_cluster_process(
                                             t_seed = tmp_t, branching_structure_seed = tmp_beta,
                                             observation_horizon = tmp_a,
                                             simulation_horizon = tmp_b,
                                             immigrant_reproduction_number = tmp_mu_imm,
                                             immigrant_dispersion_parameter = tmp_psi_imm,
                                             immigrant_gi_exp_decay_rate = tmp_eta_imm,
                                             offspring_reproduction_number = tmp_mu_off,
                                             offspring_dispersion_parameter = tmp_psi_off,
                                             offspring_gi_exp_decay_rate = tmp_eta_off,
                                             dominating_scalar = 2.5,
                                             sinusoid_coefficients = tmp_alpha,
                                             sinusoid_frequencies = omega,
                                             refactor_output = FALSE,
                                             perform_checks = FALSE
                                           ) %>%
                                             nrow()
                                         }) %>%
                                    compute_pwm_crps(tmp_size, ., perform_checks = FALSE)
                                }, future.seed=TRUE)
}
prediction_summary <- cbind(prediction_summary, crps_M4 = crps_M4)
rm(crps_M4, tmp_M)
timestamp()
print('M_4 crps complete')

# assess M_5 --------------------------------------------------------------

tmp_M <- fit_M5
n_samples <- nrow(tmp_M) 
n_sim <- 100
alpha_ind <- grepl("alpha", colnames(tmp_M))

# assess M_5 elpd ---------------------------------------------------------

prediction_summary$elpd_M5 <-
  future_sapply(as.numeric(prediction_summary$discussion), function(i) {
    tmp_df <- df %>% 
      filter(discussion == i)
    tmp_N <- nrow(tmp_df)
    tmp_t <- tmp_df$t
    tmp_beta <- tmp_df$within_discussion_parent_id
    tmp_a <- tmp_t[1] + a
    sapply(seq(1, n_samples, length.out = n_sim),
           function(j) {
             tmp_mu <- c(tmp_M$`mu[1]`[j], rep(tmp_M$`mu[2]`[j], tmp_N-1))
             tmp_psi <- c(tmp_M$`psi[1]`[j], rep(tmp_M$`psi[2]`[j], tmp_N-1))
             tmp_eta <- c(tmp_M$`eta[1]`[j], rep(tmp_M$`eta[2]`[j], tmp_N-1))
             tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
             heterogeneous_branching_point_process_marginal_likelihood(
               t = tmp_t, beta = tmp_beta, a = tmp_a,
               mu_vec = tmp_mu, psi_vec = tmp_psi, eta_vec = tmp_eta, 
               omega = omega, alpha = tmp_alpha,
               perform_checks = F
             )
           }) %>%
      logSumExp() %>% 
      subtract(log(n_sim))
  }, future.seed=TRUE)
timestamp()
print('M_5 elpd complete')

# assess M_5 crps ---------------------------------------------------------

crps_M5 <- matrix(nrow = nrow(prediction_summary), ncol = length(s))
for (h in seq_along(s)) {
  crps_M5[, h] <- future_sapply(as.numeric(prediction_summary$discussion),
                                function(i) {
                                  tmp_df <- df %>%
                                    filter(discussion == i) %>%
                                    filter(tau <= s[h])
                                  tmp_size <- prediction_summary %>%
                                    filter(discussion == i) %>%
                                    use_series(size)
                                  tmp_N <- nrow(tmp_df)
                                  tmp_t <- tmp_df$t
                                  tmp_beta <- tmp_df$within_discussion_parent_id
                                  tmp_a <- tmp_t[1] + s[h]
                                  tmp_b <- tmp_t[1] + a
                                  sapply(seq(1, n_samples, length.out = n_sim),
                                         function(j) {
                                           tmp_mu_imm <- tmp_M$`mu[1]`[j]
                                           tmp_psi_imm <- tmp_M$`psi[1]`[j]
                                           tmp_eta_imm <- tmp_M$`eta[1]`[j]
                                           tmp_mu_off <- tmp_M$`mu[2]`[j]
                                           tmp_psi_off <- tmp_M$`psi[2]`[j]
                                           tmp_eta_off <- tmp_M$`eta[2]`[j]
                                           tmp_alpha <- as.numeric(tmp_M[j, alpha_ind])
                                           simulate_gpm_cluster_process(
                                             t_seed = tmp_t, branching_structure_seed = tmp_beta,
                                             observation_horizon = tmp_a,
                                             simulation_horizon = tmp_b,
                                             immigrant_reproduction_number = tmp_mu_imm,
                                             immigrant_dispersion_parameter = tmp_psi_imm,
                                             immigrant_gi_exp_decay_rate = tmp_eta_imm,
                                             offspring_reproduction_number = tmp_mu_off,
                                             offspring_dispersion_parameter = tmp_psi_off,
                                             offspring_gi_exp_decay_rate = tmp_eta_off,
                                             dominating_scalar = 2.5,
                                             sinusoid_coefficients = tmp_alpha,
                                             sinusoid_frequencies = omega,
                                             refactor_output = FALSE,
                                             perform_checks = FALSE
                                           ) %>%
                                             nrow()
                                         }) %>%
                                    compute_pwm_crps(tmp_size, ., perform_checks = FALSE)
                                }, future.seed=TRUE)
}
prediction_summary <- cbind(prediction_summary, crps_M5 = crps_M5)
rm(crps_M5, tmp_M)
timestamp()
print('M_5 crps complete')

# save model assessment ---------------------------------------------------

saveRDS(prediction_summary, 'predictions/prediction_summary.RDS')
print('assessment complete')

