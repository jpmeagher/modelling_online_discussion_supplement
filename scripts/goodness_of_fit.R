
# packages ----------------------------------------------------------------

# devtools::install_github("jpmeagher/onlineMessageboardActivity")
library(onlineMessageboardActivity)
library(dplyr)
library(lubridate)
library(magrittr)
library(future.apply)

# future ------------------------------------------------------------------

plan("multisession", workers = 8)

# data --------------------------------------------------------------------

# observation interval
a <- 48
df <- messageboard_df %>% 
  mutate(t = as.numeric(difftime(time, dmy_hms(010419000000, tz = "Europe/London"), units = "hours"))) %>%
  group_by(discussion) %>%
  mutate(tau = t - min(t)) %>%
  filter(tau < a)


# models ------------------------------------------------------------------

fit_M1 <- readRDS("models/fit_M1.RDS") %>% as.data.frame()
fit_M2 <- readRDS("models/fit_M2.RDS") %>% as.data.frame()
fit_M3 <- readRDS("models/fit_M3.RDS") %>% as.data.frame()
fit_M4 <- readRDS("models/fit_M4.RDS") %>% as.data.frame()
fit_M5 <- readRDS("models/fit_M5.RDS") %>% as.data.frame()

# discussion size ---------------------------------------------------------

discussion_size <- df %>% 
  group_by(discussion) %>% 
  summarise(n = length(discussion)) %>% 
  ungroup() %>% 
  full_join(
    df %>% 
      filter(parent_id == 0) %>% 
      select(discussion, time) %>% 
      mutate(t = as.numeric(difftime(time, dmy_hms(010419000000, tz = "Europe/London"), units = "hours"))) %>% 
      mutate(tod = mod(t, 24))
  )


# hyperparameters ---------------------------------------------------------

# sinusoidal basis
K <- 2
day <- 24
omega <- 2 * pi * (1:K / day)

# simulate discussions ----------------------------------------------------

set.seed(101)
timestamp()
discussion_size <- discussion_size %>% 
  mutate(M1 = future_sapply(
    t, function(x){
      par <- fit_M1 %>% as.data.frame() %>% sample_n(1)
      simulate_gpm_cluster_process(
        t_seed = x, branching_structure_seed = 0,
        observation_horizon = x, simulation_horizon = x + a,
        immigrant_reproduction_number = par$`mu[1]`,
        immigrant_gi_exp_decay_rate = par$`eta[1]`,
        refactor_output = FALSE
      ) %>% 
        nrow()
    }, future.seed=TRUE
  )) %>% 
  mutate(M2 = future_sapply(
    t, function(x){
      par <- fit_M2 %>% as.data.frame() %>% sample_n(1)
      simulate_gpm_cluster_process(
        t_seed = x, branching_structure_seed = 0,
        observation_horizon = x, simulation_horizon = x + a,
        immigrant_reproduction_number = par$`mu[1]`,
        immigrant_gi_exp_decay_rate = par$`eta[1]`,
        offspring_reproduction_number = par$`mu[2]`,
        offspring_gi_exp_decay_rate = par$`eta[2]`,
        refactor_output = FALSE
      ) %>% 
        nrow()
    }, future.seed=TRUE
  )) %>%  
  mutate(M3 = future_sapply(
    t, function(x){
      par <- fit_M3 %>% as.data.frame() %>% sample_n(1)
      alpha <- par %>% 
        select(starts_with("alpha")) %>% 
        unlist()
      simulate_gpm_cluster_process(
        t_seed = x, branching_structure_seed = 0,
        observation_horizon = x, simulation_horizon = x + a,
        immigrant_reproduction_number = par$`mu[1]`,
        immigrant_gi_exp_decay_rate = par$`eta[1]`,
        offspring_reproduction_number = par$`mu[2]`,
        offspring_gi_exp_decay_rate = par$`eta[2]`,
        sinusoid_coefficients = alpha, sinusoid_frequencies = omega,
        dominating_scalar = 2.5,
        refactor_output = FALSE
      ) %>% 
        nrow()
    }, future.seed=TRUE
  )) %>%  
  mutate(M4 = future_sapply(
    t, function(x){
      par <- fit_M4 %>% as.data.frame() %>% sample_n(1)
      alpha <- par %>% 
        select(starts_with("alpha")) %>% 
        unlist()
      simulate_gpm_cluster_process(
        t_seed = x, branching_structure_seed = 0,
        observation_horizon = x, simulation_horizon = x + a,
        immigrant_reproduction_number = par$`mu[1]`,
        immigrant_gi_exp_decay_rate = par$`eta[1]`,
        immigrant_dispersion_parameter = par$`psi[1]`,
        offspring_reproduction_number = par$`mu[2]`,
        offspring_gi_exp_decay_rate = par$`eta[2]`,
        offspring_dispersion_parameter = par$`psi[2]`,
        sinusoid_coefficients = alpha, sinusoid_frequencies = omega,
        dominating_scalar = 2.5,
        refactor_output = FALSE
      ) %>% 
        nrow()
    }, future.seed=TRUE
  )) %>% 
  mutate(M5 = future_sapply(
    t, function(x){
      par <- fit_M5 %>% as.data.frame() %>% sample_n(1)
      alpha <- par %>% 
        select(starts_with("alpha")) %>% 
        unlist()
      simulate_gpm_cluster_process(
        t_seed = x, branching_structure_seed = 0,
        observation_horizon = x, simulation_horizon = x + a,
        immigrant_reproduction_number = par$`mu[1]`,
        immigrant_gi_exp_decay_rate = par$`eta[1]`,
        immigrant_dispersion_parameter = par$`psi[1]`,
        offspring_reproduction_number = par$`mu[2]`,
        offspring_gi_exp_decay_rate = par$`eta[2]`,
        sinusoid_coefficients = alpha, sinusoid_frequencies = omega,
        dominating_scalar = 2.5,
        refactor_output = FALSE
      ) %>% 
        nrow()
    }, future.seed=TRUE
  ))
timestamp()

saveRDS(discussion_size, "predictions/discussion_size.RDS")
