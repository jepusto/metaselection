#-------------------------------------------------------------------------------
# Source packages and functions
library(dplyr)
library(tidyr)
library(purrr)
library(tictoc)


source("R/data-generation.R")
source("R/step-utilities.R")
source("R/step-function.R")
source("R/hybrid-estimating.R")
source("R/selection_model.R")
source("simulation-step-bootstrap/2_estimation.R")
source("simulation-step-bootstrap/3_performance_criteria.R")
source("simulation-step-bootstrap/4_simulation_driver.R")

dat_big <- run_sim(
  iterations = 1, 
  mean_smd = 0.2,
  tau = 0.1,
  omega = 0.05,
  m = 1000, 
  cor_mu = 0.4,
  cor_sd = 0.05, 
  censor_fun = step_fun,
  n_ES_dat = wwc_es,
  steps = 0.025,
  weights = 1,
  summarize_performance = FALSE,
  filter_n_study = FALSE,
  stepfun_methods = NULL,
  comparison_methods = NULL
)

dat_small <- run_sim(
  iterations = 1, 
  mean_smd = 0.2,
  tau = 0.1,
  omega = 0.05,
  m = 1000, 
  cor_mu = 0.4,
  cor_sd = 0.05, 
  censor_fun = step_fun,
  n_ES_dat = wwc_es,
  n_multiplier = 1/3,
  steps = 0.025,
  weights = 1,
  summarize_performance = FALSE,
  filter_n_study = FALSE,
  stepfun_methods = NULL,
  comparison_methods = NULL
)

sample_size_dat <- 
  bind_rows(big = dat_big, small = dat_small, .id = "n_multiplier")

ggplot(sample_size_dat, aes(n, fill = n_multiplier)) + 
  geom_density(alpha = 0.3)

sample_size_dat %>%
  group_by(n_multiplier) %>%
  reframe(p = seq(0,100,20), N = quantile(n, seq(0,1,.2))) %>%
  pivot_wider(names_from = n_multiplier, values_from = N)
