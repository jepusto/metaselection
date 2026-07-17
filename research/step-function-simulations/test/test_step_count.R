library(tidyverse)
library(future)
library(furrr)
library(metafor)
library(metaselection)
plan(multisession)


load("research/step-function-simulations/wwc_es.RData") # do we need to change this?
source("research/step-function-simulations/1_estimation.R")
source("research/step-function-simulations/2_performance_criteria.R")
source("research/step-function-simulations/3_simulation_driver.R")


res0 <- run_sim(
  iterations = 100,
  mean_smd = 0.0,
  tau = 0.15,
  omega = 0,
  cor_mu = 0.8,
  cor_sd = 0.05,
  weight = 0.1,
  psi = 0,
  m = 60,
  n_multiplier = 3,
  steps = 0.025,
  bootstrap = "none",
  seed = 20312448,
  summarize_performance = TRUE
)

res0 %>%
  filter(param == "beta") %>%
  select(estimator, bias, var, rmse, coverage, width)

res1 <- run_sim(
  iterations = 100,
  mean_smd = 0.0,
  tau = 0.15,
  omega = 0,
  cor_mu = 0.8,
  cor_sd = 0.05,
  weight = 0.1,
  psi = 1,
  m = 60,
  n_multiplier = 3,
  steps = 0.025,
  bootstrap = "none",
  seed = 20312448,
  summarize_performance = TRUE
)

res1 %>%
  filter(param == "beta") %>%
  select(estimator, bias, var, rmse, coverage, width)

left_join(res0, res1, by = c("estimator","param")) %>%
  filter(param == "beta") %>%
  select(estimator, starts_with("bias."), starts_with("rmse."), starts_with("coverage."))
