library(dplyr)
library(tidyr)
library(purrr)
library(tictoc)
library(future)
plan(multisession, workers = 10L)


source("simulation-step-bootstrap/2_estimation.R")
source("simulation-step-bootstrap/3_performance_criteria.R")
source("simulation-step-bootstrap/4_simulation_driver.R")

load("simulation-step-bootstrap/wwc_es.RData")

datasets <- run_sim(
  iterations = 2000,
  mean_smd = 0.2,
  tau = 0.45,
  omega = 0,
  cor_mu = 0.8,
  cor_sd = 0.05,
  weights = 1,
  m = 30,
  n_multiplier = 1,
  steps = 0.025,
  bootstrap = "none",
  seed = 20312448,
  stepfun_methods = NULL,
  comparison_methods = NULL,
  summarize_performance = FALSE
)

comp_ests <- 
  datasets %>% 
  nest(.by = rep) %>%
  mutate(
    res = future_map(data, estimate_comparison_methods, .progress = TRUE)
  ) %>%
  select(-data) %>%
  unnest(res)


test_res <- run_sim(
  iterations = 20,
  mean_smd = 0.2,
  tau = 0.45,
  omega = 0,
  cor_mu = 0.8,
  cor_sd = 0.05,
  weights = 1,
  m = 30,
  n_multiplier = 1,
  steps = 0.025,
  bootstrap = "multinomial",
  seed = 20312448,
  comparison_methods = NULL,
  summarize_performance = FALSE
)

test_res$res
