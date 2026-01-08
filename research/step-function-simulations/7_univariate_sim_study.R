#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Source packages and functions ----

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(simhelpers)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Load empirical data used in simulation design ----

load("research/step-function-simulations/wwc_es.RData")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Simulation components ----

# data-generating function

gen_CHE_meta <- function(
    mean_smd, 
    tau = 0, 
    omega = 0, 
    m = 15, 
    cor_mu = 0.8, 
    cor_sd = 0.05,
    n_ES_sim = metaselection::n_ES_empirical(wwc_es)
) {
  metaselection::r_meta(
    mean_smd = mean_smd, tau = tau, omega = omega,
    m = m, cor_mu = cor_mu, cor_sd = cor_sd,
    censor_fun = metaselection::step_fun(),
    n_ES_sim = n_ES_sim,
    m_multiplier = 1L
  )
}


# Estimation function

estimate_RE <- function(dat) {
  
  RE_mod <- metafor::rma.uni(
    yi = d, vi = var_d, 
    data = dat, 
    method = c("REML","PM","HE","DL")
  )
  
  data.frame(
    beta = as.numeric(RE_mod$beta),
    tausq = RE_mod$tau2
  )
  
}


# Performance assessement
eval_tausq <- function(res, tau = 0, omega = 0) {
  tausq_total <- tau^2 + omega^2
  simhelpers::calc_absolute(
    res, 
    estimates = tausq, 
    true_param = tausq_total, 
    criteria = c("bias","rmse")
  )
}

# Simulation driver

univariate_sim <- bundle_sim(
  f_generate = gen_CHE_meta,
  f_analyze = estimate_RE,
  f_summarize = eval_tausq,
  reps_name = "iterations"
)

univariate_sim(
  iterations = 100L, 
  mean_smd = 0, 
  tau = 0.05, 
  omega = 0,
  m = 120L,
  cor_mu = 0.8
)
  
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Define experimental design parameters ----

all_params <- 
  readRDS("research/step-function-simulations/simulation_parameters.rds") 

no_selection_params <- 
  all_params %>%
  filter(weights == 1, bootstrap == "none", n_multiplier == 1) %>% 
  select(iterations, mean_smd, tau, omega, cor_mu, m, seed) %>%
  mutate(iterations = 100L)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Run simulations ----

library(future)
plan(multisession, workers = 10L)

univariate_RE_results <- 
  evaluate_by_row(
    params = no_selection_params, 
    sim_function = univariate_sim,
    .progress = TRUE
  )

plan(sequential)

saveRDS(univariate_RE_results, file = "research/step-function-simulations/univariate-RE-simulation-results.rds")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Create figure showing bias of heterogeneity estimator ----

univariate_RE_results <- readRDS(file = "research/step-function-simulations/univariate-RE-simulation-results.rds")