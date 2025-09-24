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

# Set parameters
mean_smd <- 0.2
tau <- 0.45
omega <- 0.00
m <- 15
cor_mu <- .4
cor_sd <- 0.05
wwc_es$n <- round(wwc_es$n / 3)
n_ES_sim <- n_ES_empirical(wwc_es)
weights <- 0.2
censor_fun_param <- step_fun(cut_vals = .025, weights = weights)

true_params <- data.frame(
  param = c("beta", "gamma", "zeta1"),
  true_param = c(mean_smd, log(tau^2), log(weights))
)


# Generate datasets
datasets <- future_map(
  1:1000, 
  ~ r_meta(mean_smd = mean_smd,
           tau = tau,
           omega = omega,
           m = m,
           cor_mu = cor_mu,
           cor_sd = cor_sd,
           censor_fun =  censor_fun_param,  
           n_ES_sim = n_ES_sim),
  .id = "rep",
  .options = furrr_options(seed = NULL),
  .progress = TRUE
)

# Fit selection models

data_res <- 
  tibble(
    rep = 1:length(datasets),
    data = datasets,
    CML_default = future_map(
      datasets, 
      estimate_step_models, 
      estimators = "CML", 
      steps = 0.025, 
      priors = default_priors(), 
      bootstrap = "none",
      use_jac = FALSE,
      .progress = TRUE
    ),
  )

CML_res <- 
  data_res %>%
  select(-data) %>%
  unnest(CML_default) %>%
  left_join(true_params, by = "param") %>%
  mutate(
    model = "3PSM"
  )

calc_performance(CML_res)

ggplot(CML_res) +
  aes(color = param, fill = param) + 
  facet_wrap(~ param, scales = "free") + 
  geom_vline(aes(xintercept = true_param)) + 
  geom_density(aes(Est), alpha = 0.5)

CML_wide <- 
  CML_res %>%
  select(rep, param, Est) %>%
  pivot_wider(names_from = param, values_from = Est)

ggplot(CML_wide) +
  aes(beta, zeta1) + 
  geom_point()
ggplot(CML_wide) +
  aes(beta, gamma) + 
  geom_point()
ggplot(CML_wide) +
  aes(zeta1, gamma) + 
  geom_point()

dat <- data_res$data[[1]]

estimate_step_models(
  dat = dat, 
  steps = 0.025, 
  estimators = "CML",
  priors = default_priors()
)

est <- selection_model(
  yi = d,
  sei = sd_d,
  dat = dat,
  estimators = "CML", 
  steps = 0.025, 
  priors = default_priors(), 
  bootstrap = "none",
  use_jac = TRUE,
  # vcov_type = "raw",
  optimizer = c("Rvmmin","Nelder-Mead","nlminb")
)
est
theta <- est$est$Est
step_loglik(theta = theta, yi = dat$d, sei = dat$sd_d, priors = NULL)
prior_spec <- default_priors()
debugonce(parse_step_params)
step_loglik(theta = theta, yi = dat$d, sei = dat$sd_d, priors = default_priors())

step_score(theta = theta, yi = dat$d, sei = dat$sd_d, priors = default_priors())

rma.uni(yi = d, sei = sda, data = dat) %>% 
  funnel(level=95, shade="gray55", refline=0, legend=FALSE)
rma.uni(yi = d, sei = sda, data = dat) %>%
  selmodel(type = "step", steps = .025)
