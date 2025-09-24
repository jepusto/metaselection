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
mean_smd <- 0.8
tau <- 0.05
omega <- 0.00
m <- 60
cor_mu <- .4
cor_sd <- 0.05
wwc_es$n <- round(wwc_es$n / 3)
n_ES_sim <- n_ES_empirical(wwc_es)
weights <- 0.5
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
    CML = future_map(
      datasets, 
      estimate_step_models, 
      estimators = "CML", steps = 0.025, 
      priors = default_priors(), 
      bootstrap = "none",
      .progress = TRUE
    )
  )

CML_res <- 
  data_res %>%
  select(-data) %>%
  unnest(CML) %>%
  left_join(true_params, by = "param") %>%
  mutate(
    model = "3PSM"
  )

calc_performance(CML_res)

ggplot(CML_res) +
  facet_wrap(~ param, scales = "free") + 
  aes(Est) + 
  geom_density()

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

CML_wide %>% 
  filter(beta < 0) %>%
  mutate(
    tau = exp(gamma / 2),
    lambda = exp(zeta1)
  )

dat <- filter(data_res, rep == 626)$data[[1]]

est <- selection_model(
  yi = d,
  sei = sd_d,
  dat = dat,
  estimators = "CML", 
  steps = 0.025, 
  priors = default_priors(), 
  bootstrap = "none",
  use_jac = FALSE,
  # vcov_type = "raw",
  optimizer = c("Rvmmin","Nelder-Mead","nlminb")
)
est
theta <- est$est$Est
debugonce(step_loglik)
step_loglik(theta = theta, yi = dat$d, sei = dat$sd_d, priors = NULL)
prior_spec <- default_priors()

step_loglik(theta = theta, yi = dat$d, sei = dat$sd_d, priors = default_priors())

step_score(theta = theta, yi = dat$d, sei = dat$sd_d, priors = default_priors())

rma.uni(yi = d, sei = sda, data = dat) %>% 
  funnel(level=95, shade="gray55", refline=0, legend=FALSE)
rma.uni(yi = d, sei = sda, data = dat) %>%
  selmodel(type = "step", steps = .025)
