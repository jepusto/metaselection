library(tidyverse)
library(boot)
library(tictoc)
library(furrr)
library(metaselection)

load("simulation-step-paper/wwc_es.RData") # do we need to change this?
source("simulation-step-paper/2_estimation.R")

mean_smd <- 0.8
tau <- 0.05
omega <- 0.00
m <- 60
cor_mu <- .8
cor_sd <- 0.05
censor_fun <- step_fun
wwc_es$n <- round(wwc_es$n / 3)
n_ES_sim <- n_ES_empirical(wwc_es)

censor_fun_param <- censor_fun(cut_vals = .025, weights = .05)

set.seed(20309299)

dat <- r_meta(mean_smd = mean_smd,
              tau = tau,
              omega = omega,
              m = m,
              cor_mu = cor_mu,
              cor_sd = cor_sd,
              censor_fun =  censor_fun_param,  
              n_ES_sim = n_ES_sim)

step_est <- estimate_step_models(
  dat = dat, 
  steps = .025, 
  estimators = c("step-MLE"),
  bootstrap = "multinomial",
  CML_optimizer = "BFGS"
)
step_est

dat <- r_meta(mean_smd = mean_smd,
              tau = tau,
              omega = omega,
              m = m,
              cor_mu = cor_mu,
              cor_sd = cor_sd,
              censor_fun =  censor_fun_param,  
              n_ES_sim = n_ES_sim)

debugonce(metaselection::selection_model)
estimate_step_models(
  dat = dat, 
  steps = .025, 
  estimators = c("step-MLE"),
  bootstrap = "multinomial",
  CML_optimizer = "BFGS",
  CI_type = c("large-sample","basic","percentile","bias-corrected","BCa"),
)

estimate_comparison_methods(dat = dat)

source("simulation-step-paper/3_performance_criteria.R")
source("simulation-step-paper/4_simulation_driver.R")

steps <- .025
weights <- .05
n_multiplier <- 1/3
censor_fun <- step_fun

censor_fun_param <- censor_fun(cut_vals = steps, weights = weights)
true_params <- data.frame(param = c("beta", "gamma", "zeta1"))
true_params$true_param <- c(mean_smd, log(tau^2 + omega^2), log(censor_fun_param(steps + 1e-8)))

test_res <- run_sim(
  iterations = 40L,
  mean_smd = mean_smd,
  tau = tau,
  omega = omega,
  m = m,
  cor_mu = cor_mu,
  cor_sd = cor_sd,
  steps = steps,
  weights = weights,
  n_multiplier = n_multiplier,
  stepfun_methods = c("step-MLE", "step-hybrid"),
  comparison_methods = "All",
  summarize_performance = FALSE,
  bootstrap = "none"
) |>
  merge(true_params)

View(test_res)
perf <- calc_performance(test_res)
View(perf)
perf_winz <- calc_performance(test_res, winz = 2.5)
View(perf_winz)

# 
# test_res 
# 
# sim_dat <- run_sim(
#   iterations = 100L, 
#   mean_smd = mean_smd,
#   tau = tau,
#   omega = omega,
#   m = m, 
#   cor_mu = cor_mu,
#   cor_sd = cor_sd, 
#   censor_fun = step_fun,
#   n_ES_sim = n_ES_empirical(wwc_es),
#   steps = steps,
#   weights = weights,
#   stepfun_methods = NULL,
#   comparison_methods = NULL
# )


# check the 4psm ----------------------------------------------------------

iterations <- 5
steps <- c(.025, .5)
weights <- c(.2, .3)

censor_fun_param <- censor_fun(cut_vals = steps, weights = weights)


true_params <- data.frame(param = c("beta", "gamma", "zeta1", "zeta2"))
true_params$true_param <- c(mean_smd, log(tau^2 + omega^2), log(censor_fun_param(steps + 1e-8)))

results <- run_sim(
  iterations = iterations, 
  mean_smd = mean_smd,
  tau = tau,
  omega = omega,
  m = m, 
  cor_mu = cor_mu,
  cor_sd = cor_sd, 
  censor_fun = step_fun,
  estimate_4psm = TRUE, 
  steps = steps,
  summarize_performance = FALSE,
  filter_n_study = 2L,
  stepfun_methods = c("step-MLE", "step-hybrid-Broyden"),
  ML_optimizer_control = list(all.methods = TRUE),
  comparison_methods = "FEC",
  bootstrap = NULL
) %>%
  merge(true_params)

calc_performance(results)
calc_performance(results, filter_n_study = 2L)

results <- run_sim(
  iterations = iterations, 
  mean_smd = mean_smd,
  tau = tau,
  omega = omega,
  m = m, 
  cor_mu = cor_mu,
  cor_sd = cor_sd, 
  censor_fun = step_fun,
  estimate_4psm = FALSE, 
  steps = steps,
  summarize_performance = FALSE,
  filter_n_study = 2L,
  stepfun_methods = c("step-MLE", "step-hybrid-Broyden"),
  ML_optimizer_control = list(all.methods = TRUE),
  comparison_methods = "FEC",
  R = 9
) %>%
  merge(true_params)


calc_performance(results)
calc_performance(results, filter_n_study = 2L)


results <- run_sim(
  iterations = iterations, 
  mean_smd = mean_smd,
  tau = tau,
  omega = omega,
  m = m, 
  cor_mu = cor_mu,
  cor_sd = cor_sd, 
  censor_fun = step_fun,
  estimate_4psm = FALSE, 
  steps = steps,
  summarize_performance = TRUE,
  filter_n_study = 2L,
  stepfun_methods = c("step-MLE", "step-hybrid-Broyden"),
  ML_optimizer_control = list(all.methods = TRUE),
  comparison_methods = "FEC",
  bootstrap = NULL
)
