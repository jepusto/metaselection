#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Source packages and functions ----

library(tidyverse)
library(simhelpers)
library(future)
library(metafor)

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
  m = 30L,
  cor_mu = 0.8
)
  
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Define experimental design parameters ----

all_params <- 
  readRDS("research/step-function-simulations/simulation_parameters.rds") 

no_selection_params <- 
  all_params %>%
  filter(
    weights == 1, 
    bootstrap == "none",
    mean_smd %in% c(0.0,0.8)
  ) %>% 
  select(iterations, mean_smd, tau, omega, cor_mu, n_multiplier, m, seed) %>%
  mutate(iterations = round(4800 / sqrt(m)))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Run simulations ----

plan(multisession, workers = 12L)

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

univariate_RE_results <- 
  readRDS(file = "research/step-function-simulations/univariate-RE-simulation-results.rds") %>%
  mutate(
    mu_fac = fct(as.character(mean_smd)),
    tau_fac = fct(as.character(tau), levels = c("0.05","0.15","0.3","0.45")),
    het_ratio = omega^2 / tau^2,
    het_ratio = as.character(het_ratio),
    var_fac = paste(het_ratio, tau_fac),
    cor_fac = fct(as.character(cor_mu)),
    N_factor = if_else(n_multiplier < 1, "Small", "Typical"),
    J = as.character(m),
    J = factor(J, levels = c("15", "30", "60", "90", "120")),
    J_inv_sqrt = 1 / sqrt(m),
    total_var = tau^2 + omega^2,
    relbias = bias / total_var,
    relbias_mcse = bias_mcse / total_var
  )

bias_fit <- 
  rma.uni(
    yi = bias,
    sei = bias_mcse,
    mods = ~ cor_fac : tau_fac : het_ratio : N_factor : J_inv_sqrt,
    data = univariate_RE_results
  ) 
bias_blup_est <- blup(bias_fit)

univariate_RE_blups <-
  univariate_RE_results %>%
  mutate(
    bias_blup = bias_blup_est$pred,
    bias_blup_mcse = bias_blup_est$se,
    relbias_blup = bias_blup / total_var,
    relbias_blup_mcse = bias_blup_mcse / total_var
  )

# MCSE of raw bias

ggplot(univariate_RE_results) +
  aes(J, bias_mcse, color = tau_fac, linetype = het_ratio, shape = het_ratio) + 
  geom_hline(yintercept = 0) + 
  geom_point() + geom_line(aes(group = var_fac)) + 
  facet_grid(
    mean_smd + N_factor ~ cor_mu, 
    scales = "free_y",
    labeller = label_bquote(
      cols = rho == .(cor_mu),
      rows = atop( mu == .(mean_smd), .(N_factor) ~ "sample size")
    ),
  ) +
  labs(
    x = "Number of studies (J)", 
    y = "MCSE of Bias", 
    color = expression(tau[B]),
    shape = "Heterogeneity ratio",
    linetype = "Heterogeneity ratio"
  ) + 
  theme_bw() +
  theme(legend.position = "right")


# raw bias

ggplot(univariate_RE_blups) +
  aes(J, bias, color = tau_fac, linetype = het_ratio, shape = het_ratio) + 
  geom_hline(yintercept = 0) + 
  geom_point() + geom_line(aes(group = var_fac)) + 
  facet_grid(
    mean_smd + N_factor ~ cor_mu, 
    scales = "free_y",
    labeller = label_bquote(
      cols = rho == .(cor_mu),
      rows = atop( mu == .(mean_smd), .(N_factor) ~ "sample size")
    ),
  ) +
  labs(
    x = "Number of studies (J)", 
    y = "Bias", 
    color = expression(tau[B]),
    shape = "Heterogeneity ratio",
    linetype = "Heterogeneity ratio"
  ) + 
  theme_bw() +
  theme(legend.position = "right")


# model-adjusted bias

ggplot(univariate_RE_blups) +
  aes(J, bias_blup, color = tau_fac, linetype = het_ratio, shape = het_ratio) + 
  geom_hline(yintercept = 0) + 
  geom_point() + geom_line(aes(group = var_fac)) + 
  facet_grid(
    mean_smd + N_factor ~ cor_mu, 
    scales = "free_y",
    labeller = label_bquote(
      cols = rho == .(cor_mu),
      rows = atop( mu == .(mean_smd), .(N_factor) ~ "sample size")
    ),
  ) +
  labs(
    x = "Number of studies (J)", 
    y = "Bias", 
    color = expression(tau[B]),
    shape = "Heterogeneity ratio",
    linetype = "Heterogeneity ratio"
  ) + 
  theme_bw() +
  theme(legend.position = "right")

ggsave(
  "research/step-function-simulations/univariate-RE-simulation.bmp", 
  width = 8, 
  height = 6,
  dpi = 600
)
