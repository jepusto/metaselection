#-------------------------------------------------------------------------------
# Source packages and functions
library(readr)
library(dplyr)
library(tidyr)


#-------------------------------------------------------------------------------
# Experimental Design

# Express the simulation parameters as vectors/lists

design_factors <- list(
  mean_smd = c(0.0, 0.2, 0.4, 0.8), # average effect size
  tau = c(0.05, 0.15, 0.30, 0.45), # between study heterogeneity
  het_ratio = c(0, 0.5),
  cor_mu = c(0.4, 0.8), # average correlation between outcomes         
  cor_sd = c(0.05), # sd correlation between outcomes
  weights = c(0.02, 0.05, 0.1, 0.2, 0.5, 1), # weights
  m = c(120, 90, 60, 30, 15),	# number of studies in each meta analysis
  n_multiplier = c(1/3, 1),
  batch = 1L,
  stepfun_methods = list(c("CML","CML-model","ARGL")),
  priors = "Weak"
)

params <- 
  expand_grid(!!!design_factors) |>
  mutate(
    omega = tau * sqrt(het_ratio),
    steps = 0.025,
    bootstrap = "none",
    comparison_methods = "All",
    iterations = 2400L,
    summarize_performance = TRUE
  ) %>%
  select(-het_ratio)


bootstrap_factors <- list(
  mean_smd = c(0.0, 0.2, 0.4, 0.8), # average effect size
  tau = c(0.05, 0.30, 0.45), # between study heterogeneity
  het_ratio = c(0, 0.5),
  cor_mu = 0.8, # average correlation between outcomes         
  cor_sd = 0.05, # sd correlation between outcomes
  weights = c(0.05, 0.2, 1), # weights
  m = c(60, 30, 15),	# number of studies in each meta analysis
  n_multiplier = c(1/3, 1),
  batch = 1:10,
  stepfun_methods = list(c("CML","ARGL")),
  priors = "Weak",
  bootstrap = c("multinomial","two-stage","exponential"),
  R = list(c(49,99,199,299))
)

bootstrap_params <- 
  expand_grid(!!!bootstrap_factors) %>%
  mutate(
    omega = tau * sqrt(het_ratio),
    steps = 0.025,
    comparison_methods = "None",
    iterations = 240L,
    summarize_performance = FALSE
  ) %>%
  select(-het_ratio)


big_B_factors <- list(
  mean_smd = c(0.2), # average effect size
  tau = c(0.05, 0.45), # between study heterogeneity
  het_ratio = c(0.5),
  cor_mu = c(0.8), # average correlation between outcomes         
  cor_sd = c(0.05), # sd correlation between outcomes
  weights = c(0.05, 0.2), # weights
  m = c(15),	# number of studies in each meta analysis
  n_multiplier = c(1),
  batch = 1:100,
  stepfun_methods = list(c("CML","ARGL")),
  priors = "Weak",
  bootstrap = "two-stage",
  R = list(c(49,99,199,299,1999))
)

big_B_params <- 
  expand_grid(!!!big_B_factors) |>
  mutate(
    omega = tau * sqrt(het_ratio),
    steps = 0.025,
    comparison_methods = "none",
    iterations = 24L,
    summarize_performance = FALSE,
  ) %>%
  select(-het_ratio)

all_params <- 
  bootstrap_params %>%
  anti_join(big_B_params, by = c("mean_smd","tau","cor_mu","cor_sd", "weights", "m", "n_multiplier","omega","bootstrap")) %>%
  bind_rows(big_B_params) %>%
  bind_rows(params) %>%
  mutate(
    row = row_number(),
    seed = 20250918L + row
  )

all_params %>%
  group_by(bootstrap, iterations) %>%
  count()
count(all_params)

saveRDS(all_params, file = "research/step-function-simulations/simulation_parameters.rds")

all_params %>%
  filter(bootstrap == "none") %>%
  select(row) %>%
  distinct() %>%
  write_csv("research/step-function-simulations/batches-to-run.csv", col_names = FALSE)

all_params %>%
  filter(bootstrap == "exponential") %>%
  select(row) %>%
  distinct() %>%
  write_csv("research/step-function-simulations/batches-to-run.csv", col_names = FALSE)
