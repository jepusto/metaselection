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
  priors = c("Flat","Mild","Weak")
)

params <- 
  expand_grid(!!!design_factors) |>
  mutate(
    omega = tau * sqrt(het_ratio),
    steps = 0.025,
    bootstrap = "none",
    comparison_methods = if_else(priors == "Weak", "All", "None"),
    iterations = 2000L,
    summarize_performance = TRUE,
    row = rep(1:(dplyr::n()/3L), each = 3L),
    seed = 20250918L + row
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
  batch = 1:10L,
  priors = "Weak",
  bootstrap = c("multinomial","two-stage","exponential")
)

bootstrap_params <- 
  expand_grid(!!!bootstrap_factors) %>%
  mutate(
    omega = tau * sqrt(het_ratio),
    steps = 0.025,
    comparison_methods = "None",
    iterations = 200L,
    summarize_performance = FALSE,
    row = 1e5 + 1:(dplyr::n()),
    seed = 20250918L + rep(1:(dplyr::n() / 3L), each = 3L)
  ) %>%
  select(-het_ratio)

all_params <- bind_rows(params, bootstrap_params)

all_params %>%
  count(bootstrap)

saveRDS(all_params, file = "research/step-function-simulations/simulation_parameters.rds")

all_params %>%
  filter(bootstrap == "none", priors == "Weak") %>%
  select(row) %>%
  distinct() %>%
  write_csv("research/step-function-simulations/batches-to-run.csv", col_names = FALSE)

all_params %>%
  filter(bootstrap == "multinomial") %>%
  select(row) %>%
  distinct() %>%
  write_csv("research/step-function-simulations/batches-to-run.csv", col_names = FALSE)
