#-------------------------------------------------------------------------------
# Source packages and functions
library(dplyr)
library(tidyr)

#-------------------------------------------------------------------------------
# Experimental Design

# Express the simulation parameters as vectors/lists

design_factors <- list(
  mean_smd = c(0.0, 0.1, 0.2, 0.4, 0.8), # average effect size
  tau = c(0.05, 0.15, 0.30, 0.45, 0.60), # between study heterogeneity
  het_ratio = c(0, 0.5),
  cor_mu = c(0.4, 0.8), # average correlation between outcomes         
  cor_sd = c(0.05), # sd correlation between outcomes
  weights = c(0.02, 0.05, 0.1, 0.2, 0.5, 1), # weights
  m = c(90, 60, 30, 15),	# number of studies in each meta analysis					
  batch = 1
)

params <- do.call(expand_grid, design_factors)

all_params <- 
  params |>
  mutate(
    omega = tau * sqrt(het_ratio),
    steps = 0.025,
    seed = 20240507L + 1:dplyr::n()
  ) %>%
  select(-het_ratio)

saveRDS(all_params, file = "simulation-step-function-revised/simulation_parameters.rds")
