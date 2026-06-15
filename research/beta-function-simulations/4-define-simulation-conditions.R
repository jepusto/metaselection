#-------------------------------------------------------------------------------
# Source packages and functions
library(dplyr)
library(tidyr)
library(readr)

#-------------------------------------------------------------------------------
# Experimental Design

# Express the simulation parameters as vectors/lists

design_factors <- list(
  mean_smd = c(0.0, 0.2, 0.4, 0.8), # average effect size
  tau = c(0.05, 0.15, 0.30, 0.45), # between study heterogeneity
  het_ratio = c(0, 0.5),
  cor_mu = c(0.4, 0.8), # average correlation between outcomes         
  cor_sd = c(0.05), # sd correlation between outcomes
  delta_1 = c(0.02, 0.1, 0.2, 0.5, 1), # selection parameters
  m = c(90, 60, 30, 15),	# number of studies in each meta analysis
  step_models = list(c("3PSM","4PSM"))
)

lengths(design_factors)
lengths(design_factors) |> prod()

design_conditions <- 
  expand_grid(!!!design_factors) %>%
  mutate(
    bootstrap = "none",
    batch = 1,
    iterations = 2000L,
    comparison_methods = "All",
    summarize_performance = TRUE
  )


bootstrap_factors <- list(
  mean_smd = c(0.0, 0.2, 0.4, 0.8), # average effect size
  tau = c(0.05, 0.30), # between study heterogeneity
  het_ratio = c(0, 0.5),
  cor_mu = 0.8, # average correlation between outcomes         
  cor_sd = 0.05, # sd correlation between outcomes
  delta_1 = c(0.2, 0.5, 1), # weights
  m = c(60, 30, 15),# number of studies in each meta analysis
  step_models = list(c("3PSM","4PSM")),
  batch = 1:80
)

lengths(bootstrap_factors)
lengths(bootstrap_factors[1:7]) |> prod()

bootstrap_conditions <- 
  expand_grid(!!!bootstrap_factors) %>%
  mutate(
    bootstrap = "two-stage",
    iterations = 25L,
    comparison_methods = "None",
    summarize_performance = FALSE
  )

all_params <- 
  bind_rows(design_conditions, bootstrap_conditions) |>
  mutate(
    omega = tau * sqrt(het_ratio),
    delta_2 = if_else(delta_1 == 1, 1, 0.9),
    priors = "Weak",
    seed = 20250520L + 1:dplyr::n(),
    row = row_number()
  ) %>%
  select(-het_ratio)

saveRDS(all_params, file = "research/beta-function-simulations/simulation_parameters.rds")

all_params %>% 
  filter(bootstrap == "two-stage", batch >= 61L) %>%
  select(row) %>%
  write_csv(file = "research/beta-function-simulations/batches-to-run.csv", col_names = FALSE)
