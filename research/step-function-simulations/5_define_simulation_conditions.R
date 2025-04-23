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
  batch = 1:20
)

params <- do.call(expand_grid, design_factors)

all_params <- 
  params |>
  mutate(
    omega = tau * sqrt(het_ratio),
    steps = 0.025,
    bootstrap = if_else(
      (m <= 60) & (cor_mu == 0.8) & (mean_smd %in% c(0.0, 0.2, 0.4, 0.8)) & (tau %in% c(0.05, 0.45)) & weights %in% c(0.05, 0.2, 1),
      "exponential",
      "none"
    ),
    iterations = case_match(bootstrap, "exponential" ~ 100L, "none" ~ 2000L),
    summarize_performance = bootstrap == "none",
    seed = 20250218L + 1:dplyr::n(),
    row = row_number()
  ) %>%
  filter(bootstrap == "exponential" | batch == 1L) %>%
  select(-het_ratio)

all_params <- 
  all_params %>%
  filter(bootstrap == "exponential") %>%
  mutate(bootstrap = "multinomial") %>%
  bind_rows(all_params)

all_params <- 
  all_params %>%
  filter(bootstrap == "exponential") %>%
  mutate(bootstrap = "two-stage") %>%
  bind_rows(all_params, .) %>%
  mutate(
    row = row_number(),
    comparison_methods = if_else(bootstrap %in% c("exponential","two-stage"),"None","All")
  )

all_params %>%
  count(bootstrap)
nrow(all_params)

saveRDS(all_params, file = "research/step-function-simulations/simulation_parameters.rds")

all_params %>%
  group_by(bootstrap) %>%
  sample_n(size = 1000L) %>%
  ungroup() %>%
  select(row) %>%
  write_csv("research/step-function-simulations/batches-to-run.csv", col_names = FALSE)

all_params %>%
  filter(bootstrap == "two-stage") %>%
  select(row) %>%
  write_csv("research/step-function-simulations/batches-to-run.csv", col_names = FALSE)
