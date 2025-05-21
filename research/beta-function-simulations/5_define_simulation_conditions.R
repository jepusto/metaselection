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
  delta_1 = c(0.01, 0.2, 0.5, 0.8, 1), # selection parameters
  m = c(90, 60, 30, 15),	# number of studies in each meta analysis					
  batch = 1:50
)

params <- do.call(expand_grid, design_factors)

all_params <- 
  params |>
  mutate(
    omega = tau * sqrt(het_ratio),
    delta_2 = if_else(delta_1 == 1, 1, 0.9),
    bootstrap = if_else(
      (m <= 60) & (cor_mu == 0.8) & (tau %in% c(0.05, 0.15, 0.45)) & delta_1 %in% c(0.2, 0.5, 1),
      "two-stage",
      "none"
    ),
    iterations = case_match(bootstrap, "two-stage" ~ 40L, "none" ~ 2000L),
    summarize_performance = bootstrap == "none",
    seed = 20250520L + 1:dplyr::n(),
    row = row_number()
  ) %>%
  filter(bootstrap == "two-stage" | batch == 1L) %>%
  select(-het_ratio)

saveRDS(all_params, file = "research/beta-function-simulations/simulation_parameters.rds")

all_params %>% 
  slice_sample(n = 20, by = bootstrap) %>%
  select(row) %>%
  write_csv(file = "research/beta-function-simulations/batches-to-run.csv", col_names = FALSE)
