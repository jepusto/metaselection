args <- commandArgs(trailingOnly = TRUE)

#-------------------------------------------------------------------------------
# Parse command line argument

row_to_run <- 
  args |>
  paste(collapse = " ") |>
  stringr::str_extract("batch [0-9]+") |>
  stringr::str_sub(7, -1) |>
  as.integer()

#-------------------------------------------------------------------------------
# Source packages and functions
library(dplyr)
library(tidyr)
library(purrr)
library(tictoc)

source("research/step-function-simulations/1_estimation.R")
source("research/step-function-simulations/2_performance_criteria.R")
source("research/step-function-simulations/3_simulation_driver.R")

#-------------------------------------------------------------------------------
# Load experimental design parameters

load("research/step-function-simulations/wwc_es.RData")

#-------------------------------------------------------------------------------
# Experimental Design

# Express the simulation parameters as vectors/lists

design_factors <- list(
  mean_smd = c(0.2), # average effect size
  tau = c(0.05, 0.45), # between study heterogeneity
  het_ratio = c(0.5),
  cor_mu = c(0.8), # average correlation between outcomes         
  cor_sd = c(0.05), # sd correlation between outcomes
  weights = c(0.05, 0.2), # weights
  m = c(15),	# number of studies in each meta analysis
  n_multiplier = c(1),
  batch = 1:200
)

params <- do.call(expand_grid, design_factors)

all_params <- 
  params |>
  mutate(
    omega = tau * sqrt(het_ratio),
    steps = 0.025,
    bootstrap = "two-stage",
    R = 1999L, 
    iterations = 20L,
    summarize_performance = FALSE,
    comparison_methods = "none",
    seed = 20250801L + 1:dplyr::n(),
    row = row_number()
  ) %>%
  select(-het_ratio)

saveRDS(all_params, file = "research/step-function-simulations/big_B_bootstrap_parameters.rds")

#-------------------------------------------------------------------------------
# run simulations for specified batch
res <- subset(all_params, row == row_to_run)
res$batch <- NULL
res$row <- NULL

tic()

res$res <- pmap(res, .f = run_sim)

tm <- toc(quiet = TRUE)


#-------------------------------------------------------------------------------
# Save results and details

res$run_date <- date()
res$time <- tm$toc - tm$tic

saveRDS(res, file = paste0("big_B_simulation_results_batch", row_to_run, ".rds"))