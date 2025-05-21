args <- commandArgs(trailingOnly = TRUE)

#-------------------------------------------------------------------------------
# Parse command line argument

batch <- 
  args |>
  paste(collapse = " ") |>
  stringr::str_extract("batch [0-9]+") |>
  stringr::str_sub(7, -1) |>
  as.integer()

iterations <- 
  args |> 
  paste(collapse = " ") |>
  stringr::str_extract("iterations [0-9]+") |>
  stringr::str_sub(12, -1) |>
  as.integer()

#-------------------------------------------------------------------------------
# Source packages and functions
library(dplyr)
library(tidyr)
library(purrr)
library(tictoc)


source("R/data-generation.R")
source("R/step-utilities.R")
source("R/step-function.R")
source("R/hybrid-estimating.R")
source("R/selection_model.R")
source("simulation-step-function-revised/2_estimation.R")
source("simulation-step-function-revised/3_performance_criteria.R")
source("simulation-step-function-revised/4_simulation_driver.R")

#-------------------------------------------------------------------------------
# Load experimental design parameters

load("simulation-step-function-revised/wwc_es.RData")

all_params <- readRDS("simulation-step-function-revised/simulation_parameters.rds")

#-------------------------------------------------------------------------------
# run simulations for specified batch
res <- all_params[batch,]
res$iterations <- iterations
res$batch <- NULL

tic()

res$res <- pmap(
  res, .f = run_sim, 
  stepfun_methods = c("step-MLE", "step-hybrid-Broyden", "step-hybrid-Newton","step-hybrid-rootSolve"),
  ML_optimizer_control = list(all.methods = TRUE),
  comparison_methods = NULL,
  summarize_performance = TRUE,
  filter_n_study = 2L
)

tm <- toc(quiet = TRUE)


#-------------------------------------------------------------------------------
# Save results and details

res$run_date <- date()
res$time <- tm$toc - tm$tic

saveRDS(res, file = paste0("simulation_results_batch", batch, ".rds"))