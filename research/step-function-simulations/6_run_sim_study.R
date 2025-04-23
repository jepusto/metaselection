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

source("research/step-function-simulations/2_estimation.R")
source("research/step-function-simulations/3_performance_criteria.R")
source("research/step-function-simulations/4_simulation_driver.R")

#-------------------------------------------------------------------------------
# Load experimental design parameters

load("research/step-function-simulations/wwc_es.RData")

all_params <- readRDS("research/step-function-simulations/simulation_parameters.rds")

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

saveRDS(res, file = paste0("simulation_results_batch", row_to_run, ".rds"))