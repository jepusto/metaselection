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
# Source packages, functions, and design data
library(simhelpers)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

source("research/step-function-simulations/2_performance_criteria.R")

params <- readRDS("research/step-function-simulations/simulation_parameters.rds")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# compile results from conditions with bootstraps ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

res_list <- 
  tibble(
    # file = list.files("research/step-function-simulations/batch-results", pattern = "simulation_results_batch", full.names = TRUE)
    file = list.files(pattern = "simulation_results_batch")
  ) %>%
  mutate(
    row = str_extract(file, "batch[0-9]+.rds") |> str_sub(6,-5) |> as.integer()
  )

bootstrap_files <-
  params %>%
  left_join(res_list, by = "row") %>%
  filter(bootstrap != "none") %>%
  select(-seed) %>%
  nest(batches = batch, iterations = iterations, rows = row, files = file) %>%
  mutate(
    R_max = map_dbl(R, max),
    nbatches = map_dbl(batches, nrow),
    row = row_number()
  )

file_list <- bootstrap_files$files[[row_to_run]]

cat(paste(file_list$file,"\n"))

batch_file_name <- paste0(
  "simulation_results_bootstrap_batch",
  str_match(file_list$file[[1]], "_batch(.+).rds")[,2],
  ".rds"
)

dat <- map_dfr(file_list$file, .f = readRDS)
time <- sum(dat$time)
run_date <- min(dat$run_date)

true_params <- data.frame(
  param = c("beta", "gamma", "zeta1"),
  true_param = c(unique(dat$mean_smd), log(unique(dat$tau)^2 + unique(dat$omega)^2), log(unique(dat$weight)))
)

results <-
  bind_rows(dat$res, .id = "file") %>%
  mutate(rep = as.character(as.integer(file) * 1000 + as.integer(rep))) %>%
  left_join(true_params, by = "param")

summary_res <-
  results %>%
  calc_performance() %>%
  nest(res = everything()) %>%
  mutate(
    run_date = run_date,
    time = time
  )

saveRDS(summary_res, file = batch_file_name)
