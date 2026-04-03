library(tidyverse)
library(tictoc)
library(simhelpers)
library(future)
library(furrr)

params <- readRDS("research/step-function-simulations/simulation_parameters.rds")

res_list <- tibble(
  file = list.files("research/step-function-simulations/batch-results", pattern = "simulation_results_batch", full.names = TRUE)
) %>%
  mutate(
    row = str_extract(file, "batch[0-9]+.rds") |> str_sub(6,-5) |> as.integer()
  )
nrow(res_list)


outstanding_conditions <-
  params %>%
  anti_join(res_list, by = "row")

nrow(outstanding_conditions)
outstanding_conditions %>%
  count(bootstrap, psi)

# revise_batch_file <- function(file) {
#   x <- readRDS(file)
#   if ("weights" %in% names(x)) {
#     x <- rename(x, weight = weights)
#   }
#   if (!("psi" %in% names(x))) {
#     x$psi <- 0
#     x <- relocate(x, psi, .after = n_multiplier)
#   }
#   
#   saveRDS(x, file = file)
#   return(NULL)
# }
# 
# walk(res_list$file, revise_batch_file)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Compile results from conditions with no bootstraps ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

plan(multisession, workers = 10)

tic()
no_bootstraps_res <- 
  params %>%
  left_join(res_list, by = "row") %>%
  filter(bootstrap == "none", !is.na(file)) %>%
  select(-priors, -comparison_methods) %>%
  distinct() %>%
  pull(file) %>%
  future_map_dfr(.f = readRDS) %>%
  select(-seed)
toc()

plan(sequential)

nrow(no_bootstraps_res)

# Arrange point estimation results for further analysis

res <- 
  no_bootstraps_res %>%
  select(-run_date) %>%
  unnest(res) %>%
  select(
    mean_smd:psi, priors, bootstrap, omega, 
    steps, iterations,
    model, estimator, param, 
    K_absolute:rmse_mcse, 
    K_coverage:width_mcse
  )

res %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, weight, psi, m, n_multiplier, omega, steps) %>%
  summarize(n_res = n(), .groups = "drop") %>%
  count(n_res)

res %>%
  filter(
    mean_smd == 0, tau == 0.05, cor_mu == 0.4, omega == 0,
    weight == 0.10, psi == 0, m == 60, n_multiplier == 1,
  ) %>%
  select(priors, iterations, model:width_mcse) %>%
  View()

write_rds(res, file = "research/step-function-simulations/sim-step-function-results-no-bootstraps.rds", compress = "gz", compression = 9L)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Count conditions with incomplete bootstraps ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

bootstrap_files <-
  params %>%
  inner_join(res_list, by = "row") %>%
  select(-batch, -row, -seed) %>%
  filter(bootstrap != "none") %>%
  nest(iterations = iterations, files = file) %>%
  mutate(R_max = map_dbl(R, max))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# compile results from conditions with bootstraps ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


source("research/step-function-simulations/2_performance_criteria.R")

summarize_bootstraps <- function(file_list) {
  
  batch_file_name <- paste0(
    "research/step-function-simulations/batch-results/simulation_results_bootstrap_batch",
    str_match(file_list$file[[1]], "_batch(.+).rds")[,2],
    ".rds"
  )
  
  if (file.exists(batch_file_name)) return(batch_file_name)
  
  dat <- map_dfr(file_list$file, .f = readRDS)
  time <- sum(dat$time)
  run_date <- min(dat$run_date)
  
  true_params <- data.frame(
    param = c("beta", "gamma", "zeta1"),
    true_param = c(unique(dat$mean_smd), log(unique(dat$tau)^2 + unique(dat$omega)^2), log(unique(dat$weights)))
  )
  
  results <-
    bind_rows(dat$res) %>%
    left_join(true_params, by = "param")
  
  summary_res <-
    results %>%
    calc_performance() %>%
    nest(res = everything()) %>%
    mutate(
      run_date = run_date,
      time = time
    )
  
  write_rds(summary_res, file = batch_file_name, compress = "gz", compression = 9L)
  
  return(batch_file_name)
}

# debug(calc_performance)
# summarize_bootstraps(file_list = bootstrap_files$files[[1]])
# summarize_bootstraps(file_list = bootstrap_files$files[[1296]])

plan(multisession, workers = 10L)

tic()
bootstrap_res <- 
  bootstrap_files %>%
  mutate(
    summary_file = future_map_chr(files, .f = summarize_bootstraps, .progress = TRUE)
  )
toc()

tic()
bootstrap_res <- 
  bootstrap_res %>%
  mutate(
    iterations = future_map_int(iterations, ~ sum(.x$iterations)),
    res =  future_map(summary_file, .f = read_rds, .progress = TRUE)
  ) %>%
  select(-files, -summary_file) %>%
  unnest(res)
toc()

plan(sequential)

write_rds(bootstrap_res, file = "research/step-function-simulations/sim-step-function-bootstrap-performance-results.rds", compress = "gz", compression = 9L)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Compile computation time data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# timing

timings <- 
  bind_rows(no_bootstraps_res, bootstrap_res) %>%
  select(-res, -summarize_performance) %>%
  mutate(time_hrs = time / 60^2)

write_rds(timings, file = "research/step-function-simulations/sim-step-function-timings.rds", compress = "gz", compression = 9L)

timings %>%
  mutate(bootstrap = "All") %>%
  bind_rows(timings) %>%
  group_by(bootstrap) %>%
  summarize(
    across(time_hrs, .fns = c(min = min, median = median, max = max, mean = mean, total = sum))
  ) %>%
  mutate(
    time_yrs_total = time_hrs_total / 24 / 365.25
  )
