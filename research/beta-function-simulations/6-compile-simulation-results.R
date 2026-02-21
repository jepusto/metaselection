library(tidyverse)
library(tictoc)
library(simhelpers)
library(future)
library(furrr)


params <- readRDS("research/beta-function-simulations/simulation_parameters.rds")

res_list <- tibble(
  file = list.files("research/beta-function-simulations/batch-results", pattern = "simulation_results_batch", full.names = TRUE)
) %>%
  mutate(
    row = str_extract(file, "batch[0-9]+.rds") |> str_sub(6,-5) |> as.integer()
  )
nrow(res_list)

outstanding_conditions <-
  params %>%
  anti_join(res_list, by = "row")

outstanding_conditions %>% 
  mutate(
    batch = floor((batch - 1)/ 20) + 1
  ) %>%
  count(batch)

# outstanding_conditions %>%
#   select(row) %>%
#   write_csv(file = "research/beta-function-simulations/batches-to-run.csv", col_names = FALSE)


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
    mean_smd:m, omega, delta_2, iterations,
    model, estimator, param, 
    K_absolute:rmse_mcse, 
    K_coverage:width_mcse
  )

res %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, delta_1, m, omega) %>%
  summarize(n_res = n(), .groups = "drop") %>%
  count(n_res)

res %>%
  filter(
    mean_smd == 0, tau == 0.05, cor_mu == 0.4, omega == 0, delta_1 == 0.10, m == 60,
  ) %>%
  select(iterations, model:width_mcse)

write_rds(res, file = "research/beta-function-simulations/sim-beta-function-point-estimator-results.rds", compress = "gz", compression = 9L)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# compile results from conditions with bootstraps ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

source("research/beta-function-simulations/2-performance-criteria.R")

bootstrap_files <- 
  res_list %>%
  left_join(params, by = "row") %>%
  filter(bootstrap != "none") %>%
  select(-batch, -seed, -row) %>%
  nest(iterations = iterations, files = file) %>%
  mutate(nfiles = map_int(files, nrow))

bootstrap_files %>%
  count(nfiles)


summarize_bootstraps <- function(file_list) {
  
  batch_file_name <- paste0(
    "research/beta-function-simulations/batch-results/simulation_results_bootstrap_batch",
    str_match(file_list$file[[1]], "_batch(.+).rds")[,2],
    ".rds"
  )
  
  if (file.exists(batch_file_name)) return(batch_file_name)
  
  dat <- map_dfr(file_list$file, .f = readRDS)
  time <- sum(dat$time)
  run_date <- min(dat$run_date)
  
  true_params <- data.frame(
    param = c("beta", "gamma", "zeta1", "zeta2"),
    true_param = c(unique(dat$mean_smd), log(unique(dat$tau)^2 + unique(dat$omega)^2), log(unique(dat$delta_1)), log(unique(dat$delta_2)))
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

# debug(summarize_bootstraps)
# tic()
# summarize_bootstraps(file_list = bootstrap_files$files[[1]])
# toc()

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

write_rds(bootstrap_res, file = "research/beta-function-simulations/sim-beta-function-bootstrap-performance-results.rds", compress = "gz", compression = 9L)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Compile computation time data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# timing

timings <- 
  bind_rows(no_bootstraps_res, bootstrap_res) %>%
  select(-res, -summarize_performance) %>%
  mutate(time_hrs = time / 60^2)

write_rds(timings, file = "research/beta-function-simulations/sim-beta-function-timings.rds", compress = "gz", compression = 9L)

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
