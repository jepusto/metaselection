library(tidyverse)
library(tictoc)
library(simhelpers)
library(future)
library(furrr)

rows <- 
  read_csv("research/step-function-simulations/batch-results/batches-to-run.csv", col_names = "row") %>%
  mutate(
    i = row_number()
  )

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
  anti_join(res_list, by = "row") %>%
  left_join(rows, by = "row")

nrow(outstanding_conditions)
outstanding_conditions %>%
  count(bootstrap)

#-------------------------------------------------------------------------------
# compile results from conditions with no bootstraps

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

tic()
no_bootstraps_res <- 
  res_list %>%
  filter(row < 1e5) %>%
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
    mean_smd:omega, steps, iterations,
    model, estimator, param, 
    K_absolute:rmse_mcse, 
    K_coverage:width_mcse
  )

res %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, weights, m, n_multiplier, omega, steps) %>%
  summarize(n_res = n(), .groups = "drop") %>%
  count(n_res)

res %>%
  filter(
    mean_smd == 0, tau == 0.05, cor_mu == 0.4, omega == 0,
    weights == 0.10, m == 60, n_multiplier == 1,
  ) %>%
  select(priors, iterations, model:width_mcse) %>%
  View()
write_rds(res, file = "research/step-function-simulations/sim-step-function-results-no-bootstraps.rds", compress = "gz", compression = 9L)


#-------------------------------------------------------------------------------
# Count conditions with incomplete bootstraps

bootstrap_files <-
  params %>%
  inner_join(res_list, by = "row") %>%
  select(-batch, -row, -seed) %>%
  filter(bootstrap != "none") %>%
  nest(iterations = iterations, files = file)


# count_booties <- function(file_list) {
#   dat <- map_dfr(file_list$file, .f = readRDS)
#   
#   bind_rows(dat$res) %>%
#     filter(param == "beta", model == "3PSM") %>%
#     mutate(R = map_int(boot_CIs, \(x) if (is.null(x)) 0L else max(x$bootstraps))) %>%
#     group_by(estimator) %>%
#     count(R)
# }
# 
# 
# file_list <-
#   bootstrap_files %>%
#   filter(mean_smd == 0.0, tau == 0.05, cor_mu == 0.8, m == 60, omega == 0, weights == 0.05, n_multiplier == 1/3, bootstrap == "two-stage") %>%
#   pull(files)
# count_booties(file_list[[1]])
# 
# plan(multisession, workers = 10L)
# 
# tic()
# bootstraps_counts <- 
#   bootstrap_files %>%
#   mutate(
#     iterations = future_map_int(iterations, ~ sum(.x$iterations)),
#     counts = future_map(files, count_booties, .progress = TRUE)
#   )
# toc()
# 
# plan(sequential)
# 
# bootstrap_reps <- 
#   bootstraps_counts %>%
#   select(priors, bootstrap, iterations, counts) %>%
#   unnest(counts) %>% 
#   group_by(priors, bootstrap, estimator, R) %>%
#   summarize(
#     reps = sum(n),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     status = case_when(R == 0 ~ "none", R == 299 ~ "all", TRUE ~ "some")
#   )
# 
# bootstrap_reps %>% 
#   group_by(priors, bootstrap, estimator, status) %>%
#   summarize(
#     reps = sum(reps)
#   ) 

#-------------------------------------------------------------------------------
# compile results from conditions with bootstraps


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

summarize_bootstraps(file_list = bootstrap_files$files[[1]])

# partial_files <- 
#   bootstraps_counts %>%
#   select(files, counts) %>%
#   mutate(
#     partial = map_lgl(counts, ~ any(.x$R < 299 & .x$R > 0))
#   ) %>%
#   filter(partial) %>%
#   pull(files)
# 
# tic()
# summary_dat <- summarize_bootstraps(file_list = partial_files[[1]])
# toc()

plan(multisession, workers = 10L)

tic()
bootstrap_res <- 
  bootstrap_files %>%
  mutate(
    summary_file = future_map_chr(files, .f = summarize_bootstraps, .progress = TRUE)
  )
toc()

plan(sequential)


tic()
bootstrap_res <- 
  bootstrap_res %>%
  mutate(
    iterations = map_int(iterations, ~ sum(.x$iterations)),
    res = map(summary_file, .f = read_rds, .progress = TRUE)
  ) %>%
  select(-files, -summary_file) %>%
  unnest(res)
toc()

write_rds(bootstrap_res, file = "research/step-function-simulations/sim-step-function-bootstrap-performance-results.rds", compress = "gz", compression = 9L)

#-------------------------------------------------------------------------------
# Compile computation time data

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
