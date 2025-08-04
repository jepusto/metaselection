library(tidyverse)
library(tictoc)
library(simhelpers)
library(future)
library(furrr)


params <- readRDS("research/step-function-simulations/big_B_bootstrap_parameters.rds")

res_list <- tibble(
  file = list.files("research/step-function-simulations/big-B-results", pattern = "simulation_results_batch", full.names = TRUE)
) %>%
  mutate(
    row = str_extract(file, "batch[0-9]+.rds") |> str_sub(6,-5) |> as.integer()
  )
nrow(res_list)

#-------------------------------------------------------------------------------
# compile results from conditions with bootstraps

bootstrap_files <- 
  params %>%
  inner_join(res_list, by = "row") %>%
  select(-batch, -row, -seed) %>%
  nest(iterations = iterations, files = file)

file_list <- 
  bootstrap_files %>%
  filter(tau == 0.05, weights == 0.05) %>%
  pull(files)

source("research/step-function-simulations/2_performance_criteria.R")

summarize_bootstraps <- function(file_list) {
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
    group_by(estimator, param) %>%
    summarize(
      bootstraps = mean(bootstraps, na.rm = TRUE),
      basic = calc_coverage(lower_bound = basic_lower, upper_bound = basic_upper,true_param = true_param, winz = 2.5),
      student = calc_coverage(lower_bound = student_lower, upper_bound = student_upper,true_param = true_param, winz = 2.5),
      percentile = calc_coverage(lower_bound = percentile_lower, upper = percentile_upper,true_param = true_param, winz = 2.5),
      biascorrected = calc_coverage(lower_bound = biascorrected_lower, upper_bound = biascorrected_upper,true_param = true_param, winz = 2.5),
      BCa = calc_coverage(lower_bound = BCa_lower, upper_bound = BCa_upper,true_param = true_param, winz = 2.5),
      .groups = "drop"
    ) %>%
    pivot_longer(
      c(basic,student,percentile,biascorrected,BCa), 
      names_to = "CI_type",
      values_to = "coverage"
    ) %>%
    unnest(coverage)
  
  summary_res %>%
    nest(res = everything()) %>%
    mutate(
      run_date = run_date,
      time = time
    )
  
}

summary_dat <- summarize_bootstraps(file_list = file_list[[1]])

tic()
bootstrap_res <- 
  bootstrap_files %>%
  mutate(
    iterations = map_int(iterations, ~ sum(.x$iterations)),
    all_res = map(files, .f = summarize_bootstraps, .progress = TRUE)
  ) %>%
  unnest(all_res) %>%
  select(-files, -summarize_performance)
toc()


write_rds(bootstrap_res, file = "research/step-function-simulations/sim-step-function-big-B-bootstrap-performance-results.rds", compress = "gz", compression = 9L)
