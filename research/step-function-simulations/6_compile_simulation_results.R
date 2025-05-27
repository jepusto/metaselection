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

outstanding_conditions %>%
  select(row, i, everything()) %>%
  filter(!is.na(i))

nrow(outstanding_conditions)

outstanding_conditions %>%
  select(row) %>%
  write_csv(file = "research/step-function-simulations/batches-to-run.csv", col_names = FALSE)


#-------------------------------------------------------------------------------
# compile results from conditions with no bootstraps

tic()
no_bootstraps_res <- 
  params %>%
  left_join(res_list, by = "row") %>%
  filter(bootstrap == "none", !is.na(file)) %>%
  pull(file) %>%
  future_map_dfr(.f = readRDS) %>%
  select(-seed)
toc()

params %>%
  count(bootstrap)
nrow(no_bootstraps_res)

#-------------------------------------------------------------------------------
# compile results from conditions with bootstraps

count_booties <- function(file) {
  x <- readRDS(file)
  x$res[[1]] %>%
    filter(param == "beta", model == "3PSM") %>%
    mutate(R = map_int(boot_CIs, \(x) if (is.null(x)) 0L else max(x$bootstraps))) %>%
    group_by(estimator) %>%
    count(R) %>%
    mutate(
      file = file,
      run_date = x$run_date
    )
}
  
plan(multisession, workers = 10L)

tic()
bootstraps_res <- 
  params %>%
  inner_join(res_list, by = "row") %>%
  filter(bootstrap != "none") %>%
  pull(file) %>%
  future_map_dfr(.f = count_booties, .progress = TRUE)
toc()

plan(sequential)

bootstrap_reps <- 
  bootstraps_res %>%
  group_by(run_date, file) %>%
  summarize(
    complete = all(n == 100),
    almost = all(R %in% c(0,399)),
    has_zeros = any(R == 0),
    .groups = "drop"
  ) %>%
  mutate(
    date = str_sub(run_date, 1, 10)
  )

bootstrap_reps %>%
  summarize(
    n = n(),
    complete = sum(complete),
    almost = sum(almost & !complete),
    has_zeros = sum(has_zeros)
  )


bootstrap_reps %>%
  group_by(complete, has_zeros, date) %>%
  count()

bootstrap_reps %>%
  filter(!complete, date %in% c("Sun Feb 23")) %>%
  left_join(bootstraps_res, by = join_by(run_date, file)) %>%
  group_by(estimator, R) %>%
  summarize(
    files = n(),
    reps = sum(n)
  )

#-------------------------------------------------------------------------------
# compile results from conditions with bootstraps

bootstrap_files <- 
  params %>%
  inner_join(res_list, by = "row") %>%
  select(-batch, -row, -seed) %>%
  filter(bootstrap != "none") %>%
  nest(iterations = iterations, files = file)

file_list <- 
  bootstrap_files %>%
  filter(mean_smd == 0.0, tau == 0.05, cor_mu == 0.8, m == 60, omega == 0, weights == 0.05, n_multiplier == 1/3, bootstrap == "multinomial") %>%
  pull(files)

dat <- 
  map_dfr(file_list[[1]]$file, .f = readRDS) %>%
  unnest(res)

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
    calc_performance(winz = 2.5, B_target = 1999)
  
  summary_res %>%
    nest(res = everything()) %>%
    mutate(
      run_date = run_date,
      time = time
    )
  
}

summary_dat <- summarize_bootstraps(file_list = file_list[[1]])

plan(multisession, workers = 10L)

tic()
bootstrap_res <- 
  bootstrap_files %>%
  mutate(
    iterations = future_map_int(iterations, ~ sum(.x$iterations)),
    all_res = future_map(files, .f = summarize_bootstraps, .progress = TRUE)
  ) %>%
  unnest(all_res) %>%
  select(-files, -summarize_performance)
toc()

plan(sequential)

write_rds(bootstrap_res, file = "research/step-function-simulations/sim-step-function-bootstrap-performance-results.rds", compress = "gz", compression = 9L)
bootstrap_res <- readRDS("research/step-function-simulations/sim-step-function-bootstrap-performance-results.rds")

convergence_rates <- 
  bootstrap_res %>%
  select(-run_date, -time, -comparison_methods) %>%
  unnest(res) %>%
  filter(param == "beta", model == "3PSM") %>%
  mutate(
    cnvrg_rate = K_absolute / iterations
  )

convergence_rates %>%
  arrange(cnvrg_rate) %>%
  select(cnvrg_rate, everything())
ggplot(convergence_rates) + 
  aes(cnvrg_rate, fill = estimator) + 
  geom_density(alpha = 0.5)

#-------------------------------------------------------------------------------
# organize results for further analysis

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


# point estimation results

bootstrap_point_estimator <- 
  bootstrap_res %>%
  select(-run_date) %>%
  filter(bootstrap == "multinomial") %>%
  unnest(res) %>%
  select(
    time,  
    mean_smd:omega, steps, bootstrap_condition = bootstrap,
    model, estimator, param, 
    K_absolute:rmse_mcse, 
    est_winsor_pct:var_winsor_pct_mcse
  ) %>%
  mutate(
    bootstrap_condition = "bootstrap"
  )


res_point_estimator <- 
  no_bootstraps_res %>%
  select(-run_date) %>%
  unnest(res) %>%
  select(
    time, 
    mean_smd:omega, steps, bootstrap_condition = bootstrap,
    model, estimator, param, 
    K_absolute:rmse_mcse, 
    est_winsor_pct, est_winsor_pct_mcse
  ) %>%
  bind_rows(bootstrap_point_estimator)

res_point_estimator %>%
  count(bootstrap_condition)
res_point_estimator %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, weights, m, n_multiplier, omega, steps) %>%
  summarize(n_res = n(), .groups = "drop") %>%
  count(n_res)


write_rds(res_point_estimator, file = "research/step-function-simulations/sim-step-function-point-estimator-results.rds", compress = "gz", compression = 9L)

# variance estimation and confidence interval results

nobootstrap_conf_int <- 
  no_bootstraps_res %>%
  select(-run_date, -time) %>%
  unnest(res) %>%
  select(mean_smd:steps, bootstrap_condition = bootstrap, model:param, 
         K_coverage:width_mcse,
         K_relvar:rel_rmse_var_mcse, 
         var_winsor_pct, var_winsor_pct_mcse) %>%
  mutate(
    CI_type = "large-sample"
  ) 

bootstrap_large_sample <- 
  bootstrap_res %>%
  select(-run_date, -time) %>%
  unnest(res) %>%
  select(mean_smd:steps, bootstrap_condition = bootstrap, model:param, 
         K_coverage:width_mcse, 
         K_relvar:rel_rmse_var_mcse, 
         var_winsor_pct, var_winsor_pct_mcse) %>%
  mutate(
    CI_type = "large-sample",
    bootstrap_condition = "bootstrap"
  )

bootstrap_conf_int <- 
  bootstrap_res %>%
  select(-run_date, -time) %>%
  unnest(res) %>%
  select(mean_smd:steps, bootstrap_condition = bootstrap, bootstrap_type = bootstrap, model:param, 
         bootstraps, boot_coverage, boot_coverage_mcse, boot_width, boot_width_mcse) %>%
  unnest(
    c(bootstraps, boot_coverage, boot_coverage_mcse, boot_width, boot_width_mcse),
    names_sep = "-"
  ) %>%
  pivot_longer(
    starts_with("boot_"),
    names_to = c(".value", "CI_type"),
    names_pattern = "(.+)-(.+)"
  ) %>%
  rename_with(~ str_remove(.x, "^boot_")) %>%
  mutate(
    bootstrap_condition = "bootstrap"
  )


res_conf_int <- bind_rows(
  nobootstrap_conf_int,
  bootstrap_conf_int,
  bootstrap_large_sample
)

res_conf_int %>%
  count(bootstrap_condition)

res_conf_int %>%
  group_by(model, estimator, param, bootstrap_type) %>% count()

write_rds(res_conf_int, file = "research/step-function-simulations/sim-step-function-confidence-interval-results.rds", compress = "gz", compression = 9L)
