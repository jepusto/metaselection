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

outstanding_conditions 

outstanding_conditions %>%
  select(row) %>%
  write_csv(file = "research/beta-function-simulations/batches-to-run.csv", col_names = FALSE)

#-------------------------------------------------------------------------------
# compile results from conditions with no bootstraps

tic()
no_bootstraps_res <- 
  res_list %>%
  left_join(params, by = "row") %>%
  filter(bootstrap == "none") %>%
  pull(file) %>%
  future_map_dfr(.f = readRDS) %>%
  select(-seed)
toc()

params %>%
  count(bootstrap)
nrow(no_bootstraps_res)

no_bootstraps_res %>%
  mutate(time_hrs = time / 60^2) %>%
  summarize(
    across(time_hrs, .fns = c(min = min, median = median, max = max, mean = mean, total = sum))
  ) %>%
  mutate(
    time_yrs_total = time_hrs_total / 24 / 365.25
  )

no_bootstraps_res %>% 
  select(-run_date, -time) %>%
  unnest(res) %>%
  select(
    mean_smd:omega, delta_2, bootstrap_condition = bootstrap,
    model, estimator, param, 
    K_absolute:rmse_mcse, 
    est_winsor_pct, est_winsor_pct_mcse
  )

#-------------------------------------------------------------------------------
# Check bootstrapping completeness

count_booties <- function(file) {
  x <- readRDS(file)
  x$res[[1]] %>%
    filter(param == "beta", model != "Comparison") %>%
    mutate(R = map_int(boot_CIs, \(x) if (is.null(x)) 0L else max(x$bootstraps))) %>%
    group_by(model, estimator) %>%
    count(R) %>%
    mutate(
      file = file,
      run_date = x$run_date
    )
}

plan(multisession, workers = 10L)

tic()
bootstraps_res <- 
  res_list %>%
  left_join(params, by = "row") %>%
  filter(bootstrap != "none") %>%
  pull(file) %>%
  future_map_dfr(.f = count_booties, .progress = TRUE)
toc()

plan(sequential)

bootstrap_reps <- 
  bootstraps_res %>%
  group_by(run_date, file, model) %>%
  summarize(
    complete = all(n == 25),
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
  group_by(model, complete, has_zeros) %>%
  count()

#-------------------------------------------------------------------------------
# compile results from conditions with bootstraps

bootstrap_files <- 
  res_list %>%
  left_join(params, by = "row") %>%
  select(-batch, -row, -seed) %>%
  filter(bootstrap != "none") %>%
  nest(iterations = iterations, files = file)

file_list <- 
  bootstrap_files %>%
  filter(mean_smd == 0.8, tau == 0.45, cor_mu == 0.8, m == 30, omega > 0, delta_1 == 0.2) %>%
  pull(files)

source("research/beta-function-simulations/3_performance_criteria.R")

repair_bootstraps <- function(bootstraps, B = c(49, 99, 199, 299, 399)) {
  
  if (is.null(bootstraps)) return(NULL)
  
  missing_boots <- setdiff(B, bootstraps$bootstraps)
  
  if (length(missing_boots) == 0L) return(bootstraps)
  if (length(missing_boots) == length(B)) return(NULL)
  
  bootstraps_imputed <- 
    bootstraps %>% 
    slice_max(bootstraps, n = 1L, with_ties = FALSE) %>%
    select(-bootstraps) %>%
    expand_grid(bootstraps = missing_boots)
  
  bind_rows(bootstraps, bootstraps_imputed) %>%
    filter(bootstraps %in% B)
}

summarize_bootstraps <- function(file_list) {
  dat <- map_dfr(file_list$file, .f = readRDS)
  time <- sum(dat$time)
  run_date <- min(dat$run_date)
  
  true_params <- data.frame(
    param = c("beta", "gamma", "zeta1","zeta2"),
    true_param = c(unique(dat$mean_smd), log(unique(dat$tau)^2 + unique(dat$omega)^2), log(unique(dat$delta_1)), log(unique(dat$delta_2)))
  )
  
  results <- 
    bind_rows(dat$res) %>%
    left_join(true_params, by = "param")
  
  summary_res <- 
    results %>%
    mutate(boot_CIs = map(boot_CIs, repair_bootstraps)) %>%
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

write_rds(bootstrap_res, file = "research/beta-function-simulations/sim-beta-function-bootstrap-performance-results.rds", compress = "gz", compression = 9L)
bootstrap_res <- readRDS("research/beta-function-simulations/sim-beta-function-bootstrap-performance-results.rds")

convergence_rates <- 
  bootstrap_res %>%
  select(-run_date, -time) %>%
  unnest(res) %>%
  filter(param == "beta", model == "beta") %>%
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

write_rds(timings, file = "research/beta-function-simulations/sim-beta-function-timings.rds", compress = "gz", compression = 9L)


timings %>%
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
  filter(bootstrap == "two-stage") %>%
  unnest(res) %>%
  select(
    time,  
    mean_smd:omega, delta_2, bootstrap_condition = bootstrap,
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
    mean_smd:omega, delta_2, bootstrap_condition = bootstrap,
    model, estimator, param, 
    K_absolute:rmse_mcse, 
    est_winsor_pct, est_winsor_pct_mcse
  ) %>%
  bind_rows(bootstrap_point_estimator)

params %>%
  filter(batch == 1) %>%
  count(bootstrap)
res_point_estimator %>%
  group_by(bootstrap_condition, mean_smd, tau, cor_mu, cor_sd, m, omega, delta_1) %>%
  summarize(n_res = n(), .groups = "drop") %>%
  count(bootstrap_condition)

res_point_estimator %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, m, omega, delta_1) %>%
  summarize(n_res = n(), .groups = "drop") %>%
  count(n_res)

write_rds(res_point_estimator, file = "research/beta-function-simulations/sim-beta-function-point-estimator-results.rds", compress = "gz", compression = 9L)

# variance estimation and confidence interval results

nobootstrap_conf_int <- 
  no_bootstraps_res %>%
  select(-run_date, -time) %>%
  unnest(res) %>%
  select(mean_smd:delta_2, bootstrap_condition = bootstrap, model:param, 
         K_coverage:width_mcse,
         K_relvar:rel_rmse_var_mcse, 
         var_winsor_pct, var_winsor_pct_mcse) %>%
  mutate(
    CI_type = "large-sample",
  ) 

bootstrap_large_sample <- 
  bootstrap_res %>%
  select(-run_date, -time) %>%
  unnest(res) %>%
  select(mean_smd:omega, delta_2, bootstrap_condition = bootstrap, model:param, 
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
  select(mean_smd:omega, delta_2, bootstrap_condition = bootstrap, bootstrap_type = bootstrap, model:param, 
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

params %>%
  filter(batch == 1) %>%
  count(bootstrap)

res_conf_int %>%
  group_by(bootstrap_condition, mean_smd, tau, cor_mu, cor_sd, m, omega, delta_1) %>%
  summarize(n_res = n(), .groups = "drop") %>%
  count(bootstrap_condition)

res_conf_int %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, m, omega, delta_1) %>%
  summarize(n_res = n(), .groups = "drop") %>%
  count(n_res)

write_rds(res_conf_int, file = "research/beta-function-simulations/sim-beta-function-confidence-interval-results.rds", compress = "gz", compression = 9L)
