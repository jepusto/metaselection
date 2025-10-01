library(tidyverse)
library(tictoc)
library(simhelpers)
library(future)
library(furrr)

rows <- 
  read_csv("research/step-function-simulations/penalized-results/batches-to-run.csv", col_names = "row") %>%
  mutate(
    i = row_number()
  )

params <- readRDS("research/step-function-simulations/penalized-results/simulation_parameters.rds")

res_list <- tibble(
  file = list.files("research/step-function-simulations/batch-results-penalized", pattern = "simulation_results_batch", full.names = TRUE)
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
# compile results from individual files

plan(multisession)

tic()
no_bootstraps_res <- 
  res_list %>%
  pull(file) %>%
  future_map_dfr(.f = readRDS) %>%
  select(-seed)
toc()

plan(sequential)

nrow(no_bootstraps_res)

#-------------------------------------------------------------------------------
# organize results for further analysis

# timing

timings <- 
  no_bootstraps_res %>%
  select(-res, -summarize_performance) %>%
  mutate(time_hrs = time / 60^2)

write_rds(timings, file = "research/step-function-simulations/penalized-results/sim-step-function-timings.rds", compress = "gz", compression = 9L)


timings %>%
  summarize(
    across(time_hrs, .fns = c(min = min, median = median, max = max, mean = mean, total = sum))
  ) %>%
  mutate(
    time_yrs_total = time_hrs_total / 24 / 365.25
  )


# point estimation results

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

write_rds(res, file = "research/step-function-simulations/penalized-results/sim-step-function-penalized-results-no-bootstraps.rds", compress = "gz", compression = 9L)

