library(tidyverse)
library(tictoc)
library(simhelpers)

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

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Compile results from conditions with no bootstraps ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


tic()
no_bootstraps_res <- 
  params %>%
  left_join(res_list, by = "row") %>%
  filter(bootstrap == "none", !is.na(file)) %>%
  select(-priors, -comparison_methods) %>%
  distinct() %>%
  pull(file) %>%
  map_dfr(.f = readRDS) %>%
  select(-seed)
toc()

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
  select(priors, iterations, model:width_mcse)

write_rds(res, file = "research/step-function-simulations/sim-step-function-results-no-bootstraps.rds", compress = "gz", compression = 9L)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Set up summary calculations for conditions with bootstraps ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


bootstrap_files <-
  params %>%
  filter(bootstrap != "none") %>%
  select(-seed) %>%
  mutate(file = paste0("simulation_results_batch", row, ".rds")) %>%
  nest(batches = batch, iterations = iterations, rows = row, files = file) %>%
  mutate(
    R_max = map_dbl(R, max),
    nbatches = map_dbl(batches, nrow),
    row = row_number(),
    files = map_chr(files, \(x) paste(x$file, collapse = ", "))
  )

bootstrap_files %>%
  filter(bootstrap == "exponential") %>%
  select(row, files) %>%
  write_tsv("research/step-function-simulations/bootstrap-batches-to-run.tsv", col_names = FALSE)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# compile results from conditions with bootstraps ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

bootstrap_res_list <- tibble(
  file = list.files("research/step-function-simulations/batch-results", pattern = "simulation_results_bootstrap_batch", full.names = TRUE)
) %>%
  mutate(
    row = str_extract(file, "batch[0-9]+.rds") |> str_sub(6,-5) |> as.integer()
  )

nrow(bootstrap_res_list)

tic()
bootstrap_res <- 
  bootstrap_files %>%
  mutate(row = map_int(rows, \(x) x$row[1])) %>%
  select(-batches, -rows, -files, -nbatches) %>%
  inner_join(bootstrap_res_list, by = "row") %>%
  mutate(
    iterations = map_int(iterations, \(x) sum(x$iterations)),
    res =  map(file, .f = read_rds, .progress = TRUE)
  ) %>%
  select(-file, -row) %>%
  unnest(res)
toc()

bootstrap_res %>% count(iterations)

bootstrap_res %>%
  select(-run_date, -time) %>%
  unnest(res) %>%
  filter(estimator != "CML") %>%
  select(
    mean_smd:psi, bootstrap, omega, steps, bootstrap_type = bootstrap, model:param,
    bootstraps, extrapolated, boot_coverage, boot_coverage_mcse, boot_width, boot_width_mcse
  ) %>%
  unnest(
    c(bootstraps, extrapolated, boot_coverage, boot_coverage_mcse, boot_width, boot_width_mcse),
    names_sep = "-"
  ) %>%
  pivot_longer(
    starts_with("boot_"),
    names_to = c(".value", "CI_type"),
    names_pattern = "(.+)-(.+)"
  ) %>%
  rename_with(~ str_remove(.x, "^boot_")) %>%
  summary()

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
