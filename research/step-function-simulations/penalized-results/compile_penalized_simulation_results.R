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

#-------------------------------------------------------------------------------
# compile results from individual files

tic()
no_bootstraps_res <- 
  res_list %>%
  pull(file) %>%
  future_map_dfr(.f = readRDS) %>%
  select(-seed)
toc()

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
    mean_smd:omega, steps, 
    model, estimator, param, 
    K_absolute:rmse_mcse, 
    K_coverage:width_mcse
  )

res %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, weights, m, n_multiplier, omega, steps) %>%
  summarize(n_res = n(), .groups = "drop") %>%
  count(n_res)

write_rds(res, file = "research/step-function-simulations/penalized-results/sim-step-function-penalized-results.rds", compress = "gz", compression = 9L)


res %>%
  filter(param == "beta") %>%
  arrange(mean_smd, tau, omega, cor_mu, weights, m, n_multiplier, estimator, priors) %>%


res_beta <- 
  res %>%
  filter(param == "beta") 

ggplot(res_beta) + 
  aes(factor(weights), bias, color = interaction(priors, estimator), fill = interaction(priors, estimator)) + 
  geom_boxplot() + 
  facet_grid(mean_smd ~ tau) + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "bias",
    color = "Estimator",
    fill = "Estimator"
  )

ggplot(res_beta) + 
  aes(factor(weights), rmse, color = interaction(priors, estimator), fill = interaction(priors, estimator)) + 
  geom_boxplot() + 
  facet_grid(mean_smd ~ tau) + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "RMSE",
    color = "Estimator",
    fill = "Estimator"
  )

res_beta_wide <- 
  res_beta %>%
  select(mean_smd, tau, omega, cor_mu, weights, m, n_multiplier, estimator, priors, bias, rmse) %>%
  unite("estimator_prior", estimator, priors) %>%
  pivot_wider(names_from = estimator_prior, values_from = c(bias, rmse))

res_beta_wide %>%
  ggplot() + 
  aes(bias_CML_None, bias_CML_Default, color = factor(m)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point() + 
  facet_grid(mean_smd ~ tau, scales = "free") + 
  theme_light() + 
  labs(
    x = "Flat", 
    y = "Default",
    color = "Sample size"
  )

res_beta_wide %>%
  ggplot() + 
  aes(rmse_CML_None, rmse_CML_Default, color = factor(m)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point() + 
  facet_grid(mean_smd ~ tau, scales = "free") + 
  theme_light() + 
  labs(
    x = "Flat", 
    y = "Default",
    color = "Sample size"
  )


res_beta_wide %>%
  ggplot() + 
  aes(bias_ARGL_None, bias_ARGL_Default, color = factor(m)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point() + 
  facet_grid(mean_smd ~ tau, scales = "free") + 
  theme_light() + 
  labs(
    x = "Flat", 
    y = "Default",
    color = "Sample size"
  )


res_beta_wide %>%
  ggplot() + 
  aes(rmse_ARGL_None, rmse_ARGL_Default, color = factor(m)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point() + 
  facet_grid(mean_smd ~ tau, scales = "free") + 
  theme_light() + 
  labs(
    x = "Flat", 
    y = "Default",
    color = "Sample size"
  )

res_beta_wide %>%
  ggplot() + 
  aes(rmse_CML_Default, rmse_ARGL_Default, color = factor(m)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point() + 
  facet_grid(mean_smd ~ tau, scales = "free") + 
  theme_light() + 
  labs(
    x = "CML", 
    y = "ARGL",
    color = "Sample size"
  )
