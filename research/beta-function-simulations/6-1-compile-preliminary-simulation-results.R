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


#-------------------------------------------------------------------------------
# format results for graphing

delta_selection_levels <- c(
  "Extreme (d1=0.02, d2=0.90)" = "0.02_0.90",
  "d1=0.10, d2=0.90" = "0.10_0.90",
  "d1=0.20, d2=0.90" = "0.20_0.90",
  "d1=0.50, d2=0.90" = "0.50_0.90",
  "None (d1=1.00, d2=1.00)" = "1.00_1.00"
)

delta_selection_names <- c(
  "Extreme" = "0.02_0.90",
  "Very Strong" = "0.10_0.90",
  "Strong" = "0.20_0.90",
  "Moderate" = "0.50_0.90",
  "None" = "1.00_1.00"
)

results <- 
  no_bootstraps_res %>%
  select(-run_date, -time) %>%
  unnest(res) %>%
  mutate(
    estimator = fct(estimator, levels = c("ML","FEC","CHE","PET","PEESE","PET/PEESE")),
    estimator = fct_recode(estimator, "CHE-ISCW" = "FEC", "PML" = "ML"),
    het_ratio = omega ^ 2 / tau ^ 2,
    het_ratio = as.character(het_ratio),
    J = as.character(m),
    J = factor(J, levels = c("15", "30", "60", "90", "120")), # there is no 120
    delta_combo_key = paste0(sprintf("%.2f", delta_1), "_", sprintf("%.2f", delta_2)),
    selection_strength = factor(
      delta_combo_key,
      levels = unname(delta_selection_levels),
      labels = names(delta_selection_names)
    ),
    mu_fac = fct(as.character(mean_smd)),
    tau_fac = fct(as.character(tau), levels = c("0.05","0.15","0.3","0.45","0.6")), # there is no 0.6
    convergence = K_absolute / 2000
  ) %>%
  select(-cor_sd) %>% 
  unite("model_estimator", model:estimator, na.rm = TRUE, remove = FALSE)

mu_graph_res_main <- 
  results %>%
  filter(
    param == "beta",
    model_estimator %in% c("Comparison_CHE-ISCW","Comparison_PET/PEESE","beta_PML")
  ) %>%
  droplevels() %>%
  mutate(
    method = fct_recode(estimator, "Beta" = "PML")
  )

mu_wide_res_main <- 
  mu_graph_res_main %>%
  select(mean_smd:m, mu_fac, tau_fac, het_ratio, J, bias, var, rmse, method, selection_strength) %>%
  mutate(
    rmse_trunc = pmin(rmse, 1)
  ) %>%
  pivot_wider(
    values_from = c(bias, var, rmse, rmse_trunc), 
    names_from = method
  )


mu_graph_res_ci_main <- 
  results %>%
  filter(
    param == "beta",
    !is.na(coverage),
    model_estimator %in% c("Comparison_CHE-ISCW","Comparison_PET/PEESE","beta_PML")
  ) %>%
  droplevels() %>%
  mutate(
    method = fct_recode(estimator, "Beta" = "PML")
  )

#-------------------------------------------------------------------------------
# Bias

ggplot(mu_graph_res_main) + 
  aes(x = selection_strength, y = bias, color = method, fill = method) +
  geom_hline(yintercept = 0) +
  geom_boxplot(alpha = .5, coef = Inf) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 3))+
  facet_grid(
    tau ~ mean_smd, 
    labeller = label_bquote(
      rows = tau == .(tau),
      cols = mu == .(mean_smd)
    ),
    scales = "free_y"
  ) +
  labs(
    x = "Selection probability", 
    y = "Bias", 
    color = "",
    fill = ""
  ) + 
  theme_bw() +
  theme(legend.position = "top")

ggsave("research/beta-function-simulations/graphs/bias.png", width = 12, height = 5)

#-------------------------------------------------------------------------------
# RMSE

ggplot(mu_graph_res_main) + 
  aes(x = selection_strength, y = rmse, color = method, fill = method) +
  geom_hline(yintercept = 0) +
  geom_boxplot(alpha = .5, coef = Inf) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 3))+
  facet_grid(
    tau ~ mean_smd, 
    labeller = label_bquote(
      rows = tau == .(tau),
      cols = mu == .(mean_smd)
    ),
    scales = "free_y"
  ) +
  labs(
    x = "Selection probability", 
    y = "Root Mean-Squared Error", 
    color = "",
    fill = ""
  ) + 
  theme_bw() +
  theme(legend.position = "top")

ggsave("research/beta-function-simulations/graphs/rmse.png", width = 12, height = 5)

#-------------------------------------------------------------------------------
# Coverage

ggplot(mu_graph_res_ci_main) + 
  aes(x = J, y = coverage, color = method, fill = method) +
  geom_boxplot(alpha = .5, coef = Inf) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 1.0)) + 
  scale_y_continuous(expand = expansion(c(0,0),c(0.02,0))) + 
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(
    tau ~ mean_smd, 
    labeller = label_bquote(
      rows = tau == .(tau),
      cols = mu == .(mean_smd)
    ),
    scales = "free_y"
  ) +
  labs(
    x = "Number of studies (J)", 
    y = "Coverage rate", 
    color = "",
    fill = ""
  ) + 
  theme_bw() +
  theme(legend.position = "top")

ggsave("research/beta-function-simulations/graphs/coverage.png", width = 12, height = 5)
