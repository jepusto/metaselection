library(tidyverse)

delta_selection_levels <- c(
  "Strong (d1=0.01, d2=0.90)" = "0.01_0.90",
  "d1=0.20, d2=0.90" = "0.20_0.90",
  "d1=0.50, d2=0.90" = "0.50_0.90",
  "d1=0.80, d2=0.90" = "0.80_0.90",
  "None (d1=1.00, d2=1.00)" = "1.00_1.00"
)

delta_selection_names <- c(
  "Very Strong" = "0.01_0.90",
  "Strong" = "0.20_0.90",
  "Mod" = "0.50_0.90",
  "Weak" = "0.80_0.90",
  "None" = "1.00_1.00"
)

results <- 
  #readRDS("../beta-function-simulations/sim-beta-function-point-estimator-results.rds") %>%
  readRDS("C:/GitHub/metaselection/research/beta-function-simulations/sim-beta-function-point-estimator-results.rds") %>%
  mutate(
    estimator = fct(estimator, levels = c("ML","FEC","CHE","PET","PEESE","PET/PEESE")),
    estimator = fct_recode(estimator, "CHE-ISCW" = "FEC", "CML" = "ML"),
    het_ratio = omega ^ 2 / tau ^ 2,
    het_ratio = as.character(het_ratio),
    scrmse = sqrt(m) * rmse, 
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
    convergence = K_absolute / 2000,
    winz_convergence = (1 - winsor_pct) * K_absolute / 2000
  ) %>%
  select(-cor_sd) %>% 
  unite("model_estimator", model:estimator, na.rm = TRUE, remove = FALSE)

mu_graph_res_main <- 
  results %>%
  filter(
    param == "beta",
    model_estimator %in% c("Comparison_CHE-ISCW","Comparison_PET/PEESE","beta_CML")
  ) %>%
  droplevels() %>%
  mutate(
    method = fct_recode(estimator, "Beta" = "CML")
  )

mu_wide_res_main <- 
  mu_graph_res_main %>%
  select(mean_smd:m, mu_fac, tau_fac, het_ratio, J, bias, var, rmse, scrmse, method, selection_strength) %>%
  mutate(
    rmse_trunc = pmin(rmse, 1)
  ) %>%
  pivot_wider(
    values_from = c(bias, var, rmse, rmse_trunc, scrmse), 
    names_from = method
  )

mu_graph_res_miss <- 
  results %>%
  filter(
    param == "beta",
    estimator %in% c("CML")
  ) %>%
  droplevels() %>%
  rename(method = model) %>%
  mutate(
    method = fct(method, levels = c("beta","3PSM","4PSM")) |> fct_recode("one-step" = "3PSM", "two-step" = "4PSM")
  )

mu_wide_res_miss <- 
  mu_graph_res_miss %>%
  select(mean_smd:m, mu_fac, tau_fac, het_ratio, J, bias, var, rmse, scrmse, method, selection_strength) %>%
  mutate(
    rmse_trunc = pmin(rmse, 1)
  ) %>%
  pivot_wider(
    values_from = c(bias, var, rmse, rmse_trunc, scrmse), 
    names_from = method
  )

gamma_graph_res_miss <- # gamma only estimated for the 3 CML methods
  results %>%
  filter(
    param == "gamma"
  ) %>%
  mutate(
    scrmse_trunc = pmin(scrmse, 3 / tau + rnorm(n(), sd = 0.01))
  ) %>%
  droplevels()%>%
  rename(method = model) %>%
  mutate(
    method = fct(method, levels = c("beta","3PSM","4PSM")) |> fct_recode("one-step" = "3PSM", "two-step" = "4PSM")
  )

gamma_wide_res_miss <- 
  gamma_graph_res_miss %>%
  select(mean_smd:m, mu_fac, tau_fac, het_ratio, J, bias, var, rmse, scrmse, scrmse_trunc, method) %>%
  pivot_wider(
    values_from = c(bias, var, rmse, scrmse, scrmse_trunc), 
    names_from = method
  )


results_ci <- 
  readRDS("../beta-function-simulations/sim-beta-function-confidence-interval-results.rds") %>%
  mutate(
    estimator = fct(estimator, levels = c("ML","FEC","CHE","PET","PEESE","PET/PEESE")),
    estimator = fct_recode(estimator, "CHE-ISCW" = "FEC", "CML" = "ML"),
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
    bootstrap_type = recode(bootstrap_type, .missing = "none"),
    CI_boot_method = if_else(
      CI_type == "large-sample",
      "cluster-robust", 
      paste(CI_type, " (", bootstrap_type, ")", sep = "")
    ),
    CI_boot_method = fct(
      CI_boot_method,
      levels = c("cluster-robust",
                 "percentile (two-stage)",
                 "biascorrected (two-stage)",
                 "basic (two-stage)",
                 "student (two-stage)")
    )
  ) %>%
  select(-cor_sd) %>% 
  unite("model_estimator", model:estimator, na.rm = TRUE, remove = FALSE)

mu_graph_res_ci_main <- 
  results_ci %>%
  filter(
    param == "beta",
    !is.na(coverage),
    is.na(bootstraps) | bootstraps == 1999,
    model_estimator %in% c("Comparison_CHE-ISCW","Comparison_PET/PEESE","beta_CML")
  ) %>%
  droplevels() %>%
  mutate(
    method = fct_recode(estimator, "Beta" = "CML")
  )

mu_graph_res_ci_miss <- 
  results_ci %>%
  filter(
    param == "beta",
    !is.na(coverage),
    is.na(bootstraps) | bootstraps == 1999,
    estimator %in% c("CML")
  ) %>%
  droplevels() %>%
  rename(method = model) %>%
  mutate(
    method = fct(method, levels = c("beta","3PSM","4PSM")) |> fct_recode("one-step" = "3PSM", "two-step" = "4PSM")
  )

gamma_graph_res_ci_main <- 
  results_ci %>%
  filter(
    param == "gamma",
    !is.na(coverage),
    is.na(bootstraps) | bootstraps == 1999,
    model_estimator %in% c("Comparison_CHE-ISCW","Comparison_PET/PEESE","beta_CML")
  ) %>%
  droplevels()

gamma_graph_res_ci_miss <- 
  results_ci %>%
  filter(
    param == "gamma",
    !is.na(coverage),
    is.na(bootstraps) | bootstraps == 1999,
    estimator %in% c("CML")
  ) %>%
  droplevels() %>%
  rename(method = model) %>%
  mutate(
    method = fct(method, levels = c("beta","3PSM","4PSM")) |> fct_recode("one-step" = "3PSM", "two-step" = "4PSM")
  )


RMSE_comparison_plot <- function(data, x_method, y_method, col_factor = J, col_lab = "Number of studies (J)", legend_rows = 1L) {
  
  y_lab <- paste0("RMSE ratio (",y_method, " / ", x_method, ")")
  x_var <- sym(paste("scrmse", x_method, sep = "_"))
  y_var <- sym(paste("scrmse", y_method, sep = "_"))
  
  ggplot(data) + 
    aes(x = selection_strength, y = {{y_var}} / {{x_var}}, shape = {{col_factor}}, color = {{col_factor}}) +
    geom_hline(yintercept = 1) + 
    geom_point(alpha = .5, position = position_jitter(width = 0.2)) +
    expand_limits(y = 0.5) + 
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 3))+
    scale_y_continuous(transform = "log2") + 
    scale_color_brewer(palette = "Dark2", guide = guide_legend(nrow=legend_rows)) +
    facet_grid(tau ~ mean_smd, 
               labeller = label_bquote(rows = tau == .(tau),
                                       cols = mu == .(mean_smd)),
               scales = "free_y"
    ) +
    labs(
      x = "Selection probability",
      y = y_lab,
      shape = col_lab, color = col_lab
    ) + 
    theme_bw() +
    theme(legend.position = "top")
}

coverage_plot <- function(data) {
  
  ggplot(data, aes(x = J, y = coverage, color = CI_type, fill = CI_type)) +
    geom_boxplot(alpha = .5, coef = Inf) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    scale_y_continuous(limits = c(NA, 1), expand = expansion(0,0)) + 
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
      color = "Method",
      fill = "Method"
    ) + 
    theme_bw() +
    theme(legend.position = "top")
}
