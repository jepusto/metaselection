library(tidyverse)


selection_levels <- c(
  "0.02 (Strong)" = 0.02,
  "0.05" = 0.05,
  "0.10" = 0.10,
  "0.20" = 0.20,
  "0.50" = 0.50,
  "1.00 (None)" = 1.00
)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Compile convergence results including the non-contingent CML estimator ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

convergence_results <- 
  readRDS(here::here("research","step-function-simulations","sim-step-function-results-no-bootstraps.rds")) %>%
  filter(
    estimator %in% c("ARGL","CML","CHE-ISCW","PET/PEESE"), 
    priors == "Weak",
    param == "beta"
  ) %>% 
  select(-param, -priors) %>%
  mutate(
    estimator = fct(estimator, levels = c("CML","ARGL","CHE","CHE-ISCW","PET","PEESE","PET/PEESE")),
    estimator = fct_recode(estimator, "PML" = "CML"),
    N_factor = fct(if_else(n_multiplier < 1, "Small", "Typical")),
    psi_fac = as.character(psi) |> fct() |> fct_recode("Independent" = "0", "Dependent" = "1"),
    weight_num = weight,
    weight = as.character(weight),
    het_ratio = omega ^ 2 / tau ^ 2,
    het_ratio = as.character(het_ratio),
    scrmse = sqrt(m) * rmse, 
    J = as.character(m),
    J = factor(J, levels = c("15", "30", "60", "90", "120")),
    weight = factor(
      weight, 
      levels = selection_levels,
      labels = names(selection_levels)
    ),
    mu_fac = fct(as.character(mean_smd)),
    tau_fac = fct(as.character(tau), levels = c("0.05","0.15","0.3","0.45","0.6")),
    convergence = 100 * K_absolute / iterations
  )

convergence_results %>%
  filter(model == "3PSM") %>%
  group_by(estimator) %>%
  summarize(
    min_convergence = min(convergence)
  )

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Compile point estimator performance results ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

results <- 
  readRDS(here::here("research", "step-function-simulations","sim-step-function-results-no-bootstraps.rds")) %>%
  filter(estimator != "CML", priors == "Weak") %>%
  mutate(
    estimator = case_match(estimator, 'CML-fallback' ~ "CML", .default = estimator),
    estimator = fct(estimator, levels = c("CML","CML-model","ARGL","CHE","CHE-ISCW","PET","PEESE","PET/PEESE")),
    estimator = fct_recode(estimator, "PML" = "CML", "PML-model" = "CML-model"),
    N_factor = fct(if_else(n_multiplier < 1, "Small", "Typical")),
    psi_fac = as.character(psi) |> fct() |> fct_recode("Independent" = "0", "Dependent" = "1"),
    weight_num = weight,
    weight = as.character(weight),
    het_ratio = omega ^ 2 / tau ^ 2,
    het_ratio = as.character(het_ratio),
    scrmse = sqrt(m) * rmse, 
    J = as.character(m),
    J = factor(J, levels = c("15", "30", "60", "90", "120")),
    weight = factor(
      weight, 
      levels = selection_levels,
      labels = names(selection_levels)
    ),
    mu_fac = fct(as.character(mean_smd)),
    tau_fac = fct(as.character(tau), levels = c("0.05","0.15","0.3","0.45","0.6")),
    convergence = K_absolute / iterations
  ) %>%
  select(-cor_sd, -n_multiplier)

mu_graph_res <- 
  results %>%
  filter(
    param == "beta",
    estimator %in% c("CHE-ISCW","PET/PEESE","PML","ARGL")
  )

mu_wide_res <- 
  mu_graph_res %>%
  select(mean_smd:psi, psi_fac, mu_fac, tau_fac, N_factor, het_ratio, J, estimator, bias, var, rmse, scrmse) %>%
  mutate(
    rmse_trunc = pmin(rmse, 1)
  ) %>%
  pivot_wider(
    values_from = c(bias, var, rmse, rmse_trunc, scrmse), 
    names_from = estimator
  )

tau2_graph_res <- 
  results %>%
  filter(
    param == "tau2"
  ) %>%
  mutate(
    total_var = tau^2 + omega^2,
    rbias = bias / total_var,
    rbias_mcse = bias_mcse / total_var,
    rvar = var / total_var^2,
    rvar_mcse = var_mcse / total_var^2,
    rrmse = rmse / total_var,
    rrmse_mcse = rmse_mcse / total_var
  )

tau2_wide_res <- 
  tau2_graph_res %>%
  filter(estimator %in% c("ARGL","PML","CHE")) %>%
  select(mean_smd:psi, psi_fac, mu_fac, tau_fac, N_factor, het_ratio, J, estimator, rbias, rvar, rrmse) %>%
  pivot_wider(
    values_from = c(rbias, rvar, rrmse), 
    names_from = estimator
  )

zeta_graph_res <- 
  results %>%
  filter(
    param == "zeta1"
  )

zeta_wide_res <- 
  zeta_graph_res %>%
  select(mean_smd:psi, psi_fac, mu_fac, tau_fac, N_factor, het_ratio, J, estimator, bias, var, rmse) %>%
  pivot_wider(
    values_from = c(bias, var, rmse), 
    names_from = estimator
  )

lambda_graph_res <- 
  results %>%
  filter(
    param == "lambda1"
  ) %>%
  mutate(
    rbias = bias / weight_num,
    rbias_mcse = bias_mcse / weight_num,
    rvar = var / weight_num^2,
    rvar_mcse = var_mcse / weight_num^2,
    rrmse = rmse / weight_num,
    rrmse_mcse = rmse_mcse / weight_num
  )

lambda_wide_res <- 
  lambda_graph_res %>%
  select(mean_smd:psi, psi_fac, mu_fac, tau_fac, N_factor, het_ratio, J, estimator, rbias, rvar, rrmse) %>%
  pivot_wider(
    values_from = c(rbias, rvar, rrmse), 
    names_from = estimator
  )

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Compile confidence interval performance results ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

nobootstrap_conf_int <- 
  readRDS(here::here("research", "step-function-simulations","sim-step-function-results-no-bootstraps.rds")) %>%
  filter(estimator != "CML") %>%
  select(mean_smd:psi, bootstrap, omega, steps, model:param, K_coverage:width_mcse) %>%
  mutate(
    CI_type = "large-sample"
  ) 

bootstrap_conf_int <- 
  readRDS(here::here("research", "step-function-simulations","sim-step-function-bootstrap-performance-results.rds")) %>%
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
  rename_with(~ str_remove(.x, "^boot_"))


results_ci <- 
  bind_rows(
    nobootstrap_conf_int,
    bootstrap_conf_int
  ) %>%
  group_by(mean_smd, tau, omega, cor_mu, cor_sd, weight, m, n_multiplier, psi, steps) %>%
  mutate(
    bootstrap_condition = any(!is.na(bootstrap_type))
  ) %>%
  ungroup() %>%
  mutate(
    estimator = case_match(estimator, 'CML-fallback' ~ "CML", .default = estimator),
    estimator = fct(estimator, levels = c("CML","CML-model","ARGL","CHE","CHE-ISCW","PET","PEESE","PET/PEESE")),
    estimator = fct_recode(estimator, "PML" = "CML", "PML-model" = "CML-model"),
    N_factor = fct(if_else(n_multiplier < 1, "Small", "Typical")),
    psi_fac = as.character(psi) |> fct() |> fct_recode("Independent" = "0", "Dependent" = "1"),
    weight_num = weight,
    weight = as.character(weight),
    het_ratio = omega ^ 2 / tau ^ 2,
    het_ratio = as.character(het_ratio),
    J = as.character(m),
    J = factor(J, levels = c("15", "30", "60", "90", "120")),  # 120 and 200 not there 
    weight = factor(
      weight, 
      levels = selection_levels,
      labels = names(selection_levels)
    ),
    mu_fac = fct(as.character(mean_smd)),
    tau_fac = fct(as.character(tau), levels = c("0.05","0.15","0.3","0.45","0.6")),
    bootstrap_type = recode(bootstrap_type, .missing = "none"),
    CI_boot_method = case_when(
      estimator == "PML-model" ~ "model-based",
      CI_type == "large-sample" ~ "cluster-robust", 
      TRUE ~ paste(CI_type, " (", bootstrap_type, ")", sep = "")
    ),
    CI_boot_method = fct(
      CI_boot_method,
      levels = c("cluster-robust","model-based",
                 "percentile (two-stage)","percentile (multinomial)","percentile (exponential)",
                 "biascorrected (two-stage)","biascorrected (multinomial)","biascorrected (exponential)",
                 "BCa (two-stage)","BCa (multinomial)","BCa (exponential)",
                 "basic (two-stage)","basic (multinomial)","basic (exponential)",
                 "student (two-stage)","student (multinomial)","student (exponential)")
    )
  ) %>%
  select(-cor_sd, -n_multiplier)


mu_graph_res_ci <- 
  results_ci %>%
  filter(
    param == "beta",
    !is.na(coverage),
    is.na(bootstraps) | (bootstraps == 1999 & extrapolated),
    estimator %in% c("CHE-ISCW","PET/PEESE","PML","ARGL")
  )


gamma_graph_res_ci <- 
  results_ci %>%
  filter(
    param == "gamma",
    !is.na(coverage),
    is.na(bootstraps) | (bootstraps == 1999 & extrapolated),
    estimator %in% c("CHE-ISCW","PML","ARGL")
  )

zeta_graph_res_ci <- 
  results_ci %>%
  filter(
    param == "zeta1",
    !is.na(coverage),
    is.na(bootstraps) | (bootstraps == 1999 & extrapolated),
    estimator %in% c("PML","ARGL")
  )


performance_plot <- function(data, measure, label) {
  ggplot(data) + 
    aes(x = weight, y = {{measure}}, color = estimator, fill = estimator) +
    geom_hline(yintercept = 0) +
    geom_boxplot(alpha = .5, coef = Inf) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 3))+
    facet_grid(
      tau ~ mean_smd, 
      labeller = label_bquote(
        rows = tau[B] == .(tau),
        cols = mu == .(mean_smd)
      ),
      scales = "free_y"
    ) +
    labs(
      x = "Selection probability", 
      y = label, 
      color = "Estimator",
      fill = "Estimator"
    ) + 
    theme_bw() +
    theme(legend.position = "top")
}

  
RMSE_comparison_plot <- function(
    data, 
    x_method, y_method, 
    measure = "rmse", 
    col_factor = psi_fac, 
    col_lab = "Selection process", 
    legend_rows = 1L
) {
  
  y_lab <- paste0("RMSE ratio (",y_method, " / ", x_method, ")")
  x_var <- sym(paste(measure, x_method, sep = "_"))
  y_var <- sym(paste(measure, y_method, sep = "_"))
  
  ggplot(data) + 
    aes(x = weight, y = {{y_var}} / {{x_var}}, shape = {{col_factor}}, color = {{col_factor}}, fill = {{col_factor}}) +
    geom_hline(yintercept = 1) + 
    geom_point(alpha = .5, position = position_jitterdodge()) +
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
      color = col_lab,
      shape = col_lab,
      fill = col_lab,
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
      weight_num ~ mean_smd, 
      labeller = label_bquote(
        rows = lambda[1] == .(weight_num),
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


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Comparison of extrapolation versus bootstraps with large B ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

big_B_boot_CIs <- 
  results_ci %>%
  filter(
    param == "beta",
    estimator %in% c("PML","ARGL"),
    bootstrap_type == "two-stage",
    CI_type %in% c("percentile","basic","BCa")
  ) %>%
  group_by(mean_smd, tau, cor_mu, psi, weight, m, N_factor, het_ratio) %>%
  filter(any(bootstraps == 1999 & !extrapolated)) %>%
  mutate(
    cover_lo = coverage - qnorm(0.975) * coverage_mcse,
    cover_hi = coverage + qnorm(0.975) * coverage_mcse,
    tau_lab = paste("tau ==", tau),
    lambda_lab = paste("lambda ==", weight),
    CI_lab = paste("CI type:", CI_type)
  )


bootstrap_Bs <- unique(results_ci$bootstraps)[-1]
