library(tidyverse)

selection_levels <- c(
  "0.02 (Strong)" = 0.02,
  "0.05" = 0.05,
  "0.10" = 0.10,
  "0.20" = 0.20,
  "0.50" = 0.50,
  "1.00 (None)" = 1.00
)


results <- 
  readRDS("../step-function-simulations/sim-step-function-point-estimator-results.rds") %>%
  mutate(
    estimator = fct(estimator, levels = c("CML","ARGL","CHE","CHE-ISCW","PET","PEESE","PET/PEESE")),
    N_factor = fct(if_else(n_multiplier < 1, "Small", "Typical")),
    weights = as.character(weights),
    het_ratio = omega ^ 2 / tau ^ 2,
    het_ratio = as.character(het_ratio),
    scrmse = sqrt(m) * rmse, 
    J = as.character(m),
    J = factor(J, levels = c("15", "30", "60", "90", "120")),
    weights = factor(
      weights, 
      levels = selection_levels,
      labels = names(selection_levels)
    ),
    mu_fac = fct(as.character(mean_smd)),
    tau_fac = fct(as.character(tau), levels = c("0.05","0.15","0.3","0.45","0.6")),
    convergence = K_absolute / 2000,
    winz_convergence = (1 - winsor_pct) * K_absolute / 2000
  ) %>%
  select(-cor_sd, -n_multiplier)

mu_graph_res <- 
  results %>%
  filter(
    param == "beta",
    estimator %in% c("CHE-ISCW","PET/PEESE","CML","ARGL")
  )

mu_wide_res <- 
  mu_graph_res %>%
  select(mean_smd:m, mu_fac, tau_fac, N_factor, het_ratio, J, estimator, bias, var, rmse, scrmse) %>%
  mutate(
    rmse_trunc = pmin(rmse, 1)
  ) %>%
  pivot_wider(
    values_from = c(bias, var, rmse, rmse_trunc, scrmse), 
    names_from = estimator
  )

gamma_graph_res <- 
  results %>%
  filter(
    param == "gamma"
  ) %>%
  mutate(
    scrmse_trunc = pmin(scrmse, 3 / tau + rnorm(n(), sd = 0.01))
  )

gamma_wide_res <- 
  gamma_graph_res %>%
  select(mean_smd:m, mu_fac, tau_fac, N_factor, het_ratio, J, estimator, bias, var, rmse, scrmse,scrmse_trunc) %>%
  pivot_wider(
    values_from = c(bias, var, rmse, scrmse, scrmse_trunc), 
    names_from = estimator
  )

zeta_graph_res <- 
  results %>%
  filter(
    param == "zeta1"
  )

zeta_wide_res <- 
  zeta_graph_res %>%
  select(mean_smd:m, mu_fac, tau_fac, N_factor, het_ratio, J, estimator, bias, var, rmse, scrmse) %>%
  pivot_wider(
    values_from = c(bias, var, rmse, scrmse), 
    names_from = estimator
  )


results_ci <- 
  readRDS("../step-function-simulations/sim-step-function-confidence-interval-results.rds") %>%
  mutate(
    estimator = fct(estimator, levels = c("CML","ARGL","CHE","CHE-ISCW","PET","PEESE","PET/PEESE")),
    N_factor = fct(if_else(n_multiplier < 1, "Small", "Typical")),
    weights = as.character(weights),
    het_ratio = omega ^ 2 / tau ^ 2,
    het_ratio = as.character(het_ratio),
    J = as.character(m),
    J = factor(J, levels = c("15", "30", "60", "90", "120")),  # 120 and 200 not there 
    weights = factor(
      weights, 
      levels = selection_levels,
      labels = names(selection_levels)
    ),
    mu_fac = fct(as.character(mean_smd)),
    tau_fac = fct(as.character(tau), levels = c("0.05","0.15","0.3","0.45","0.6")),
    bootstrap_type = recode(bootstrap_type, .missing = "none"),
    CI_boot_method = if_else(
      CI_type == "large-sample",
      "cluster-robust", 
      paste(CI_type, " (", bootstrap_type, ")", sep = "")
    ),
    CI_boot_method = fct(
      CI_boot_method,
      levels = c("cluster-robust",
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
    is.na(bootstraps) | bootstraps == 1999,
    estimator %in% c("CHE-ISCW","PET/PEESE","CML","ARGL")
  )

gamma_graph_res_ci <- 
  results_ci %>%
  filter(
    param == "gamma",
    !is.na(coverage),
    is.na(bootstraps) | bootstraps == 1999,
    estimator %in% c("CHE-ISCW","CML","ARGL")
  )

zeta_graph_res_ci <- 
  results_ci %>%
  filter(
    param == "zeta1",
    !is.na(coverage),
    is.na(bootstraps) | bootstraps == 1999,
    estimator %in% c("CML","ARGL")
  )


RMSE_comparison_plot <- function(data, x_method, y_method, col_factor = J, col_lab = "Number of studies (J)", legend_rows = 1L) {
  
  y_lab <- paste0("RMSE ratio (",y_method, " / ", x_method, ")")
  x_var <- sym(paste("scrmse", x_method, sep = "_"))
  y_var <- sym(paste("scrmse", y_method, sep = "_"))
  
  ggplot(data) + 
    aes(x = weights, y = {{y_var}} / {{x_var}}, shape = {{col_factor}}, color = {{col_factor}}) +
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


#-------------------------------------------------------------------------------
# Comparison of extrapolation versus bootstraps with large B

boot_real <- 
  results_ci %>%
  filter(
    mean_smd %in% c(0.2), 
    m == 15, 
    bootstrap_condition == "bootstrap", 
    param == "beta",
    estimator %in% c("CML","ARGL"),
    bootstrap_type == "two-stage",
    weights %in% c("0.05","0.20"),
    N_factor == "Typical",
    het_ratio == 0.5,
  ) %>%
  mutate(
    cover_lo = coverage - qnorm(0.975) * coverage_mcse,
    cover_hi = coverage + qnorm(0.975) * coverage_mcse,
    tau_lab = paste("tau ==", tau),
    lambda_lab = paste("lambda ==", weights),
    CI_lab = paste("CI type:", CI_type)
  )

big_B_bootstraps <- 
  read_rds(file = "../step-function-simulations/sim-step-function-big-B-bootstrap-performance-results.rds") %>%
  unnest(res) %>%
  filter(
    param == "beta",
  ) %>%
  mutate(
    bootstraps = 1999L,
    estimator = fct(estimator, levels = c("CML","ARGL","CHE","CHE-ISCW","PET","PEESE","PET/PEESE")),
    cover_lo = coverage - qnorm(0.975) * coverage_mcse,
    cover_hi = coverage + qnorm(0.975) * coverage_mcse,
    weights = as.character(weights),
    weights = factor(
      weights, 
      levels = selection_levels,
      labels = names(selection_levels)
    ),
    tau_lab = paste("tau ==", tau),
    lambda_lab = paste("lambda ==", weights),
    CI_lab = paste("CI type:", CI_type)
  )


boot_compare <- 
  bind_rows(
    extrapolated = boot_real,
    direct = big_B_bootstraps,
    .id = "coverage_estimator"
  ) %>%
  filter(
    bootstraps == 1999L
  ) %>%
  select(tau, weights, estimator, CI_type, coverage_estimator, coverage, coverage_mcse) %>%
  pivot_wider(
    names_from = coverage_estimator,
    values_from = c(coverage, coverage_mcse)
  ) %>%
  mutate(
    diff = 100 * (coverage_extrapolated - coverage_direct),
    diff_mcse = 100 * sqrt(coverage_mcse_extrapolated^2 + coverage_mcse_direct^2)
  ) %>%
  summarize(
    mean = mean(diff), 
    min = min(diff),
    max = max(diff),
    rmse = sqrt(mean(diff^2))
  )
