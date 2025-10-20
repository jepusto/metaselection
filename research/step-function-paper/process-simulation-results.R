library(tidyverse)


selection_levels <- c(
  "0.02 (Strong)" = 0.02,
  "0.05" = 0.05,
  "0.10" = 0.10,
  "0.20" = 0.20,
  "0.50" = 0.50,
  "1.00 (None)" = 1.00
)

#-------------------------------------------------------------------------------
# Compile convergence results including the non-contingent CML estimator

convergence_results <- 
  readRDS("../step-function-simulations/sim-step-function-results-no-bootstraps.rds") %>%
  filter(
    estimator %in% c("ARGL","CML","CHE-ISCW","PET/PEESE"), 
    priors == "Weak",
    param == "beta"
  ) %>% 
  select(-param, -priors) %>%
  mutate(
    estimator = fct(estimator, levels = c("CML","ARGL","CHE","CHE-ISCW","PET","PEESE","PET/PEESE")),
    N_factor = fct(if_else(n_multiplier < 1, "Small", "Typical")),
    weights_num = weights,
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
    convergence = 100 * K_absolute / iterations
  )

convergence_results %>%
  filter(model == "3PSM") %>%
  group_by(estimator) %>%
  summarize(
    min_convergence = min(convergence)
  )

#-------------------------------------------------------------------------------
# Compile point estimator performance results 
# Using CML-fallback estimator and ARGL estimator, both with weak priors

results <- 
  readRDS("../step-function-simulations/sim-step-function-results-no-bootstraps.rds") %>%
  filter(estimator != "CML", priors == "Weak") %>%
  mutate(
    estimator = case_match(estimator, 'CML-fallback' ~ "CML", .default = estimator),
    estimator = fct(estimator, levels = c("CML","ARGL","CHE","CHE-ISCW","PET","PEESE","PET/PEESE")),
    N_factor = fct(if_else(n_multiplier < 1, "Small", "Typical")),
    weights_num = weights,
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
    convergence = K_absolute / iterations
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

tau2_graph_res <- 
  results %>%
  filter(
    param == "tau2"
  ) %>%
  mutate(
    rbias = bias / tau^2,
    rbias_mcse = bias_mcse / tau^2,
    rvar = var / tau^4,
    rvar_mcse = var_mcse / tau^4,
    rrmse = rmse / tau^2,
    rrmse_mcse = rmse_mcse / tau^2,
    scrrmse = scrmse / tau^2
  )

tau2_wide_res <- 
  tau2_graph_res %>%
  filter(estimator %in% c("ARGL","CML","CHE")) %>%
  select(mean_smd:m, mu_fac, tau_fac, N_factor, het_ratio, J, estimator, rbias, rvar, rrmse, scrmse) %>%
  pivot_wider(
    values_from = c(rbias, rvar, rrmse, scrmse), 
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

lambda_graph_res <- 
  results %>%
  filter(
    param == "lambda1"
  ) %>%
  mutate(
    rbias = bias / weights_num,
    rbias_mcse = bias_mcse / weights_num,
    rvar = var / weights_num^2,
    rvar_mcse = var_mcse / weights_num^2,
    rrmse = rmse / weights_num,
    rrmse_mcse = rmse_mcse / weights_num,
    scrrmse = scrmse / weights_num
  )

lambda_wide_res <- 
  lambda_graph_res %>%
  select(mean_smd:m, mu_fac, tau_fac, N_factor, het_ratio, J, estimator, rbias, rvar, rrmse, scrrmse) %>%
  pivot_wider(
    values_from = c(rbias, rvar, rrmse, scrrmse), 
    names_from = estimator
  )

#-------------------------------------------------------------------------------
# Compile confidence interval performance results 
# Using CML-fallback estimator and ARGL estimator, both with weak priors

nobootstrap_conf_int <- 
  readRDS("../step-function-simulations/sim-step-function-results-no-bootstraps.rds") %>%
  filter(estimator != "CML", priors == "Weak") %>%
  select(mean_smd:steps, model:param, 
         K_coverage:width_mcse) %>%
  mutate(
    CI_type = "large-sample"
  ) 

bootstrap_conf_int <- 
  readRDS("../step-function-simulations/sim-step-function-bootstrap-performance-results.rds") %>%
  select(-run_date, -time) %>%
  unnest(res) %>%
  filter(estimator != "CML", priors == "Weak") %>%
  select(mean_smd:steps, bootstrap_type = bootstrap, model:param, 
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
  rename_with(~ str_remove(.x, "^boot_"))


results_ci <- 
  bind_rows(
    nobootstrap_conf_int,
    bootstrap_conf_int
  ) %>%
  group_by(mean_smd, tau, omega, cor_mu, cor_sd, weights, m, n_multiplier, steps, priors) %>%
  mutate(
    bootstrap_condition = any(!is.na(bootstrap_type))
  ) %>%
  ungroup() %>%
  mutate(
    estimator = case_match(estimator, 'CML-fallback' ~ "CML", .default = estimator),
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
    param == "beta",
    estimator %in% c("CML","ARGL"),
    bootstraps < 1999L,
    bootstrap_type == "two-stage",
    weights %in% c("0.05","0.20"),
    N_factor == "Typical",
    het_ratio == 0.5,
    CI_type %in% c("percentile","basic","BCa")
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
    CI_type %in% c("percentile","basic","BCa")
  ) %>%
  mutate(
    bootstraps = if_else(is.na(bootstraps), 1999L, bootstraps),
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
  )

big_B_CML <- filter(big_B_bootstraps, estimator == "CML")
big_B_ARGL <- filter(big_B_bootstraps, estimator == "ARGL")



bootstraps <- unique(results_ci$bootstraps)[-1]

# boot_real %>%
#   filter(estimator == "CML") %>%
# ggplot() + 
#   aes(bootstraps, coverage, color = CI_type) + 
#   geom_hline(yintercept = 0.95, linetype = "dashed") + 
#   geom_point() + 
#   geom_smooth(method = "lm", formula = y ~ x, fullrange = TRUE, se = FALSE) + 
#   geom_pointrange(
#     data = big_B_CML,
#     aes(ymin = cover_lo, ymax = cover_hi),
#     shape = "square"
#   ) + 
#   facet_grid(tau ~ weights, labeller = "label_both") + 
#   scale_x_continuous(transform = "reciprocal", breaks = bootstraps) + 
#   theme_minimal() + 
#   labs(x = "B", y = "Coverage rate", color = "")
