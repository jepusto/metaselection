source(here::here("research","step-function-paper","process-simulation-results.R"))


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Bias ----

mu_graph_res %>%
  filter(psi_fac == "Independent") %>%
  performance_plot(measure = bias, label = "Bias")


mu_graph_res %>%
  filter(psi_fac == "Dependent") %>%
  performance_plot(measure = bias, label = "Bias")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# RMSE ----

mu_graph_res %>%
  filter(psi_fac == "Independent") %>%
  performance_plot(measure = rmse, label = "Root Mean-Squared Error")

mu_graph_res %>%
  filter(psi_fac == "Dependent") %>%
  performance_plot(measure = rmse, label = "Root Mean-Squared Error")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# RMSE comparisons ----

# PML vs. ARGL

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


mu_wide_res %>%
  RMSE_comparison_plot("PML","ARGL")


mu_wide_res %>%
  filter(psi_fac == "Independent") %>%
  RMSE_comparison_plot("PML","ARGL")

mu_wide_res %>%
  filter(psi_fac == "Dependent") %>%
  RMSE_comparison_plot("PML","ARGL")

# PML vs. CHE-ISCW

mu_wide_res %>%
  RMSE_comparison_plot("PML","CHE-ISCW")

mu_wide_res %>%
  filter(psi_fac == "Independent") %>%
  RMSE_comparison_plot("PML","CHE-ISCW")

mu_wide_res %>%
  filter(psi_fac == "Dependent") %>%
  RMSE_comparison_plot("PML","CHE-ISCW")

# ARGL vs. CHE-ISCW

mu_wide_res %>%
  RMSE_comparison_plot("ARGL","CHE-ISCW")

mu_wide_res %>%
  filter(psi_fac == "Independent") %>%
  RMSE_comparison_plot("ARGL","CHE-ISCW")

mu_wide_res %>%
  filter(psi_fac == "Dependent") %>%
  RMSE_comparison_plot("ARGL","CHE-ISCW")

# PML vs. CHE-ISCW

mu_wide_res %>%
  RMSE_comparison_plot("PML","PET/PEESE")

mu_wide_res %>%
  filter(psi_fac == "Independent") %>%
  RMSE_comparison_plot("PML","PET/PEESE")

mu_wide_res %>%
  filter(psi_fac == "Dependent") %>%
  RMSE_comparison_plot("PML","PET/PEESE")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Large-sample CI coverage ----

mu_graph_res_ci %>%
  filter(
    CI_type %in% c("large-sample"), 
    psi_fac == "Independent"
  ) %>%
  ggplot(aes(x = J, y = coverage, color = estimator, fill = estimator)) +
  geom_boxplot(alpha = .5, coef = Inf) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 1.0)) + 
  scale_y_continuous(expand = expansion(c(0,0),c(0.02,0))) + 
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(
    tau ~ mean_smd, 
    labeller = label_bquote(
      rows = tau[B] == .(tau),
      cols = mu == .(mean_smd)
    ),
    scales = "free_y"
  ) +
  labs(
    x = "Number of studies (J)", 
    y = "Coverage rate", 
    color = "Estimator",
    fill = "Estimator"
  ) + 
  theme_bw() +
  theme(legend.position = "top")

mu_graph_res_ci %>%
  filter(
    CI_type %in% c("large-sample"), 
    psi_fac == "Dependent"
  ) %>%
  ggplot(aes(x = J, y = coverage, color = estimator, fill = estimator)) +
  geom_boxplot(alpha = .5, coef = Inf) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 1.0)) + 
  scale_y_continuous(expand = expansion(c(0,0),c(0.02,0))) + 
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(
    tau ~ mean_smd, 
    labeller = label_bquote(
      rows = tau[B] == .(tau),
      cols = mu == .(mean_smd)
    ),
    scales = "free_y"
  ) +
  labs(
    x = "Number of studies (J)", 
    y = "Coverage rate", 
    color = "Estimator",
    fill = "Estimator"
  ) + 
  theme_bw() +
  theme(legend.position = "top")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Bootstrap CI coverage ----

mu_graph_res_ci %>%
  filter(
    bootstrap_condition,
    estimator %in% c("PML"),
    CI_type %in% c("large-sample","percentile"),
    psi_fac == "Independent"
  ) %>%
  ggplot(aes(x = J, y = coverage, color = CI_boot_method, fill = CI_boot_method)) +
  geom_boxplot(alpha = .5, coef = Inf) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_y_continuous(limits = c(0.65, 1.0), breaks = seq(0.65,0.95,0.05), expand = expansion(0,0)) +
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


mu_graph_res_ci %>%
  filter(
    bootstrap_condition,
    estimator %in% c("PML"),
    CI_type %in% c("large-sample","percentile"),
    psi_fac == "Dependent"
  ) %>%
  ggplot(aes(x = J, y = coverage, color = CI_boot_method, fill = CI_boot_method)) +
  geom_boxplot(alpha = .5, coef = Inf) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_y_continuous(limits = c(0.65, 1.0), breaks = seq(0.65,0.95,0.05), expand = expansion(0,0)) +
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
