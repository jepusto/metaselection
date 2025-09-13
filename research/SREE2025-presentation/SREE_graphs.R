library(tidyverse)

setwd("C:/GitHub/metaselection/research/SREE2025-presentation")

source("process-simulation-results.R")

## Beta-Density Selection Model Compared to CHE-ISCW and PET/PEESE ##

# Bias - all
bias_main_all <- ggplot(mu_graph_res_main) + 
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

ggsave("bias_main_all.png", plot = bias_main_all, width = 10, height = 6, dpi = 300)

# Bias - mu = 0, .2, .8 and tau = .15, .45
bias_main_sub <- ggplot(mu_graph_res_main %>% filter(mean_smd %in% c(0, 0.2, 0.8), tau %in% c(0.15, 0.45))) + 
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

ggsave("bias_main_sub.png", plot = bias_main_sub, width = 10, height = 6, dpi = 300)

# Scaled RMSE - all
rmse_main_all <- ggplot(mu_graph_res_main) + 
  aes(x = selection_strength, y = scrmse, color = method, fill = method) +
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
    y = "Scaled Root Mean-Squared Error", 
    color = "",
    fill = ""
  ) + 
  theme_bw() +
  theme(legend.position = "top")

ggsave("rmse_main_all.png", plot = rmse_main_all, width = 10, height = 6, dpi = 300)

# Scaled RMSE - mu = 0, .2, .8 and tau = .15, .45
rmse_main_sub <- ggplot(mu_graph_res_main %>% filter(mean_smd %in% c(0, 0.2, 0.8), tau %in% c(0.15, 0.45))) + 
  aes(x = selection_strength, y = scrmse, color = method, fill = method) +
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
    y = "Scaled Root Mean-Squared Error", 
    color = "",
    fill = ""
  ) + 
  theme_bw() +
  theme(legend.position = "top")

ggsave("rmse_main_sub.png", plot = rmse_main_sub, width = 10, height = 6, dpi = 300)

# Coverage - all
coverage_main_all <- ggplot(mu_graph_res_ci_main %>% filter(CI_type %in% c("large-sample"))) +
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

ggsave("coverage_main_all.png", plot = coverage_main_all, width = 10, height = 6, dpi = 300)

# Coverage - mu = 0, .2, .8 and tau = .15, .45
coverage_main_sub <- ggplot(mu_graph_res_ci_main %>% filter(CI_type %in% c("large-sample"), mean_smd %in% c(0, 0.2, 0.8), tau %in% c(0.15, 0.45))) +
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

ggsave("coverage_main_sub.png", plot = coverage_main_sub, width = 10, height = 6, dpi = 300)

# Coverage - bootstrap
coverage_main_boot_all <- ggplot(mu_graph_res_ci_main %>% filter(bootstrap_condition == "bootstrap", estimator %in% c("CML"))) +
  aes(x = J, y = coverage, color = CI_boot_method, fill = CI_boot_method) +
  geom_boxplot(alpha = .5, coef = Inf) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  #scale_y_continuous(limits = c(0.55, 1.0), breaks = seq(0.55,1.0,0.05), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(0.55, 1.0)) +
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
    color = "Bootstrap Method",
    fill = "Bootstrap Method"
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave("coverage_main_boot_all.png", plot = coverage_main_boot_all, width = 10, height = 6, dpi = 300)

# Coverage - bootstrap, mu = 0, .2 and tau = .15, .45
coverage_main_boot_sub <- ggplot(mu_graph_res_ci_main %>% filter(bootstrap_condition == "bootstrap", estimator %in% c("CML"), CI_type %in% c("large-sample", "percentile"), mean_smd %in% c(0, 0.2), tau %in% c(0.15, 0.45))) +
  aes(x = J, y = coverage, color = CI_type, fill = CI_type) +
  geom_boxplot(alpha = .5, coef = Inf) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  #scale_y_continuous(limits = c(0.55, 1.0), breaks = seq(0.55,1.0,0.05), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(0.55, 1.0)) +
  scale_color_brewer(
    palette = "Dark2",
    labels = c("large-sample" = "RVE large-sample", "percentile" = "Cluster bootstrapping")
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    labels = c("large-sample" = "RVE large-sample", "percentile" = "Cluster bootstrapping")
  ) +
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

ggsave("coverage_main_boot_sub.png", plot = coverage_main_boot_sub, width = 10, height = 6, dpi = 300)


## Beta-Density Compared to One-Step and Two-Step Selection Models ##

# Bias - all
bias_miss_all <- ggplot(mu_graph_res_miss) + 
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

ggsave("bias_miss_all.png", plot = bias_miss_all, width = 10, height = 6, dpi = 300)

# Bias - mu = 0, .2, .8 and tau = .15, .45
bias_miss_sub <- ggplot(mu_graph_res_main %>% filter(mean_smd %in% c(0, 0.2, 0.8), tau %in% c(0.15, 0.45))) + 
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

ggsave("bias_miss_sub.png", plot = bias_miss_sub, width = 10, height = 6, dpi = 300)

# Scaled RMSE - all
rmse_miss_all <- ggplot(mu_graph_res_miss) + 
  aes(x = selection_strength, y = scrmse, color = method, fill = method) +
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
    y = "Scaled Root Mean-Squared Error", 
    color = "",
    fill = ""
  ) + 
  theme_bw() +
  theme(legend.position = "top")

ggsave("rmse_miss_all.png", plot = rmse_miss_all, width = 10, height = 6, dpi = 300)

# Scaled RMSE - mu = 0, .2, .8 and tau = .15, .45
rmse_miss_sub <- ggplot(mu_graph_res_main %>% filter(mean_smd %in% c(0, 0.2, 0.8), tau %in% c(0.15, 0.45))) + 
  aes(x = selection_strength, y = scrmse, color = method, fill = method) +
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
    y = "Scaled Root Mean-Squared Error", 
    color = "",
    fill = ""
  ) + 
  theme_bw() +
  theme(legend.position = "top")

ggsave("rmse_miss_sub.png", plot = rmse_miss_sub, width = 10, height = 6, dpi = 300)

# Coverage - all
coverage_miss_all <- ggplot(mu_graph_res_ci_miss %>% filter(CI_type %in% c("large-sample"))) +
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

ggsave("coverage_miss_all.png", plot = coverage_miss_all, width = 10, height = 6, dpi = 300)

# Coverage - mu = 0, .2, .8 and tau = .15, .45
coverage_miss_sub <- ggplot(mu_graph_res_ci_miss %>% filter(CI_type %in% c("large-sample"), mean_smd %in% c(0, 0.2, 0.8), tau %in% c(0.15, 0.45))) +
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

ggsave("coverage_miss_sub.png", plot = coverage_miss_sub, width = 10, height = 6, dpi = 300)
