source("research/step-function-paper/process-simulation-results.R")

# non scaled RMSE
ggplot(mu_graph_res) + 
  aes(x = weights, y = rmse, color = estimator, fill = estimator) +
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
    y = "Root Mean-Squared Error", 
    color = "Estimator",
    fill = "Estimator"
  ) + 
  theme_bw() +
  theme(legend.position = "top")

# coverage by rho

mu_graph_res_ci %>%
  filter(
    CI_type %in% c("large-sample")
  ) %>%
  ggplot(aes(x = J, y = coverage, color = estimator, fill = estimator)) +
  geom_boxplot(alpha = .5, coef = Inf) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 1.0)) + 
  scale_y_continuous(expand = expansion(c(0,0),c(0.02,0))) + 
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(
    cor_mu ~ mean_smd, 
    labeller = label_bquote(
      rows = rho == .(cor_mu),
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


# the boostrap one on only has rho = .8 condition
