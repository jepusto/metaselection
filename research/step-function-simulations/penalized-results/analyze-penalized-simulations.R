library(tidyverse)

res <- read_rds(file = "research/step-function-simulations/penalized-results/sim-step-function-penalized-results-no-bootstraps.rds")

res_beta_3PSM <- 
  res %>%
  filter(param == "beta", model == "3PSM", estimator %in% c("CML-fallback","ARGL")) %>%
  mutate(
    convergence = 100 * K_absolute / iterations,
    priors = factor(priors, levels = c("None","Weaker","Default")),
    estimator = case_match(estimator, 'CML-fallback' ~ "CML", .default = estimator)
  )

#-------------------------------------------------------------------------------
# Convergence

res_convergence <- 
  res %>%
  filter(param == "beta", model == "3PSM") %>%
  filter(estimator %in% c("ARGL","CML")) %>%
  mutate(
    convergence = 100 * K_absolute / iterations,
    priors = factor(priors, levels = c("None","Weaker","Default"))
  )

ggplot(res_convergence) + 
  aes(factor(weights), convergence, color = interaction(estimator, priors), fill = interaction(estimator, priors)) + 
  geom_boxplot() + 
  facet_grid(mean_smd ~ tau) + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "bias",
    color = "Estimator",
    fill = "Estimator"
  )

res_convergence %>%
  filter(estimator != "CML" | priors != "None") %>%
ggplot() +
  aes(factor(weights), convergence, color = interaction(estimator, priors), fill = interaction(estimator, priors)) + 
  geom_boxplot() + 
  facet_grid(tau ~ mean_smd, scales = "free_y") + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "bias",
    color = "Estimator",
    fill = "Estimator"
  )

res_convergence %>%
  filter(estimator == "ARGL") %>%
  ggplot() +
  aes(factor(weights), convergence, color = priors, fill = priors) + 
  geom_boxplot() + 
  facet_grid(tau ~ mean_smd, scales = "free_y") + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "bias",
    color = "Estimator",
    fill = "Estimator"
  )

#-------------------------------------------------------------------------------
# Bias

ggplot(res_beta_3PSM) + 
  aes(factor(weights), bias, color = interaction(estimator, priors), fill = interaction(estimator, priors)) + 
  geom_boxplot() + 
  facet_grid(mean_smd ~ tau) + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "bias",
    color = "Estimator",
    fill = "Estimator"
  )


res_beta_3PSM %>%
  filter(tau == 0.45) %>%
  ggplot() + 
  aes(factor(weights), bias, color = interaction(estimator, priors), fill = interaction(estimator, priors)) + 
  geom_boxplot() + 
  facet_grid(mean_smd ~ m) + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "bias",
    color = "Estimator",
    fill = "Estimator"
  )

res_beta_3PSM %>%
  filter(mean_smd == 0.00) %>%
  ggplot() + 
  aes(factor(weights), bias, color = interaction(estimator, priors), fill = interaction(estimator, priors)) + 
  geom_boxplot() + 
  facet_grid(tau ~ m) + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "bias",
    color = "Estimator",
    fill = "Estimator"
  )

res_beta_3PSM %>%
  filter(mean_smd == 0.00, tau == 0.45) %>%
  ggplot() + 
  aes(factor(weights), bias, color = interaction(estimator, priors), fill = interaction(estimator, priors)) + 
  geom_boxplot() + 
  facet_grid(n_multiplier ~ m, scales = "free_y") + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "bias",
    color = "Estimator",
    fill = "Estimator"
  )

#-------------------------------------------------------------------------------
# RMSE

ggplot(res_beta_3PSM) + 
  aes(factor(weights), rmse, color = interaction(estimator, priors), fill = interaction(estimator, priors)) + 
  geom_boxplot() + 
  facet_grid(mean_smd ~ tau) + 
  coord_cartesian(ylim = c(0,0.5)) + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "RMSE",
    color = "Estimator",
    fill = "Estimator"
  )

res_beta_3PSM %>%
  filter(mean_smd == 0.00, tau == 0.45) %>%
  ggplot() + 
  aes(factor(weights), rmse, color = interaction(estimator, priors), fill = interaction(estimator, priors)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0,0.5)) + 
  facet_grid(n_multiplier ~ m, scales = "free_y") + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "RMSE",
    color = "Estimator",
    fill = "Estimator"
  )


res_beta_3PSM_wide <- 
  res_beta_3PSM %>%
  select(mean_smd, tau, omega, cor_mu, weights, m, n_multiplier, estimator, priors, bias, rmse) %>%
  unite("estimator_prior", estimator, priors) %>%
  pivot_wider(names_from = estimator_prior, values_from = c(bias, rmse))

res_beta_3PSM_wide %>%
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

res_beta_3PSM_wide %>%
  ggplot() + 
  aes(bias_CML_Weaker, bias_CML_Default, color = factor(m)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point() + 
  facet_grid(mean_smd ~ tau, scales = "free") + 
  theme_light() + 
  labs(
    x = "Weaker", 
    y = "Default",
    color = "Sample size"
  )

res_beta_3PSM_wide %>%
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

res_beta_3PSM_wide %>%
  ggplot() + 
  aes(rmse_CML_None, rmse_CML_Weaker, color = factor(m)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point() + 
  facet_grid(mean_smd ~ tau, scales = "free") + 
  theme_light() + 
  labs(
    x = "Flat", 
    y = "Weaker",
    color = "Sample size"
  )

res_beta_3PSM_wide %>%
  ggplot() + 
  aes(rmse_CML_Weaker, rmse_CML_Default, color = factor(m)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point() + 
  facet_grid(mean_smd ~ tau, scales = "free") + 
  theme_light() + 
  labs(
    x = "Weaker", 
    y = "Default",
    color = "Sample size"
  )

res_beta_3PSM_wide %>%
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


res_beta_3PSM_wide %>%
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

res_beta_3PSM_wide %>%
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

#-------------------------------------------------------------------------------
# Coverage

ggplot(res_beta_3PSM) + 
  aes(factor(weights), coverage, color = interaction(estimator, priors), fill = interaction(estimator, priors)) + 
  geom_hline(yintercept = 0.95, linetype = "dashed") + 
  geom_boxplot() + 
  facet_grid(mean_smd ~ tau) + 
  theme_light() + 
  labs(
    x = "Selection weight", 
    y = "CI Coverage",
    color = "Estimator",
    fill = "Estimator"
  )
