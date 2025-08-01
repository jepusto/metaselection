source("process-simulation-results.R")

results_ci %>%
  filter(
    mean_smd %in% c(0.2), 
    m == 15, 
    bootstrap_condition == "bootstrap", 
    param == "beta",
    estimator == "CML",
    bootstrap_type == "two-stage",
    bootstraps == 1999,
    weights %in% c("0.05","0.20"),
    N_factor == "Typical"
  ) %>% 
  ggplot() + 
  aes(CI_type, coverage, color = het_ratio, shape == het_ratio) + 
  geom_point() + 
  facet_grid(weights ~ tau) + 
  theme_minimal()

boot_real <- 
  results_ci %>%
  filter(
    mean_smd %in% c(0.2), 
    m == 15, 
    bootstrap_condition == "bootstrap", 
    param == "beta",
    estimator == "CML",
    bootstraps < 1999L,
    bootstrap_type == "two-stage",
    weights %in% c("0.05","0.20"),
    N_factor == "Typical",
    het_ratio == 0.5
  )

boot_extrapolations <- 
  results_ci %>%
  filter(
    mean_smd %in% c(0.2), 
    m == 15, 
    bootstrap_condition == "bootstrap", 
    param == "beta",
    estimator == "CML",
    bootstraps == 1999L,
    bootstrap_type == "two-stage",
    weights %in% c("0.05","0.20"),
    N_factor == "Typical",
    het_ratio == 0.5
  )


bootstraps <- unique(results_ci$bootstraps)[-1]

ggplot(boot_real) + 
  aes(bootstraps, coverage, color = CI_type, shape == CI_type) + 
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ x, fullrange = TRUE) + 
  geom_point(data = boot_extrapolations, shape = "square") + 
  facet_grid(weights ~ tau) + 
  scale_x_continuous(transform = "reciprocal", breaks = bootstraps) + 
  theme_minimal()
