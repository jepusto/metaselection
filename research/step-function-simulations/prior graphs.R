library(tidyverse)


p_beta <- function(x, beta_mean = 0, beta_precision = 1 / 2, beta_L = 2) {
  -beta_precision * abs(x - beta_mean)^beta_L
}

x <- seq(-6, 6, length.out = 500)

tibble(
  beta = x,
  log_sq = p_beta(x),
  log_cube = p_beta(x, beta_precision = 1 / 4, beta_L = 3),
  log_quart = p_beta(x, beta_precision = 1 / 8, beta_L = 4)
) %>%
  ggplot() + 
  geom_line(aes(beta, log_sq)) + 
  geom_line(aes(beta, log_cube), color = "red") + 
  geom_line(aes(beta, log_quart), color = "purple") + 
  scale_x_continuous(expand = expansion(0,0)) + 
  scale_y_continuous(expand = expansion(0,0)) + 
  coord_cartesian(xlim = c(-4, 4), ylim = c(-18,0)) + 
  labs(
    x = expression(beta), 
    y = expression(log-prior(beta))
  )

p_gamma <- function(x, tau_mode = 0.20, gamma_alpha = 1) {
  gamma_lambda <- gamma_alpha / tau_mode
  gamma_alpha * x / 2 - gamma_lambda * exp(x / 2)
}

tibble(
  gamma = seq(log(0.01^2), log(0.50^2), length.out = 500),
) %>%
  mutate(
    tau = exp(gamma / 2),
    log_prior = p_gamma(gamma),
    log_prior = log_prior - max(log_prior)
  ) %>%
  ggplot() + 
  aes(tau, log_prior) + 
  geom_line() + 
  scale_x_continuous(
    transform = "log", 
    breaks = c(0.01, seq(0.05, 0.50, 0.05)),
    expand = expansion(0,0)
  ) + 
  scale_y_continuous(expand = expansion(0,0)) + 
  labs(
    x = expression(tau), 
    y = expression(log-prior(gamma))
  )


p_zeta <- function(x, lambda_mode = 0.8, zeta_precision = 1 / 2, zeta_L = 2) {
  zeta_mean <- log(lambda_mode)
  -zeta_precision * abs(x - zeta_mean)^zeta_L
}

tibble(
  zeta = seq(log(0.01), log(16), length.out = 500),
) %>%
  mutate(
    lambda = exp(zeta),
    log_sq = p_zeta(zeta),
    log_cube = p_zeta(zeta, zeta_precision = 1 / 7.377, zeta_L = 3),
    log_quart = p_zeta(zeta, zeta_precision = 1 / 27, zeta_L = 4),
    across(starts_with("log_"), ~ .x - max(.x))
  ) %>%
  ggplot() + 
  geom_line(aes(lambda, log_sq)) + 
  geom_line(aes(lambda, log_cube), color = "red") + 
  geom_line(aes(lambda, log_quart), color = "purple") + 
  scale_x_continuous(
    transform = "log", 
    breaks = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1, 1.5, 2, 4, 8, 16),
    expand = expansion(0,0)
  ) + 
  scale_y_continuous(expand = expansion(0,0)) + 
  coord_cartesian(ylim = c(-20,0)) + 
  labs(
    x = expression(lambda), 
    y = expression(log-prior(zeta))
  )
