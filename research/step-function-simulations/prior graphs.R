library(tidyverse)


p_beta <- function(x, beta_mean = 0, beta_sd = 1) {
  -0.5 * (x - beta_mean)^2 / beta_sd^2
}

x <- seq(-2, 2, length.out = 500)

tibble(
  beta = x,
  log_prior = p_beta(x) - max(p_beta(x))
) %>%
  ggplot() + 
  aes(beta, log_prior) + 
  geom_line() + 
  scale_x_continuous(expand = expansion(0,0)) + 
  scale_y_continuous(expand = expansion(0,0)) + 
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


p_zeta <- function(x, lambda_mode = 0.8, zeta_alpha = 1) {
  zeta_lambda <- zeta_alpha / lambda_mode
  zeta_alpha * x - zeta_lambda * exp(x)
}

tibble(
  zeta = seq(log(0.01), log(4), length.out = 500),
) %>%
  mutate(
    lambda = exp(zeta),
    log_prior = p_zeta(zeta),
    log_prior2 = - 0.5 * zeta^2 + exp(zeta) * log(0.8) / 0.8,
    across(starts_with("log_prior"), ~ .x - max(.x))
  ) %>%
  ggplot() + 
  geom_line(aes(lambda, log_prior)) + 
  geom_line(aes(lambda, log_prior2), color = "green") + 
  scale_x_continuous(
    transform = "log", 
    breaks = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1, 1.5, 2, 4),
    expand = expansion(0,0)
  ) + 
  scale_y_continuous(expand = expansion(0,0)) + 
  labs(
    x = expression(lambda), 
    y = expression(log-prior(zeta))
  )
