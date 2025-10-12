library(tidyverse)

Lp_dist <- function(mu = 0, precision = 1, L = 2, ...) {
  Const <- gamma(1 / L) * 2 / (L * precision^(1 / L))
  f <- function(x) exp(-precision * abs(x - mu)^L) / Const
  g <- function(x) sapply(x, \(y) integrate(f, lower = -Inf, upper = y, ...)$value)
  
  return(list(pdf = f, cdf = g))
}

#-------------------------------------------------------------------------------
# Prior for beta (mean effect)

beta_default <- Lp_dist(mu = 0, precision = 1 / 2, L = 2)
beta_whoops <- Lp_dist(mu = 0, precision = 2, L = 2)
beta_weak <- Lp_dist(mu = 0, precision = 1 / 16, L = 4)

x <- seq(-4, 4, length.out = 500)

tibble(
  beta = x,
  default = beta_default$pdf(x),
  whoops = beta_whoops$pdf(x),
  weak = beta_weak$pdf(x)
) %>%
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_area(aes(beta, default), fill = "red", alpha = 0.5) + 
  geom_area(aes(beta, weak), fill = "purple", alpha = 0.5) + 
  scale_x_continuous(expand = expansion(0,0)) + 
  labs(
    x = expression(beta), 
    y = expression(p(beta))
  )

x_bench <- c(0.25, 0.5, 0.8, seq(1, 4, 0.5))

tibble(
  beta = x_bench,
  d_default = beta_default$pdf(x_bench) / beta_default$pdf(0),
  d_weak = beta_weak$pdf(x_bench) / beta_weak$pdf(0),
  c_default = 2 * (beta_default$cdf(x_bench) - 0.5),
  c_weak = 2 * (beta_weak$cdf(x_bench) - 0.5),
)

#-------------------------------------------------------------------------------
# Prior for zeta (log-selection)

z <- seq(log(0.01), log(16), length.out = 500)
lambda_bench <- c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.8, 1, 1.5, 2, 4, 8, 16)
z_bench <- log(lambda_bench)

zeta_default <- Lp_dist(mu = log(0.8), precision = 1 / 2, L = 2)
zeta_weak <- Lp_dist(mu = log(0.5), precision = 1 / 16, L = 4)

tibble(
  zeta = z,
  lambda = exp(z),
  default = zeta_default$pdf(z),
  weak = zeta_weak$pdf(z)
) %>%
  ggplot() + 
  geom_vline(xintercept = 1, linetype = "dashed") + 
  geom_area(aes(lambda, default), fill = "red", alpha = 0.5) + 
  geom_area(aes(lambda, weak), fill = "purple", alpha = 0.5) + 
  scale_x_continuous(
    transform = "log", 
    breaks = lambda_bench,
    expand = expansion(0,0)
  ) + 
  labs(
    x = expression(lambda), 
    y = expression(p(lambda))
  )

tibble(
  lambda = lambda_bench,
  zeta = z_bench,
  d_default = zeta_default$pdf(z_bench) / zeta_default$pdf(0),
  d_weak = zeta_weak$pdf(z_bench) / zeta_weak$pdf(0),
  c_default = zeta_default$cdf(z_bench),
  c_weak = zeta_weak$cdf(z_bench)
)

#-------------------------------------------------------------------------------
# Prior for gamma (log-heterogeneity)

loggamma_dist <- function(mode, alpha, ...) {
  lambda <- alpha / mode
  f <- function(x) dgamma(exp(x / 2), shape = alpha + 1, rate = lambda)
  g <- function(x) pgamma(exp(x / 2), shape = alpha, rate = lambda)
  return(list(pdf = f, cdf = g))
}


g <- seq(log(0.01^2), log(1^2), length.out = 500)
tau_bench <- c(0.01, 0.05, seq(0.1, 1, 0.1))
g_bench <- 2 * log(tau_bench)

gamma_default <- loggamma_dist(mode = 0.2, alpha = 1)
gamma_weak <- loggamma_dist(mode = 0.2, alpha = 0.5)


tibble(
  gamma = g,
  tau = exp(g / 2),
  default = gamma_default$pdf(g),
  weak = gamma_weak$pdf(g)
) %>%
  ggplot() + 
  geom_area(aes(tau, default), fill = "red", alpha = 0.5) + 
  geom_area(aes(tau, weak), fill = "purple", alpha = 0.5) + 
  scale_x_continuous(
    transform = "log", 
    breaks = tau_bench,
    expand = expansion(0,0)
  ) + 
  labs(
    x = expression(tau), 
    y = expression(p(tau))
  )

tibble(
  tau = tau_bench,
  gamma = g_bench,
  d_default = gamma_default$pdf(g_bench) / gamma_default$pdf(2 * log(0.2)),
  d_weak = gamma_weak$pdf(g_bench) / gamma_weak$pdf(2 * log(0.2)),
  c_default = gamma_default$cdf(g_bench),
  c_weak = gamma_weak$cdf(g_bench)
)
