
ggplot() + 
  scale_x_continuous(limits = c(-3,3)) + 
  geom_function(fun = dnorm, args = list(mean = 0, sd = 1, log = TRUE))

tau_mode <- 0.15
tau_sd <- 0.30
tau_lambda <- (tau_mode + sqrt(tau_mode^2 + 4 * tau_sd^2)) / (2 * tau_sd^2)
tau_alpha <- tau_lambda * tau_mode + 1
tau_alpha / tau_lambda
(tau_alpha - 1) / tau_lambda
sqrt(tau_alpha) / tau_lambda

ggplot() + 
  scale_x_continuous(limits = c(0,1)) + 
  geom_function(fun = dgamma, args = list(shape = tau_alpha, rate = tau_lambda, log = TRUE))

ggplot() + 
  scale_x_continuous(limits = c(0,3)) + 
  geom_function(fun = dgamma, args = list(shape = 2, rate = 1, log = TRUE))

ggplot() + 
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10, 0.20, 0.50, 1, 1.5, 2), limits = c(0.005,2), transform = "log") + 
  geom_function(fun = dlnorm, args = list(meanlog = 4, sdlog = 2, log = FALSE))

ggplot() + 
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10, 0.20, 0.50, 1, 1.5, 2), limits = c(0.005,2), transform = "log") + 
  geom_function(fun = dnorm, args = list(mean = 1, sd = 1, log = TRUE))
