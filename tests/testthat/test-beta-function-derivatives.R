# Make up some fake data

set.seed(20240207)
K <- 34
sei <- 1 / (2 * sqrt(rchisq(K, df = 8)))
yi <- 0.2 + rnorm(K, sd = .05) + rnorm(K, sd = sei)
X <- cbind(rep(1,K), rnorm(K))
colnames(X) <- c("intercept","X")
u_cats <- sample(LETTERS[1:4], size = K, replace = TRUE)
U <- model.matrix(~ 0 + u_cats)

test_that("beta_loglik() divides up the parameter vector appropriately.", {
  
  # intercept-only
  theta_A <- c(0.2, log(0.05), log(0.36), log(1.3))
  
  expect_equal(
    beta_loglik(theta = theta_A, yi = yi, sei = sei),
    beta_loglik(yi = yi, sei = sei, steps = c(.025, .975), 
                beta = theta_A[1], gamma = theta_A[2], zeta = theta_A[3:4])
  )
  
  # with predictors
  theta_B <- c(0.2, 0.04, log(seq(0.02, 0.08, length.out = 4)), log(0.36), log(1.3))
  
  expect_equal(
    beta_loglik(theta = theta_B, yi = yi, sei = sei, steps = c(.01, .99),
                X = X, U = U),
    beta_loglik(yi = yi, sei = sei, steps = c(.01, .99), 
                X = X, U = U,
                beta = theta_B[1:2], gamma = theta_B[3:6], zeta = theta_B[7:8])
  )
  
})

test_that("beta_score() divides up the parameter vector appropriately.", {
  
  # intercept-only
  theta_A <- c(0.2, log(0.05), log(0.36), log(1.3))

  A1 <- beta_score(theta = theta_A, yi = yi, sei = sei)
  A2 <- beta_score(yi = yi, sei = sei, steps = c(.025, .975),
                   beta = theta_A[1], gamma = theta_A[2], zeta = theta_A[3:4])
  expect_equal(A1, A2)
  expect_equal(length(theta_A), length(A2))
  
  # with predictors
  theta_B <- c(0.2, 0.04, log(seq(0.02, 0.08, length.out = 4)), log(0.36), log(1.3))
  
  B1 <- beta_score(theta = theta_B, yi = yi, sei = sei, steps = c(.01, .99),
                   X = X, U = U)
  B2 <- beta_score(yi = yi, sei = sei, steps= c(.01, .99), 
                   X = X, U = U, 
                   beta = theta_B[1:2], gamma = theta_B[3:6], zeta = theta_B[7:8])
  expect_equal(B1, B2)
  expect_equal(length(theta_B), length(B2))
  
})


test_that("beta_hessian() divides up the parameter vector appropriately.", {
  
  intcpt <- matrix(1, nrow = K, ncol = 1)
  colnames(intcpt) <- "intrcpt"
  
  # intercept-only
  theta_A <- c(0.2, log(0.05), log(0.36), log(1.3))

  A1 <- beta_hessian(theta = theta_A, yi = yi, sei = sei, steps = c(.025, .975))
  A2 <- beta_hessian(yi = yi, sei = sei, steps = c(.025, .975), 
                     beta = theta_A[1], gamma = theta_A[2], zeta = theta_A[3:4])
  A3 <- beta_hessian(theta = theta_A, yi = yi, sei = sei, X = intcpt, U = intcpt)
  expect_equal(length(theta_A), nrow(A2))
  expect_equal(length(theta_A), ncol(A2))
  expect_equal(A1, A2)
  expect_equal(A1, A3, ignore_attr = TRUE)

  # with predictors for mu
  theta_B <- c(0.2, 0.04, log(0.15), log(0.36), log(1.3))
  
  B1 <- beta_hessian(theta = theta_B, yi = yi, sei = sei, steps = c(.01, .88),
                     X = X)
  B2 <- beta_hessian(yi = yi, sei = sei, steps = c(.01, .88), 
                     X = X,
                     beta = theta_B[1:2], gamma = theta_B[3], zeta = theta_B[4:5])
  expect_equal(length(theta_B), nrow(B2))
  expect_equal(length(theta_B), ncol(B2))
  expect_equal(B1, B2)
  
  # with predictors for mu and tau
  theta_C <- c(0.2, 0.04, log(seq(0.02, 0.08, length.out = 4)), log(0.36), log(1.3))
  
  C1 <- beta_hessian(theta = theta_C, yi = yi, sei = sei, steps = c(.01, .88),
                     X = X, U = U)
  C2 <- beta_hessian(yi = yi, sei = sei, steps = c(.01, .88), 
                     X = X, U = U,
                     beta = theta_C[1:2], gamma = theta_C[3:6], zeta = theta_C[7:8])
  expect_equal(length(theta_C), nrow(C2))
  expect_equal(length(theta_C), ncol(C2))
  expect_equal(C1, C2)
  
})

test_that("E_Y_f and E_Y_f_vec work.", {
  
  integral_val <- function(sei, mu, eta, lambda, alpha, integrand) {
    bounds <- sei * qnorm(alpha, lower.tail = FALSE)
    res <- integrate(integrand, sei = sei, mu = mu, eta = eta, lambda = lambda, lower = bounds[2], upper = bounds[1])
    res$value  
  }
  
  eta <- 0.2^2 + sei^2
  
  integrand_A <- function(Y, sei, mu, eta, lambda) {
    A <- (pnorm(-Y/sei)^lambda[1]) * (pnorm(Y/sei)^lambda[2])
    B <- dnorm((Y - mu) / sqrt(eta)) / sqrt(eta)
    C <- (Y - mu) / sqrt(eta)
    return(A * B * C)
  }
  
  A0 <- mapply(
    integral_val, 
    sei = sei, eta = eta, 
    MoreArgs = list(integrand = integrand_A, mu = 0.2, lambda = c(0, 0), alpha = c(.01, .96))
  )
  
  A1 <- mapply(
    E_Y_f, sei = sei, eta = eta, 
    MoreArgs = list(f_exp = "(Y - mu) / sqrt(eta)", mu = 0.2, lambda = c(0, 0), alpha = c(.01, .96))
  )
  
  A2 <- E_Y_f_vec(
    f_exp = "(Y - mu) / sqrt(eta)", 
    sei = sei, mu = 0.2, eta = eta, 
    lambda = c(0,0), alpha = c(.01, .96)
  )
  
  expect_equal(A0, A1)
  expect_equal(A0, A2)
  
  
  integrand_B <- function(Y, sei, mu, eta, lambda) {
    A <- (pnorm(-Y/sei)^lambda[1]) * (pnorm(Y/sei)^lambda[2])
    B <- dnorm((Y - mu) / sqrt(eta)) / sqrt(eta)
    C <- ((Y - mu)^2 / eta - 1)
    return(A * B * C)
  }
  
  B0 <- mapply(
    integral_val, 
    sei = sei, eta = eta, 
    MoreArgs = list(integrand = integrand_B, mu = 0.2, lambda = c(0, 0), alpha = c(.01, .96))
  )
  
  B1 <- mapply(
    E_Y_f, sei = sei, eta = eta, 
    MoreArgs = list(f_exp = "((Y - mu)^2 / eta - 1)", mu = 0.2, lambda = c(0, 0), alpha = c(.01, .96))
  )
  
  B2 <- E_Y_f_vec(
    f_exp = "((Y - mu)^2 / eta - 1)", 
    sei = sei, mu = 0.2, eta = eta, 
    lambda = c(0,0), alpha = c(.01, .96)
  )
  
  expect_equal(B0, B1)
  expect_equal(B0, B2)
})

test_that("beta_score() is an unbiased estimating equation for models with no covariates.", {
  
  skip_if_not_installed("DescTools")
  verbose <- FALSE

  check_beta_score_hessian_bias(
    mean_smd = 0.3, 
    tau = 0.1, 
    m = 1000, 
    mean_N = 120,
    m_multiplier = 1,
    verbose = verbose,
    fit_mod = verbose,
    seed = 20240209
  )

  check_beta_score_hessian_bias(
    mean_smd = 0.1, 
    tau = 0.15, 
    m = 1000, 
    steps = c(.01, .99), 
    lambdas = c(1, 5),
    m_multiplier = 5,
    verbose = verbose,
    fit_mod = verbose,
    seed = 202404020
  )
  
  check_beta_score_hessian_bias(
    mean_smd = 0.0, 
    tau = 0.01, 
    m = 1000, 
    steps = c(.01, .94), 
    lambdas = c(1, 1),
    m_multiplier = 1,
    verbose = verbose,
    fit_mod = verbose,
    seed = 20240211
  )

  # this one isn't working  
  check_beta_score_hessian_bias(
    mean_smd = 0.02, 
    tau = 0.25, 
    m = 1000, 
    mean_N = 120,
    steps = c(.02, .98),
    lambdas = c(0.8, 1.9),
    m_multiplier = 2,
    verbose = verbose,
    fit_mod = verbose,
    seed = 20240419
  )
  
})

test_that("beta_score() is an unbiased estimating equation for models with covariates.", {
  
  skip_if_not_installed("DescTools")
  verbose <- FALSE

  check_beta_score_hessian_bias(
    mean_smd = c(0.3, 0.6), 
    tau = 0.1, 
    m = 1000, 
    mean_N = 120,
    m_multiplier = 1,
    verbose = verbose,
    seed = 20240420,
    hessian_threshold = 2
  )
  
  check_beta_score_hessian_bias(
    mean_smd = c(0.1, 0.26), 
    tau = 0.15, 
    m = 1000, 
    steps = c(.01, .99), 
    lambdas = c(1, 5),
    m_multiplier = 5,
    verbose = verbose,
    seed = 20240420
  )
  
  check_beta_score_hessian_bias(
    mean_smd = c(0.0, 0.2, 0.3), 
    tau = 0.3, 
    m = 1000, 
    steps = c(.01, .99), 
    lambdas = c(1, 1),
    m_multiplier = 1,
    verbose = verbose,
    seed = 20240421
  )
  
  check_beta_score_hessian_bias(
    mean_smd = c(0.0, 0.02, 0.04, 0.05), 
    tau = 0.5, 
    m = 1000, 
    steps = c(.05,0.5),
    lambdas = c(0.8, 0.9),
    m_multiplier = 5,
    verbose = verbose,
    seed = 20240212
  )
  
})


test_that("beta_score and beta_hessian agree with numerical derivatives.", {
  
  set.seed(20240418)
  
  dat <- r_meta(
    mean_smd = 0.2, 
    tau = 0.2, 
    omega = 0, 
    m = 50, 
    cor_mu = 0, 
    cor_sd = 0.01, 
    censor_fun = step_fun(cut_vals = c(.025, .50), weights = c(0.5, 0.2)), 
    n_ES_sim = n_ES_param(40, 1) 
  )
  
  beta_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "beta",
    steps = c(.025, .500),
    estimator = "ML",
    crit = qnorm(0.92)
  )
  beta_derivs$selmod_fit$est
  beta_derivs$score_diff_over_range
  expect_lt(max(beta_derivs$score_diff_over_range), 1e-3)
  round(beta_derivs$hess_diff_over_range, 5)
  expect_lt(max(beta_derivs$hess_diff_over_range), 1e-3)
  
  beta_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "beta",
    steps = c(.025, .975),
    estimator = "ML",
    params = c(1L, 3L, 4L)
  )
  beta_derivs$selmod_fit$est
  beta_derivs$score_diff_over_range
  expect_lt(max(beta_derivs$score_diff_over_range), 5e-4)
  round(beta_derivs$hess_diff_over_range, 5)
  expect_lt(max(beta_derivs$hess_diff_over_range), 5e-4)
  
  set.seed(20240418)
  
  dat <- r_meta(
    mean_smd = 0.2, 
    tau = 0.3, 
    omega = 0, 
    m = 50, 
    cor_mu = 0, 
    cor_sd = 0.01, 
    censor_fun = beta_fun(delta_1 = 1.3, delta_2 = 0.7, trunc_1 = .025, trunc_2 = .5), 
    n_ES_sim = n_ES_param(40, 1) 
  )
  
  beta_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "beta",
    steps = c(.025, .500),
    estimator = "ML"
  )
  beta_derivs$selmod_fit$est
  beta_derivs$score_diff_over_range
  expect_lt(max(beta_derivs$score_diff_over_range), 1e-3)
  round(beta_derivs$hess_diff_over_range, 5)
  expect_lt(max(beta_derivs$hess_diff_over_range), 5e-3)
  
  
  beta_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "beta",
    steps = c(.025, .975),
    estimator = "ML",
    crit = c(2, 1, 1, 0.5)
  )
  
  beta_derivs$selmod_fit$est
  beta_derivs$score_diff_over_range
  expect_lt(max(beta_derivs$score_diff_over_range), 1e-3)
  round(beta_derivs$hess_diff_over_range, 5)
  expect_lt(max(beta_derivs$hess_diff_over_range), 1e-3)
  
  
  # library(tidyverse)
  # ggplot(beta_derivs$data, aes(x = val)) +
  #   facet_wrap(~ param, scales = "free") +
  #   geom_hline(yintercept = 0) +
  #   geom_line(aes(y = d_num), color = "blue", linewidth = 2) +
  #   geom_line(aes(y = score), color = "green", linewidth = 1) +
  #   theme_minimal()
  # 
  # beta_derivs_hessian_comparison <-
  #   beta_derivs$data %>%
  #   select(param, val, starts_with("h_num"), starts_with("hess")) %>%
  #   pivot_longer(c(starts_with("h_num"), starts_with("hess")),
  #                names_to = "hp", values_to = "hess") %>%
  #   separate(hp, c("version", "hdim"), sep = "\\.") %>%
  #   group_by(param, val, version) %>%
  #   mutate(
  #     hdim = colnames(beta_derivs$hess_diff_over_range),
  #     param = paste(param, hdim, sep = "-")
  #   ) %>%
  #   pivot_wider(names_from = version, values_from = hess)
  # 
  # ggplot(beta_derivs_hessian_comparison, aes(x = val)) +
  #   facet_wrap(~  param, scales = "free") +
  #   geom_line(aes(y = h_num), color = "blue", linewidth = 2) +
  #   geom_line(aes(y = hess), color = "green", linewidth = 1) +
  #   theme_minimal()
  
})
