# Make up some fake data

set.seed(20230817)
K <- 25
sei <- 1 / (2 * sqrt(rchisq(K, df = 8)))
yi <- 0.2 + rnorm(K, sd = .05) + rnorm(K, sd = sei)
X <- cbind(rep(1,K), rnorm(K))
colnames(X) <- c("intercept","X")
u_cats <- sample(LETTERS[1:4], size = K, replace = TRUE)
U <- model.matrix(~ 0 + u_cats)
z_cats <- sample(LETTERS[1:2], size = K, replace = TRUE)
Z <- model.matrix(~ 0 + z_cats)
z0_cats <- sample(LETTERS[1:3], size = K, replace = TRUE)
Z0 <- model.matrix(~ 0 + z0_cats)[,-1]

test_that("step_weighted_logpartlik() divides up the parameter vector appropriately.", {
  
  # intercept-only, 3PSM
  theta_A <- c(0.2, log(0.05), log(0.36))
  expect_equal(
    step_weighted_logpartlik(theta = theta_A, yi = yi, sei = sei, steps = .07),
    step_weighted_logpartlik(yi = yi, sei = sei, steps = .07, 
                             beta = theta_A[1], gamma = theta_A[2], zeta = theta_A[3])
  )
  
  # intercept-only, 4PSM
  theta_E <- c(0.2, log(0.05), log(0.36), log(0.54))
  expect_equal(
    step_weighted_logpartlik(theta = theta_E, yi = yi, sei = sei, steps = c(.07,.50)),
    step_weighted_logpartlik(yi = yi, sei = sei, steps = c(.07,.50),
                             beta = theta_E[1], gamma = theta_E[2], zeta = theta_E[3:4])
  )
  
  # with predictors, one step
  theta_B <- c(0.2, 0.04, log(seq(0.02, 0.08, length.out = 4)), log(0.36), log(0.54))
  expect_equal(
    step_weighted_logpartlik(theta = theta_B, yi = yi, sei = sei, steps = .07,
                             X = X, U = U, Z = Z),
    step_weighted_logpartlik(yi = yi, sei = sei, steps = .07, 
                             X = X, U = U, Z = Z,
                             beta = theta_B[1:2], gamma = theta_B[3:6], zeta = theta_B[7:8])
  )
  
  # with predictors, multiple steps
  theta_C <- c(0.2, 0.04, log(seq(0.02, 0.08, length.out = 4)), 
               log(seq(0.2, 0.5, length.out = 6)))
  expect_equal(
    step_weighted_logpartlik(theta = theta_C, yi = yi, sei = sei, steps = c(.02,0.08,0.20),
                             X = X, U = U, Z = list(Z, Z, Z)),
    step_weighted_logpartlik(yi = yi, sei = sei, steps = c(.02,0.08,0.20), 
                             X = X, U = U, Z = list(Z, Z, Z),
                             beta = theta_C[1:2], gamma = theta_C[3:6], zeta = theta_C[7:12])
  )
  
  # with predictors for lambda0
  expect_equal(
    step_weighted_logpartlik(theta = theta_C, yi = yi, sei = sei, steps = c(.02,0.08),
                             X = X, U = U, Z0 = Z0, Z = list(Z, Z)),
    step_weighted_logpartlik(yi = yi, sei = sei, steps = c(.02,0.08), 
                             X = X, U = U, Z0 = Z0, Z = list(Z, Z),
                             beta = theta_C[1:2], gamma = theta_C[3:6], zeta0 = theta_C[7:8], zeta = theta_C[9:12])
  )
  
  theta_D <- c(0.2, log(0.05), log(0.36), log(1.2), 0)
  expect_equal(
    step_weighted_logpartlik(theta = theta_D, yi = yi, sei = sei, steps = .07, Z0 = Z0),
    step_weighted_logpartlik(yi = yi, sei = sei, steps = .07, Z0 = Z0, 
                             beta = theta_D[1], gamma = theta_D[2], zeta0 = theta_D[3:4], zeta = theta_D[5])
  )
  
})


test_that("step_hybrid_score() divides up the parameter vector appropriately.", {
  
  # intercept-only, 3PSM
  theta_A <- c(0.2, log(0.05), log(0.36))
  A1 <- step_hybrid_score(theta = theta_A, yi = yi, sei = sei, steps = .07)
  A2 <- step_hybrid_score(yi = yi, sei = sei, steps = .07, 
                   beta = theta_A[1], gamma = theta_A[2], zeta = theta_A[3])
  expect_equal(A1, A2)
  expect_equal(length(theta_A), length(A2))
  
  # intercept-only, 4PSM
  theta_E <- c(0.2, log(0.05), log(0.36), log(0.54))
  E1 <- step_hybrid_score(theta = theta_E, yi = yi, sei = sei, steps = c(.07,.50))
  E2 <- step_hybrid_score(yi = yi, sei = sei, steps = c(.07,.50), 
                   beta = theta_E[1], gamma = theta_E[2], zeta = theta_E[3:4])
  expect_equal(E1, E2)
  expect_equal(length(theta_E), length(E2))
  
  # with predictors, one step
  theta_B <- c(0.2, 0.04, log(seq(0.02, 0.08, length.out = 4)), log(0.36), log(0.54))
  B1 <- step_hybrid_score(theta = theta_B, yi = yi, sei = sei, steps = .07,
                   X = X, U = U, Z = Z)
  B2 <- step_hybrid_score(yi = yi, sei = sei, steps = .07, 
                   X = X, U = U, Z = Z,
                   beta = theta_B[1:2], gamma = theta_B[3:6], zeta = theta_B[7:8])
  expect_equal(B1, B2)
  expect_equal(length(theta_B), length(B2))
  
  # with predictors, multiple steps
  theta_C <- c(0.2, 0.04, log(seq(0.02, 0.08, length.out = 4)), 
               log(seq(0.2, 0.5, length.out = 6)))
  C1 <- step_hybrid_score(theta = theta_C, yi = yi, sei = sei, steps = c(.02,0.08,0.20),
                   X = X, U = U, Z = list(Z, Z, Z))
  C2 <- step_hybrid_score(yi = yi, sei = sei, steps = c(.02,0.08,0.20), 
                   X = X, U = U, Z = list(Z, Z, Z),
                   beta = theta_C[1:2], gamma = theta_C[3:6], zeta = theta_C[7:12])
  expect_equal(C1, C2)
  expect_equal(length(theta_C), length(C2))
  
  # with predictors for lambda0
  C3 <- step_hybrid_score(theta = theta_C, yi = yi, sei = sei, steps = c(.02,0.08),
                   X = X, U = U, Z0 = Z0, Z = list(Z, Z))
  C4 <- step_hybrid_score(yi = yi, sei = sei, steps = c(.02,0.08), 
                   X = X, U = U, Z0 = Z0, Z = list(Z, Z),
                   beta = theta_C[1:2], gamma = theta_C[3:6], 
                   zeta0 = theta_C[7:8], zeta = theta_C[9:12])
  expect_equal(C3, C4)
  expect_equal(length(theta_C), length(C4))
  
  
  theta_D <- c(0.2, log(0.05), log(0.36), log(1.2), 0)
  D1 <- step_hybrid_score(theta = theta_D, yi = yi, sei = sei, steps = .07, Z0 = Z0)
  D2 <- step_hybrid_score(yi = yi, sei = sei, steps = .07, Z0 = Z0, 
                   beta = theta_D[1], gamma = theta_D[2], 
                   zeta0 = theta_D[3:4], zeta = theta_D[5])
  expect_equal(D1, D2)
  expect_equal(length(theta_D), length(D2))
  
})

intcpt <- matrix(1, nrow = K, ncol = 1)
colnames(intcpt) <- "intrcpt"

test_that("step_hybrid_jacobian() divides up the parameter vector appropriately.", {
  
  theta_A <- c(0.2, log(0.05), log(0.36))
  A1 <- step_hybrid_jacobian(theta = theta_A, yi = yi, sei = sei, steps = .07)
  A2 <- step_hybrid_jacobian(yi = yi, sei = sei, steps = .07, 
                     beta = theta_A[1], gamma = theta_A[2], zeta = theta_A[3])
  A3 <- step_hybrid_jacobian(theta = theta_A, yi = yi, sei = sei, steps = .07, X = intcpt, U = intcpt)
  A4 <- step_hybrid_jacobian(theta = theta_A, yi = yi, sei = sei, steps = .07, Z = intcpt)
  A5 <- step_hybrid_jacobian(theta = theta_A, yi = yi, sei = sei, steps = .07, Z = list(intcpt))
  expect_equal(length(theta_A), nrow(A2))
  expect_equal(length(theta_A), ncol(A2))
  expect_equal(A1, A2)
  expect_equal(A1, A3, ignore_attr = TRUE)
  expect_equal(A1, A4, ignore_attr = TRUE)
  expect_equal(A1, A5, ignore_attr = TRUE)
  
  # intercept-only, 4PSM
  theta_E <- c(0.2, log(0.05), log(0.36), log(0.54))
  E1 <- step_hybrid_jacobian(theta = theta_E, yi = yi, sei = sei, steps = c(.07,.50))
  E2 <- step_hybrid_jacobian(yi = yi, sei = sei, steps = c(.07,.50), 
                     beta = theta_E[1], gamma = theta_E[2], zeta = theta_E[3:4])
  E3 <- step_hybrid_jacobian(theta = theta_E, yi = yi, sei = sei, steps = c(.07,.50), X = intcpt)
  E4 <- step_hybrid_jacobian(theta = theta_E, yi = yi, sei = sei, steps = c(.07,.50), Z = list(intcpt, intcpt))
  expect_equal(length(theta_E), nrow(E2))
  expect_equal(length(theta_E), ncol(E2))
  expect_equal(E1, E2)
  expect_equal(E1, E3, ignore_attr = TRUE)
  expect_equal(E1, E4, ignore_attr = TRUE)
  
  # with predictors, one step
  theta_B <- c(0.2, 0.04, log(seq(0.02, 0.08, length.out = 4)), log(0.36), log(0.54))
  B1 <- step_hybrid_jacobian(theta = theta_B, yi = yi, sei = sei, steps = .07,
                     X = X, U = U, Z = Z)
  B2 <- step_hybrid_jacobian(yi = yi, sei = sei, steps = .07, 
                     X = X, U = U, Z = Z,
                     beta = theta_B[1:2], gamma = theta_B[3:6], zeta = theta_B[7:8])
  expect_equal(length(theta_B), nrow(B2))
  expect_equal(length(theta_B), ncol(B2))
  expect_equal(B1, B2)
  
  # with predictors, multiple steps
  theta_C <- c(0.2, 0.04, log(seq(0.02, 0.08, length.out = 4)), 
               log(seq(0.2, 0.5, length.out = 6)))
  C1 <- step_hybrid_jacobian(theta = theta_C, yi = yi, sei = sei, steps = c(.02,0.08,0.20),
                    X = X, U = U, Z = list(Z, Z, Z))
  C2 <- step_hybrid_jacobian(yi = yi, sei = sei, steps = c(.02,0.08,0.20), 
                     X = X, U = U, Z = list(Z, Z, Z),
                     beta = theta_C[1:2], gamma = theta_C[3:6], zeta = theta_C[7:12])
  expect_equal(C1, C2)
  expect_equal(length(theta_C), nrow(C2))
  expect_equal(length(theta_C), ncol(C2))
  
  # with predictors for lambda0
  F1 <- step_hybrid_jacobian(theta = theta_C, yi = yi, sei = sei, steps = c(.02,0.08),
                     X = X, U = U, Z0 = Z0, Z = list(Z, Z))
  F2 <- step_hybrid_jacobian(yi = yi, sei = sei, steps = c(.02,0.08), 
                     X = X, U = U, Z0 = Z0, Z = list(Z, Z),
                     beta = theta_C[1:2], gamma = theta_C[3:6], 
                     zeta0 = theta_C[7:8], zeta = theta_C[9:12])
  expect_equal(length(theta_C), nrow(F2))
  expect_equal(length(theta_C), ncol(F2))
  expect_equal(F1, F2)

  theta_D <- c(0.2, log(0.05), log(0.36), log(1.2), 0)
  D1 <- step_hybrid_jacobian(theta = theta_D, yi = yi, sei = sei, steps = .07, Z0 = Z0)
  D2 <- step_hybrid_jacobian(yi = yi, sei = sei, steps = .07, Z0 = Z0, 
                     beta = theta_D[1], gamma = theta_D[2], 
                     zeta0 = theta_D[3:4], zeta = theta_D[5])
  D3 <- step_hybrid_jacobian(theta = theta_D, yi = yi, sei = sei, steps = .07, Z0 = Z0, Z = intcpt)
  expect_equal(length(theta_D), nrow(D2))
  expect_equal(length(theta_D), ncol(D2))
  expect_equal(D1, D2)
  expect_equal(D1, D3, ignore_attr = TRUE)
  
})


test_that("step_hybrid_score() is an unbiased estimating equation for models with no covariates.", {
  
  skip_if_not_installed("DescTools")
  verbose <- FALSE
  
  check_step_score_hessian_bias(
    mean_smd = 0.3, 
    tau = 0.1, 
    m = 1000, 
    steps = .025, 
    weights = 1,
    m_multiplier = 1,
    verbose = verbose,
    seed = 20230622,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )
  
  check_step_score_hessian_bias(
    mean_smd = 0.1, 
    tau = 0.15, 
    m = 1000, 
    steps = .025, 
    weights = 0.05,
    m_multiplier = 5,
    verbose = verbose,
    seed = 20230623,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )
  
  check_step_score_hessian_bias(
    mean_smd = 0.0, 
    tau = 0.01, 
    m = 1000, 
    steps = c(.025,0.5),
    weights = c(1, 1),
    m_multiplier = 1,
    verbose = verbose,
    seed = 20230624,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )
  
  check_step_score_hessian_bias(
    mean_smd = 0.03, 
    tau = 0.5, 
    m = 1000, 
    steps = c(.025,0.5),
    weights = c(0.1, 0.4),
    m_multiplier = 5,
    verbose = verbose,
    seed = 20230625,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )

})

test_that("step_hybrid_score() is an unbiased estimating equation for models with covariates.", {
  
  skip_if_not_installed("DescTools")
  verbose <- FALSE
  
  check_step_score_hessian_bias(
    mean_smd = seq(0, 0.4, 0.13), 
    tau = 0.2, 
    m = 1000, 
    steps = .025, 
    weights = 1,
    m_multiplier = 1,
    mean_N = 150,
    verbose = verbose,
    seed = 20230815,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )
  
  check_step_score_hessian_bias(
    mean_smd = c(0, 0.3), 
    tau = c(0.50,0.20), 
    m = 1000, 
    steps = .025, 
    weights = 0.20,
    m_multiplier = 5,
    mean_N = 100,
    verbose = verbose,
    seed = 20230727,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )
  
  check_step_score_hessian_bias(
    mean_smd = c(0.8, 0.3), 
    tau = c(0.50,0.10), 
    m = 1000, 
    steps = .025, 
    weights = list(0.20,1),
    m_multiplier = 5,
    mean_N = 75,
    verbose = verbose,
    seed = 20230728,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )
  
  check_step_score_hessian_bias(
    mean_smd = c(0.4, 0.1, -0.2), 
    tau = 0.1, 
    m = 1000, 
    steps = .5, 
    weights = list(0.20,0.50,0.10),
    m_multiplier = 5,
    mean_N = 120,
    verbose = verbose,
    seed = 20230729,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )
  
  check_step_score_hessian_bias(
    mean_smd = 0.2, 
    tau = c(0.1,0.05), 
    m = 1000, 
    steps = .05, 
    weights = list(0.10,0.05),
    m_multiplier = 5,
    verbose = verbose,
    seed = 20230730,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )
  
  check_step_score_hessian_bias(
    mean_smd = c(0.25, 0.15, -0.15), 
    tau = 0.2, 
    m = 1000, 
    steps = c(.025,.100), 
    weights = list(c(0.10,0.50), c(0.05, 0.20), c(0.20, 0.30)),
    m_multiplier = 5,
    verbose = verbose,
    seed = 20240403,
    hessian_threshold = NULL,
    score_type = "hybrid"
  )
  
})


test_that("step_hybrid_score and step_hybrid_jacobian agree with numerical derivatives.", {
  
  which_hess <- cbind(
    c(1,2,2,3,3,4,4),
    c(2,1,2,3,4,3,4)
  )
  which_hess3 <- which_hess[1:4,]
  which_hess5 <- rbind(which_hess, cbind(
    c(3,4,5,5,5),
    c(5,5,3,4,5)
  ))
  
  set.seed(20240425)
  
  dat <- r_meta(
    mean_smd = 0.2, 
    tau = 0.4, 
    omega = 0, 
    m = 50, 
    cor_mu = 0, 
    cor_sd = 0.01, 
    censor_fun = step_fun(cut_vals = c(.025, .50), weights = c(0.5, 0.2)), 
    n_ES_sim = n_ES_param(40, 1) 
  )
  
  step_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "step",
    steps = c(.025, .500),
    estimator = "hybrid",
    optimizer = c("nleqslv","rootSolve"),
    crit = qnorm(0.975),
    N = 200
  )
  step_derivs$selmod_fit$est
  step_derivs$score_diff_over_range
  expect_lt(max(step_derivs$score_diff_over_range[1:2]), 1e-4)
  round(step_derivs$hess_diff_over_range, 5)
  expect_lt(max(step_derivs$hess_diff_over_range[which_hess]), 1e-3)
  
  step_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "step",
    steps = c(.05, .75),
    estimator = "hybrid",
    optimizer = c("nleqslv","rootSolve"),
    N = 200,
    crit = 2
  )
  step_derivs$selmod_fit$est
  step_derivs$score_diff_over_range
  expect_lt(max(step_derivs$score_diff_over_range[1:2]), 1e-3)
  round(step_derivs$hess_diff_over_range, 5)
  expect_lt(max(step_derivs$hess_diff_over_range[which_hess]), 5e-3)
  
  
  step_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "step",
    steps = c(.05),
    estimator = "hybrid",
    optimizer = c("nleqslv","rootSolve"),
    N = 200,
    crit = 2
  )
  step_derivs$selmod_fit$est
  step_derivs$score_diff_over_range
  expect_lt(max(step_derivs$score_diff_over_range[1:2]), 1e-4)
  round(step_derivs$hess_diff_over_range, 5)
  expect_lt(max(step_derivs$hess_diff_over_range[which_hess3]), 1e-4)
  
  
  step_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "step",
    steps = c(.025, .05, .50),
    estimator = "hybrid",
    optimizer = c("nleqslv","rootSolve"),
    N = 200,
    crit = 2
  )
  step_derivs$selmod_fit$est
  step_derivs$score_diff_over_range
  expect_lt(max(step_derivs$score_diff_over_range[1:2]), 1e-4)
  round(step_derivs$hess_diff_over_range, 5)
  expect_lt(max(step_derivs$hess_diff_over_range[which_hess5]), 1e-4)
  
  set.seed(20240425)
  dat <- r_meta(
    mean_smd = 0.1, 
    tau = 0.3, 
    omega = 0, 
    m = 50, 
    cor_mu = 0, 
    cor_sd = 0.01, 
    censor_fun = beta_fun(delta_1 = 0.6, delta_2 = 1.7, trunc_1 = .025, trunc_2 = .975), 
    n_ES_sim = n_ES_param(40, 1) 
  )
  
  
  step_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "step",
    steps = c(.02),
    estimator = "hybrid",
    optimizer = c("nleqslv","rootSolve"),
  )
  step_derivs$selmod_fit$est
  step_derivs$score_diff_over_range
  expect_lt(max(step_derivs$score_diff_over_range[1:2]), 1e-4)
  round(step_derivs$hess_diff_over_range, 5)
  expect_lt(max(step_derivs$hess_diff_over_range[which_hess3]), 1e-3)
  
  
  step_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "step",
    steps = c(.025, .500),
    estimator = "hybrid",
    optimizer = c("nleqslv","rootSolve"),
  )
  step_derivs$selmod_fit$est
  step_derivs$score_diff_over_range
  expect_lt(max(step_derivs$score_diff_over_range[1:2]), 1e-4)
  round(step_derivs$hess_diff_over_range, 5)
  expect_lt(max(step_derivs$hess_diff_over_range[which_hess]), 1e-3)
  
  step_derivs <- check_all_derivatives(
    data = dat, 
    yi = d, sei = sd_d, 
    selection_type = "step",
    steps = c(.025, .975),
    estimator = "hybrid",
    optimizer = c("nleqslv","rootSolve"),
    crit = c(2,2,2,1e7)
  )
  step_derivs$selmod_fit$est
  step_derivs$score_diff_over_range
  expect_lt(max(step_derivs$score_diff_over_range[1:2]), 3e-4)
  round(step_derivs$hess_diff_over_range, 5)
  expect_lt(max(step_derivs$hess_diff_over_range[which_hess]), 1e-2)
  
  # library(tidyverse)
  # 
  # ggplot(step_derivs$data, aes(x = val)) +
  #   facet_wrap(~ param, scales = "free") +
  #   geom_hline(yintercept = 0) +
  #   geom_line(aes(y = d_num), color = "blue", linewidth = 2) +
  #   geom_line(aes(y = score), color = "green", linewidth = 1) +
  #   theme_minimal()
  # 
  # step_derivs_hessian_comparison <-
  #   step_derivs$data %>%
  #   select(param, val, starts_with("h_num"), starts_with("hess")) %>%
  #   pivot_longer(c(starts_with("h_num"), starts_with("hess")),
  #                names_to = "hp", values_to = "hess") %>%
  #   separate_wider_delim(hp, names = c("version", "hdim"), delim = ".", too_many = "merge") %>%
  #   group_by(param, val, version) %>%
  #   mutate(
  #     hdim = colnames(step_derivs$hess_diff_over_range),
  #     param = paste(param, hdim, sep = "-")
  #   ) %>%
  #   pivot_wider(names_from = version, values_from = hess)
  # 
  # ggplot(step_derivs_hessian_comparison, aes(x = val)) +
  #   facet_wrap(~  param, scales = "free") +
  #   geom_line(aes(y = h_num), color = "blue", linewidth = 2) +
  #   geom_line(aes(y = hess), color = "green", linewidth = 1) +
  #   theme_minimal()
  
})

