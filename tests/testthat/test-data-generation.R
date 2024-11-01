test_that("beta_wt_fun() is always between 0 and 1.", {
  
  plot <- FALSE
  beta_range <- eval_beta_range(c(1, 1), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)

  beta_range <- eval_beta_range(c(1.1, 1), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)

  beta_range <- eval_beta_range(c(1, 1.5), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)
  
  beta_range <- eval_beta_range(c(0.7, 1), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)
  
  beta_range <- eval_beta_range(c(1, 0.3), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)

  beta_range <- eval_beta_range(c(0.8, 1.2), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)

  beta_range <- eval_beta_range(c(1.4, 0.6), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)

  beta_range <- eval_beta_range(c(0.7, 1.2), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)
  
  beta_range <- eval_beta_range(c(1.4, 0.5), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)
  
  beta_range <- eval_beta_range(c(0.8, 0.8), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)
  
  beta_range <- eval_beta_range(c(0.6, 0.9), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)

  beta_range <- eval_beta_range(c(4.8, 4.8), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)
  
  beta_range <- eval_beta_range(c(2, 7), c(.025, .975), plot = plot)$range
  expect_gte(beta_range[1], 0)
  expect_lte(beta_range[2], 1)
  
})

test_that("r_meta generates target number of studies with single parameter values", {
  
  check_target_m(
    mean_smd = 0.3, 
    tau = 0.1, 
    omega = 0.1, 
    m = 10L, 
    cor_mu = 0.6, 
    censor_fun = step_fun(), 
    n_ES_sim = n_ES_param(40, 3), 
    m_multiplier = 10
  )
  
  check_target_m(
    mean_smd = 0.1, 
    tau = 0.1, 
    omega = 0.0, 
    m = 15L, 
    cor_mu = 0.4, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.001, 0.100)), 
    n_ES_sim = n_ES_param(40, 9), 
    m_multiplier = 12
  )
  
  check_target_m(
    mean_smd = -0.2, 
    tau = 0.5, 
    omega = 0.2, 
    m = 6L, 
    cor_mu = 0.8, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.01, 0.20)), 
    n_ES_sim = n_ES_param(190, 8), 
    m_multiplier = 1
  )
  
  check_target_m(
    mean_smd = 0.1, 
    tau = 0.2, 
    omega = 0.05, 
    m = 16L, 
    cor_mu = 0.8, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.01, 0.20)), 
    n_ES_sim = n_ES_param(20, 1), 
    m_multiplier = 2
  )

  check_target_m(
    mean_smd = 0.0, 
    tau = 0.0, 
    omega = 0.01, 
    m = 70L, 
    cor_mu = 0.3, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.01, 0.20)), 
    n_ES_sim = n_ES_param(10, 1, min_N = 6), 
    m_multiplier = 1.5
  )
  
})


test_that("r_meta generates target number of studies with multiple parameter values", {
  
  check_target_m(
    mean_smd = c(0.1, 0.3, 0.5), 
    tau = 0.1, 
    omega = 0.1, 
    m = 3L, 
    cor_mu = 0.6, 
    censor_fun = step_fun(cut_vals = 0.4, weights = 0.1), 
    n_ES_sim = n_ES_param(40, 3), 
    m_multiplier = 10
  )
  
  check_target_m(
    mean_smd = 0.1, 
    tau = c(0.0, 0.1), 
    omega = 0.0, 
    m = 15L, 
    cor_mu = 0.4, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.001, 0.100)), 
    n_ES_sim = n_ES_param(40, 9), 
    m_multiplier = 12
  )
  
  check_target_m(
    mean_smd = -0.2, 
    tau = 0.5, 
    omega = 0.2, 
    m = c(6L,3L,2L,3L),
    cor_mu = 0.8, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.01, 0.20)), 
    n_ES_sim = n_ES_param(190, 8), 
    m_multiplier = 1
  )
  
  check_target_m(
    mean_smd = c(0.1,0.2,0.3), 
    tau = c(0.1,0.2), 
    omega = 0.05, 
    m = c(16L,4L), 
    cor_mu = 0.8, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.01, 0.20)), 
    n_ES_sim = n_ES_param(20, 1), 
    m_multiplier = 1
  )
  
  check_target_m(
    mean_smd = 0.0, 
    tau = 0.0, 
    omega = 0.01, 
    m = 70L, 
    cor_mu = 0.3, 
    censor_fun = list(
      step_fun(cut_vals = c(.025, .500), weights = c(0.01, 0.20)),
      step_fun(cut_vals = c(.025), weights = c(0.20)),
      step_fun(cut_vals = c(.025), weights = c(1))
    ), 
    n_ES_sim = n_ES_param(10, 1, min_N = 6), 
    m_multiplier = 1.5
  )
    
})



test_that("r_meta generates ES estimates with expected structure", {
  
  skip_on_cran()
  skip_if_not_installed("metafor")
  verbose <- FALSE
  
  set.seed(20230621)
  check_model_structure(
    mean_smd = 0.0, tau = 0.21, omega = 0.06, 
    m = 180, 
    n_ES_sim = n_ES_param(100, 3),
    verbose = verbose
  )
  
  set.seed(20240505)
  check_model_structure(
    mean_smd = 0.3, tau = 0.08, omega = 0.01, 
    cor_mu = 0.3,
    m = 250, 
    n_ES_sim = n_ES_param(140, 2),
    verbose = verbose
  )

  set.seed(20230623)
  check_model_structure(
    mean_smd = 0.5, tau = 0.8, omega = 0.4, 
    cor_mu = 0.4,
    m = 200, 
    n_ES_sim = n_ES_param(1000, 2),
    verbose = verbose
  )
  
  set.seed(20230624)
  check_model_structure(
    mean_smd = 0.1, tau = 0.1, omega = 0.1, 
    cor_mu = 0,
    m = 250, 
    n_ES_sim = n_ES_param(420, 2),
    verbose = verbose
  )

})


test_that("r_meta_categories generates ES estimates with expected structure", {
  
  skip_on_cran()
  skip_if_not_installed("metafor")
  verbose <- FALSE

  set.seed(20230621)
  check_model_structure(
    mean_smd = seq(0, 0.5, 0.15),
    tau = 0.21, omega = 0, 
    cor_mu = 0,
    m = 80, 
    n_ES_sim = n_ES_param(100, 5),
    verbose = verbose,
    check_varcomp = FALSE
  )
  
  set.seed(20230622)
  check_model_structure(
    mean_smd = c(0.3, 0.07),
    tau = c(0.08, 0.28),
    omega = 0.01, 
    cor_mu = 0,
    m = 100, 
    n_ES_sim = n_ES_param(70, 2),
    verbose = verbose
  )
  
  set.seed(20230623)
  check_model_structure(
    mean_smd = c(0.3, 0.07),
    tau = 0.8, 
    omega = c(0.06, 0.23), 
    cor_mu = 0,
    m = 75, 
    n_ES_sim = n_ES_param(140, 2),
    verbose = verbose
  )
  
  set.seed(20230624)
  check_model_structure(
    mean_smd = c(0.10, 0.13),
    tau = c(0.20, 0.13),
    omega = c(0.1, 0.08), 
    cor_mu = 0,
    m = 90, 
    n_ES_sim = n_ES_param(160, 2),
    verbose = verbose
  )
  
})

