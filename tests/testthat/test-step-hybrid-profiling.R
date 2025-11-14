# Make up some fake data

set.seed(20240502)

dat <- r_meta_categories(
  mean_smd = c(0, 0.2, 0.4),
  tau = c(0.2, 0.05, 0.15),
  omega = 0,
  m = 50,
  cor_mu = 0.6,
  cor_sd = 0.01,
  censor_fun = step_fun(cut_vals = .025, weights = 0.4),
  n_ES_sim = n_ES_param(mean_N = 40, mean_ES = 1)
)

pred_mat <- model.matrix(~ X, data = dat)

test_that("parse_step_params works when profiling beta.", {
  
  verbose <- FALSE
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = .025, verbose = verbose, score_tol = 2e-8)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.025, .05), verbose = verbose, score_tol = Inf, jac_tol = 1e-3)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.1, .5, .7), tol = 1e-3, score_tol = Inf, jac_tol = 1e-3, verbose = verbose)

  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = .05, X = pred_mat, score_tol = 1e-6, verbose = verbose)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.025, .500), X = pred_mat, tol = c(5e-3, 1e-4, 1e-4), score_tol = Inf, jac_tol = 1e-3, verbose = verbose)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.1, .5, .7), X = pred_mat, tol = c(5e-3, 1e-4, 1e-4), score_tol = 1e-7, verbose = verbose)
  
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = .05, U = pred_mat, score_tol = 5e-7, verbose = verbose)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.03, .40), U = pred_mat, score_tol = Inf, jac_tol = 1e-3, verbose = verbose)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.1, .5, .7), U = pred_mat, tol = c(1e-4, 1e-8, 1e-8), score_tol = Inf, jac_tol = 1e-3, verbose = verbose)

  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = .025, Z = list(pred_mat), tol = c(1e-4, 1e-8, 1e-8), score_tol = Inf, jac_tol = 1e-3, verbose = verbose)
  
})

test_that("parse_step_params works when profiling beta with non-flat priors.", {
  
  prior_spec <- define_priors(beta_mean = 0.3, beta_precision = 0.3^2 / 2)
  verbose <- FALSE
  
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = .02, priors = prior_spec, verbose = verbose, score_tol = 1e-7, jac_tol = 5e-3)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.025, .05), priors = prior_spec, verbose = verbose, score_tol = Inf, jac_tol = 5e-3)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.1, .5, .7), priors = prior_spec, tol = 1e-3, score_tol = Inf, jac_tol = 5e-3, verbose = verbose)
  
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = .05, X = pred_mat, priors = prior_spec, score_tol = 5e-7, jac_tol = 5e-3, verbose = verbose)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.025, .500), X = pred_mat, priors = prior_spec, tol = c(5e-3, 1e-4, 1e-4), score_tol = Inf, jac_tol = 5e-3, verbose = verbose)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.1, .6), X = pred_mat, priors = prior_spec, tol = c(5e-3, 1e-4, 1e-4), score_tol = 1e-6, jac_tol = 5e-3, verbose = verbose)
  
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = .05, U = pred_mat, priors = prior_spec, verbose = verbose, score_tol = 5e-3, jac_tol = 1e-2)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.025, .500), U = pred_mat, priors = prior_spec, score_tol = Inf, jac_tol = 1e-2, verbose = verbose)
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = c(.1, .7), U = pred_mat, priors = prior_spec, tol = c(2e-4, 1e-8, 1e-8), score_tol = Inf, jac_tol = 1e-2, verbose = verbose)
  
  check_profiling_equivalence(yi = dat$d, sei = dat$sd_d, steps = .025, Z = list(pred_mat), priors = prior_spec, tol = c(1e-4, 1e-8, 1e-8), score_tol = Inf, jac_tol = 1e-3, verbose = verbose)
  
})
