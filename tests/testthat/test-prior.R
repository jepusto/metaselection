
test_that("default_priors() works for default priors.", {

  res <- default_priors()
  expect_s3_class(res, "selmodel_prior")
  expect_identical(names(res), c("log_prior","score_prior","hessian_prior"))
  
  lp <- res$log_prior(beta = rnorm(1,0,5), gamma = rnorm(1,0,2), zeta = rnorm(1,0,2))
  expect_identical(length(lp), 1L)
  
  sp <- res$score_prior(beta = rnorm(1,0,5), gamma = rnorm(1,0,2), zeta = rnorm(1,0,2))
  expect_identical(length(sp), 3L)
  
  sh <- res$hessian_prior(beta = rnorm(3,0,5), gamma = rnorm(1,0,2), zeta = rnorm(2,0,2))
  expect_identical(dim(sh), c(6L, 6L))
  expect_identical(sh[lower.tri(sh)], rep(0,15))
  expect_lt(max(diag(sh)),0)
  
})


test_that("priors = NULL is equivalent to defining flat priors.", {
  
  dat <- r_meta_categories(
    mean_smd = c(0, 0.1, 0.2),
    tau = 0.3,
    omega = 0,
    m = 12,
    cor_mu = 0.6,
    cor_sd = 0.01,
    censor_fun = step_fun(cut_vals = .025, weights = 0.4),
    n_ES_sim = n_ES_param(mean_N = 40, mean_ES = 1)
  )
  
  prior_spec <- default_priors(beta_precision = 1e-6, tau_alpha = 1e-6, lambda_precision = 1e-6)
  
  cml_flat <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "step",
    steps = .025,
    priors = NULL,
    estimator = "CML"
  )
  
  cml_prior <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "step",
    steps = .025,
    priors = prior_spec,
    estimator = "CML"
  )
  
  expect_equal(cml_flat$est$Est, cml_prior$est$Est, tolerance = 1e-5)
  expect_equal(cml_flat$est$SE, cml_prior$est$SE, tolerance = 1e-5)
  
  argl_flat <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "step",
    steps = .025,
    priors = NULL,
    estimator = "ARGL"
  )
  
  argl_prior <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "step",
    steps = .025,
    priors = prior_spec,
    estimator = "ARGL"
  )
  
  expect_equal(argl_flat$est$Est, argl_prior$est$Est, tolerance = 1e-5)
  expect_equal(argl_flat$est$SE, argl_prior$est$SE, tolerance = 1e-5)
  
  
  fargl_flat <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "step",
    steps = .025,
    priors = NULL,
    estimator = "ARGL-full"
  )
  
  fargl_prior <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "step",
    steps = .025,
    priors = prior_spec,
    estimator = "ARGL-full"
  )
  
  expect_equal(fargl_flat$est$Est, fargl_prior$est$Est, tolerance = 1e-5)
  expect_equal(fargl_flat$est$SE, fargl_prior$est$SE, tolerance = 1e-5)
  
  

})

test_that("Score contributions sum to total when accounting for priors.", {
  
  dat <- r_meta_categories(
    mean_smd = c(0, 0.1, 0.2),
    tau = 0.2,
    omega = 0,
    m = 12,
    cor_mu = 0.6,
    cor_sd = 0.01,
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.4, 0.2)),
    n_ES_sim = n_ES_param(mean_N = 40, mean_ES = 1)
  )
  
  prior_spec <- default_priors(beta_mean = 0.1, beta_precision = 25)
  
  cml_prior <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "step",
    steps = c(.025,.500),
    priors = prior_spec,
    estimator = "CML"
  )
  
  cml_score <- step_score(
    theta = cml_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec
  )
  
  cml_score_cont <- step_score(
    theta = cml_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec,
    contributions = TRUE
  )
  
  expect_identical(dim(cml_score_cont), c(nrow(dat), length(cml_score)))
  expect_lt(max(abs(colSums(cml_score_cont) - cml_score)), 1e-14)
  
  cml_hyscore <- step_hybrid_score(
    theta = cml_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec
  )
  
  cml_hyscore_cont <- step_hybrid_score(
    theta = cml_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec,
    contributions = TRUE
  )
  
  expect_identical(dim(cml_hyscore_cont), c(nrow(dat), length(cml_hyscore)))
  expect_lt(max(abs(colSums(cml_hyscore_cont) - cml_hyscore)), 1e-14)
  
  
  argl_prior <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "step",
    steps = c(.025,.500),
    priors = prior_spec,
    estimator = "ARGL"
  )
  
  argl_score <- step_score(
    theta = argl_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec
  )
  
  argl_score_cont <- step_score(
    theta = argl_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec,
    contributions = TRUE
  )
  
  expect_identical(dim(argl_score_cont), c(nrow(dat), length(argl_score)))
  expect_lt(max(abs(colSums(argl_score_cont) - argl_score)), 1e-14)
  
  argl_hyscore <- step_hybrid_score(
    theta = argl_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec
  )
  
  argl_hyscore_cont <- step_hybrid_score(
    theta = argl_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec,
    contributions = TRUE
  )
  
  expect_identical(dim(argl_hyscore_cont), c(nrow(dat), length(argl_hyscore)))
  expect_lt(max(abs(colSums(argl_hyscore_cont) - argl_hyscore)), 1e-14)
  
  
  fargl_prior <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "step",
    steps = c(.025,.500),
    priors = prior_spec,
    estimator = "ARGL-full"
  )
  
  fargl_score <- step_score(
    theta = fargl_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec
  )
  
  fargl_score_cont <- step_score(
    theta = fargl_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec,
    contributions = TRUE
  )
  
  expect_identical(dim(fargl_score_cont), c(nrow(dat), length(fargl_score)))
  expect_lt(max(abs(colSums(fargl_score_cont) - fargl_score)), 1e-14)
  
  fargl_hyscore <- step_hybrid_score(
    theta = fargl_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec
  )
  
  fargl_hyscore_cont <- step_hybrid_score(
    theta = fargl_prior$est$Est, 
    yi = dat$d, sei = dat$sd_d, steps = c(.025,.500), 
    priors = prior_spec,
    contributions = TRUE
  )
  
  expect_identical(dim(fargl_hyscore_cont), c(nrow(dat), length(fargl_hyscore)))
  expect_lt(max(abs(colSums(fargl_hyscore_cont) - fargl_hyscore)), 1e-14)
  
})
