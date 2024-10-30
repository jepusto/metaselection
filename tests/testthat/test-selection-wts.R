set.seed(20241030)

dat <- r_meta(
  mean_smd = 0, 
  tau = .1, omega = .01,
  m = 100, 
  cor_mu = .4, cor_sd = 0.001, 
  censor_fun = step_fun(cut_vals = c(.025, .05, .50, .9), weights = c(0.6, 0.4, 0.2, 0.1)), 
  n_ES_sim = n_ES_param(40, 3)
)


test_that("selection_wts() handles default values properly.", {
  
  check_selection_weights(dat, steps = .025)
  check_selection_weights(dat, steps = c(.05, .50))
  check_selection_weights(dat, steps = c(.025, .05, .50))
  check_selection_weights(dat, steps = c(.025, .05, .50, .9))
  
  check_selection_weights(dat, steps = .025, bootstrap = "multinomial")
  check_selection_weights(dat, steps = c(.05, .50), bootstrap = "exp")
  
})


test_that("selection_wts() handles ref_pval.", {
  
  step_fit <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(.025, .50),
    estimator = "hybrid"
  )
  
  wts_raw <- selection_wts(step_fit)
  wts_ref <- selection_wts(step_fit, ref_pval = runif(1, min = .1, max = .9))
  expect_equal(diff(range(wts_raw$wt / wts_ref$wt)), 0)

  
  step_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(.025, .50),
    estimator = "hybrid",
    bootstrap = "exp", 
    R = 4
  )
  
  wts_raw <- selection_wts(step_boot)
  wts_ref <- selection_wts(step_boot, ref_pval = runif(1, min = .1, max = .9))
  expect_equal(diff(range(wts_raw$wts$wt / wts_ref$wts$wt)), 0)
  expect_equal(max(tapply(wts_raw$boot_wts$wt / wts_ref$boot_wts$wt, wts_raw$boot_wts$rep, \(x) diff(range(x)))), 0)

  suppressWarnings(
    beta_fit <- selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      cluster = studyid,
      selection_type = "beta",
      steps = c(.025, .5)
    )
  )
  
  wts_raw <- selection_wts(beta_fit)
  wts_ref <- selection_wts(beta_fit, ref_pval = runif(1, min = .1, max = .9))
  expect_equal(diff(range(wts_raw$wt / wts_ref$wt)), 0)
  
  
  suppressWarnings(
    beta_boot <- selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      cluster = studyid,
      selection_type = "beta",
      steps = c(.025, .50),
      bootstrap = "multinomial", 
      R = 4
    )
  )
  
  wts_raw <- selection_wts(beta_boot)
  wts_ref <- selection_wts(beta_boot, ref_pval = runif(1, min = .1, max = .9))
  expect_equal(diff(range(wts_raw$wts$wt / wts_ref$wts$wt)), 0)
  expect_equal(max(tapply(wts_raw$boot_wts$wt / wts_ref$boot_wts$wt, wts_raw$boot_wts$rep, \(x) diff(range(x)))), 0)
  
})
