test_that("selection_wts() handles default values properly.", {
  
  set.seed(20241030)
  
  dat <- r_meta(
    mean_smd = 0, 
    tau = .1, omega = .01,
    m = 100, 
    cor_mu = .4, cor_sd = 0.001, 
    censor_fun = step_fun(cut_vals = c(.025, .05, .50, .9), weights = c(0.6, 0.4, 0.2, 0.1)), 
    n_ES_sim = n_ES_param(40, 3)
  )
  
  check_selection_weights(dat, steps = .025)
  check_selection_weights(dat, steps = c(.05, .50))
  check_selection_weights(dat, steps = c(.025, .05, .50))
  check_selection_weights(dat, steps = c(.025, .05, .50, .9))
  
  check_selection_weights(dat, steps = .025, bootstrap = "multinomial")
  check_selection_weights(dat, steps = c(.05, .50), bootstrap = "exp")
  
})
