
test_that("predict.step.selmodel() works for the original dataset and newdata.", {
  
  set.seed(20241031)
  
  dat <- r_meta_categories(
    mean_smd = c(0, 0.2, 0.4),
    tau = c(0.2, 0.05, 0.5),
    omega = 0,
    m = 100,
    cor_mu = 0.6,
    cor_sd = 0.01,
    censor_fun = step_fun(cut_vals = .025, weights = 0.3),
    n_ES_sim = n_ES_param(mean_N = 40, mean_ES = 1)
  )
  
  PSM3_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = 0.025
  )
  
  PSM4_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500)
  )
  
  X_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    mean_mods = ~ X
  )
  
  check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    mean_mods = ~ 0 + X,
    check_subset = FALSE
  ) |> 
    expect_equal(X_preds, tolerance = 1e-6)
  
  U_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    var_mods = ~ X
  )
  
  check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    var_mods = ~ 0 + X,
    check_subset = FALSE
  ) |> 
    expect_equal(U_preds, tolerance = 1e-6)
  
  Z_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    sel_mods = ~ X
  )

  check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    sel_mods = ~ 0 + X,
    check_subset = FALSE
  ) |> 
    expect_equal(Z_preds, tolerance = 1e-6)
  
  
  XUZ_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    mean_mods = ~ 0 + X,
    var_mods = ~ X,
    sel_mods = ~ X
  )

  A_preds <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    subset = X == "A",
    check_subset = FALSE
  ) |> predict()
  
  B_preds <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    subset = X == "B",
    check_subset = FALSE
  ) |> predict()
  
  C_preds <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = c(.025, .500),
    subset = X == "C"
  ) |> predict()

  ABC_preds <- rbind(A_preds, B_preds, C_preds)
  
  expect_equal(ABC_preds, XUZ_preds)

})


test_that("predict.beta.selmodel() works for the original dataset and newdata.", {
  
  skip_on_cran()
  
  set.seed(20241031)
  
  dat <- r_meta_categories(
    mean_smd = c(0.1, 0.2, 0.6),
    tau = c(0.2, 0.05, 0.5),
    omega = 0,
    m = 20,
    cor_mu = 0.6,
    cor_sd = 0.01,
    censor_fun = beta_fun(delta_1 = 0.6, delta_2 = 0.9, trunc_1 = 0.025, trunc_2 = 0.500),
    n_ES_sim = n_ES_param(mean_N = 40, mean_ES = 1)
  )
  
  beta_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "beta",
    estimator = "ML",
    steps = c(.025, .500)
  )
  
  X_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "beta",
    estimator = "ML",
    steps = c(.025, .500),
    mean_mods = ~ X
  )
  
  check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "beta",
    estimator = "ML",
    steps = c(.025, .500),
    mean_mods = ~ 0 + X,
    check_subset = FALSE
  ) |> 
    expect_equal(X_preds, tolerance = 1e-5)
  
  U_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "beta",
    estimator = "ML",
    steps = c(.025, .500),
    var_mods = ~ X
  )
  
  check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "beta",
    estimator = "ML",
    steps = c(.025, .500),
    var_mods = ~ 0 + X,
    check_subset = FALSE
  ) |> 
    expect_equal(U_preds, tolerance = 5e-3)
  
  XU_preds <- check_predictions(
    data = dat,
    yi = d,
    sei = sd_d,
    selection_type = "beta",
    estimator = "ML",
    steps = c(.025, .500),
    mean_mods = ~ 0 + X,
    var_mods = ~ X
  )
  
})
