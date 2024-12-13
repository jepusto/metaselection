
test_that("build_model_frame works properly.", {

  k <- 10L
  dat <- data.frame(
    smd = rnorm(k),
    se = 1 / (3 + rpois(k, 20)),
    pvalues = runif(k),
    aweights = rpois(k, 5),
    x1 = rnorm(k),
    v1 = rbeta(k, 3, 2),
    z1 = sample(LETTERS[1:3], size = k, replace = TRUE),
    z2 = runif(k),
    sz1 = sample(LETTERS[5:6], size = k, replace = TRUE),
    sz2 = rnorm(k),
    cl = sample(LETTERS[1:4], size = k, replace = TRUE)
  )
  dat$v <- dat$se^2
  
  build_model_frame(data = dat, yi = smd, sei = se) |>
    check_dims(k, 2L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues) |>
    check_dims(k, 3L)
  build_model_frame(data = dat, yi = smd, sei = se, ai = aweights) |>
    check_dims(k, 3L)
  build_model_frame(data = dat, yi = smd, vi = v) |>
    check_dims(k, 2L)
  build_model_frame(data = dat, yi = smd, vi = v, pi = pvalues) |>
    check_dims(k, 3L)
  build_model_frame(data = dat, yi = smd, vi = v, ai = aweights) |>
    check_dims(k, 3L)
  build_model_frame(data = dat, yi = smd, vi = v, pi = pvalues, ai = aweights) |>
    check_dims(k, 4L)
  
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, subset = NULL) |>
    check_dims(k, 4L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, subset = cl != "D") |>
    check_dims(sum(dat$cl != "D"), 4L)  
  build_model_frame(data = dat, yi = smd, vi = v, pi = pvalues, ai = aweights, subset = NULL) |>
    check_dims(k, 4L)
  build_model_frame(data = dat, yi = smd, vi = v, pi = pvalues, ai = aweights, subset = cl != "D") |>
    check_dims(sum(dat$cl != "D"), 4L)  
  
  build_model_frame(data = dat, yi = smd, sei = se, cluster = cl) |>
    check_dims(k, 3L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, cluster = cl) |>
    check_dims(k, 4L)
  build_model_frame(data = dat, yi = smd, sei = se, ai = aweights, cluster = cl) |>
    check_dims(k, 4L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, cluster = cl) |>
    check_dims(k, 5L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, cluster = cl, subset = sz1 == "E") |>
    check_dims(sum(dat$sz1 == "E"), 5L)
  
  build_model_frame(data = dat, yi = smd, sei = se, mean_mods = ~ 0 + x1) |>
    check_dims(k, 3L)
  build_model_frame(data = dat, yi = smd, sei = se, ai = aweights, var_mods = ~ v1) |>
    check_dims(k, 4L)
  build_model_frame(data = dat, yi = smd, sei = se, sel_mods = ~ z1 + z2 + z1:z2) |>
    check_dims(k, 4L)
  build_model_frame(data = dat, yi = smd, sei = se, sel_mods = list(~ z1, ~ z1 + z2)) |>
    check_dims(k, 4L)
  build_model_frame(data = dat, yi = smd, sei = se, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 4L)
  build_model_frame(data = dat, yi = smd, sei = se, sel_zero_mods = ~ 0 + sz1 * sz2, subset = aweights > 4) |>
    check_dims(sum(dat$aweights > 4), 4L)
  build_model_frame(data = dat, yi = smd, vi = v, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 4L)
  build_model_frame(data = dat, yi = smd, vi = v, sel_zero_mods = ~ 0 + sz1 * sz2, subset = aweights > 4) |>
    check_dims(sum(dat$aweights > 4), 4L)
  
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, mean_mods = ~ 0 + x1, var_mods = ~ v1) |>
    check_dims(k, 5L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, mean_mods = ~ 0 + x1, sel_mods = ~ z1 + z2 + z1:z2) |>
    check_dims(k, 6L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, mean_mods = ~ 0 + x1, sel_mods = list(~ z1, ~ z1 + z2)) |>
    check_dims(k, 6L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, mean_mods = ~ 0 + x1, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 6L)
  
  build_model_frame(data = dat, yi = smd, sei = se, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2) |>
    check_dims(k, 5L)
  build_model_frame(data = dat, yi = smd, sei = se, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2)) |>
    check_dims(k, 5L)
  build_model_frame(data = dat, yi = smd, sei = se, var_mods = ~ v1, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 5L)

  build_model_frame(data = dat, yi = smd, sei = se, ai = aweights, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 7L)
  build_model_frame(data = dat, yi = smd, sei = se, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 6L)

  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2) |>
    check_dims(k, 8L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2)) |>
    check_dims(k, 8L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 8L)

  build_model_frame(data = dat, yi = smd, sei = se, mean_mods = ~ 0 + x1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 7L)
  build_model_frame(data = dat, yi = smd, sei = se, mean_mods = ~ 0 + x1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 7L)

  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 8L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 8L)

  build_model_frame(data = dat, yi = smd, sei = se, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 8L)
  build_model_frame(data = dat, yi = smd, sei = se, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 8L)
  
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 9L)
  build_model_frame(data = dat, yi = smd, sei = se, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 9L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 10L)
  
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, cluster = cl, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 10L)
  build_model_frame(data = dat, yi = smd, sei = se, ai = aweights, cluster = cl, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 10L)
  build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, cluster = cl, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2) |>
    check_dims(k, 11L)
  
  # Error if yi or (sei and vi) or data arguments missing
  expect_error(
    build_model_frame(data = dat, sei = se, pi = pvalues, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2)     
  )
  expect_error(
    build_model_frame(data = dat, vi = v, pi = pvalues, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2)     
  )
  expect_error(
    build_model_frame(data = dat, yi = smd, pi = pvalues, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2) 
  )
  expect_error(
    build_model_frame(yi = smd, sei = se, pi = pvalues, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2)     
  )
  expect_error(
    build_model_frame(yi = smd, vi = v, sei = se, pi = pvalues, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = ~ z1 + z2 + z1:z2, sel_zero_mods = ~ 0 + sz1 * sz2)     
  )
  
  # Error if variable names do not occur in data
  expect_error(
    build_model_frame(data = dat, yi = di, sei = se, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2)
  )
  expect_error(
    build_model_frame(data = dat, yi = smd, sei = sei, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2)
  )
  expect_error(
    build_model_frame(data = dat, yi = smd, vi = Var_i, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2)
  )
  expect_error(
    build_model_frame(data = dat, yi = smd, sei = se, pi = pval, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2)
  )
  expect_error(
    build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = awt, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2)
  )
  expect_error(
    build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2)
  )
  expect_error(
    build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1a, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2)
  )
  expect_error(
    build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1a, ~ z1 + z2), sel_zero_mods = ~ 0 + sz1 * sz2)
  )
  expect_error(
    build_model_frame(data = dat, yi = smd, sei = se, pi = pvalues, ai = aweights, mean_mods = ~ 0 + x1, var_mods = ~ v1, sel_mods = list(~ z1, ~ z1 + z2), sel_zero_mods = ~ 0 + szb1 * sz2)
  )
})


test_that("selection_model() returns results of correct dimension when estimator = 'ML'.", {
  
  # Generate some independent data
  
  set.seed(20230523)
  dat <- r_meta(
    mean_smd = 0.3, tau = 0.1, omega = 0,
    m = 1000, cor_mu = 0.6, cor_sd = 0.001, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.5, 0.2)), 
    n_ES_sim = n_ES_param(40, 1)
  ) 
  dat$X1 <- rnorm(nrow(dat))
  dat$Z1 <- sample(LETTERS[1:2], size = nrow(dat), replace = TRUE)

  methods <- c("BFGS","CG","Nelder-Mead","nlminb","Rvmmin","bobyqa")
  
  m1_est <- selection_model(data = dat, yi = d, sei = sda, steps = c(.025, .500), 
                            vcov_type = "raw", 
                            optimizer_control = list(all.methods = TRUE))
  expect_s3_class(m1_est, "data.frame")
  expect_true(all(methods %in% rownames(m1_est)))
  
  m1a <- selection_model(data = dat, yi = d, sei = sda, steps = c(.025, .500), 
                         optimizer = methods)
  
  expect_equal(nrow(m1a$est), 4L)
  expect_equal(nrow(m1a$est), nrow(m1a$vcov))
  
  m1b <- selection_model(data = dat, yi = d, vi = Va, steps = c(.025, .500), 
                         cluster = studyid,
                         optimizer = methods)
  expect_equal(m1a$est, m1b$est)
  expect_equal(m1a$vcov, m1b$vcov)
  
  m2_est <- selection_model(data = dat, yi = d, vi = Va, 
                            mean_mods = ~ X1, sel_mods = ~ Z1, 
                            steps = c(.025, .500), 
                            vcov_type = "raw", 
                            optimizer_control = list(all.methods = TRUE))
  expect_s3_class(m2_est, "data.frame")
  expect_true(all(methods %in% rownames(m2_est)))
  
  m2a <- selection_model(data = dat, yi = d, vi = Va, 
                         mean_mods = ~ X1, sel_mods = ~ Z1, 
                         steps = c(.025, .500), 
                         optimizer = methods)
  expect_equal(nrow(m2a$est), 7L)
  expect_equal(nrow(m2a$est), nrow(m2a$vcov))

  m2b <- selection_model(data = dat, yi = d, sei = sda, 
                         cluster = studyid,
                         mean_mods = ~ X1, sel_mods = ~ Z1, 
                         steps = c(.025, .500), 
                         optimizer = methods)

  expect_equal(m2a$est, m2b$est)
  expect_equal(m2a$vcov, m2b$vcov)
  
})

test_that("selection_model() returns results of correct dimension when estimator = 'hybrid'.", {
  
  # Generate some independent data
  
  set.seed(20230920)
  dat <- r_meta(
    mean_smd = 0.3, tau = 0.1, omega = 0,
    m = 1000, cor_mu = 0.6, cor_sd = 0.001, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.5, 0.2)), 
    n_ES_sim = n_ES_param(40, 1)
  ) 
  dat$X1 <- rnorm(nrow(dat))
  dat$Z1 <- sample(LETTERS[1:2], size = nrow(dat), replace = TRUE)
  
  
  m1_est <- selection_model(data = dat, yi = d, sei = sda, steps = c(.025, .500), 
                            vcov_type = "raw", 
                            estimator = "hybrid")
  expect_identical(names(m1_est), c("est","max_method","info"))
  expect_identical(names(m1_est$info$nleqslv), c("fvec","termcd","message","scalex","nfcnt","njcnt","iter","jac","f_norm"))
  
  m1a <- selection_model(data = dat, yi = d, sei = sda, steps = c(.025, .500), 
                         estimator = "hybrid")
  
  expect_equal(nrow(m1a$est), 4L)
  expect_equal(nrow(m1a$est), nrow(m1a$vcov))
  
  m1b <- selection_model(data = dat, yi = d, vi = Va, steps = c(.025, .500), 
                         cluster = studyid, estimator = "hybrid")
  expect_equal(m1a$est, m1b$est)
  expect_equal(m1a$vcov, m1b$vcov)
  
  m2_est <- selection_model(data = dat, yi = d, vi = Va, 
                            mean_mods = ~ X1, sel_mods = ~ Z1, 
                            steps = c(.025, .500), 
                            vcov_type = "raw", 
                            estimator = "hybrid")
  expect_identical(names(m2_est), c("est","max_method","info"))
  expect_identical(names(m2_est$info$nleqslv), c("fvec","termcd","message","scalex","nfcnt","njcnt","iter","jac","f_norm"))
  
  m2a <- selection_model(data = dat, yi = d, vi = Va, 
                         mean_mods = ~ X1, sel_mods = ~ Z1, 
                         steps = c(.025, .500),
                         estimator = "hybrid")
  expect_equal(nrow(m2a$est), 7L)
  expect_equal(nrow(m2a$est), nrow(m2a$vcov))
  
  m2b <- selection_model(data = dat, yi = d, sei = sda, 
                         cluster = studyid,
                         mean_mods = ~ X1, sel_mods = ~ Z1, 
                         steps = c(.025, .500), 
                         estimator = "hybrid")
  
  expect_equal(m2a$est, m2b$est)
  expect_equal(m2a$vcov, m2b$vcov)
  
})


test_that("selection_model() works with the subset argument.", {
  
  set.seed(20240410)
  dat <- r_meta(
    mean_smd = 0.3, tau = 0.15, omega = 0,
    m = 200, cor_mu = 0.6, cor_sd = 0.001, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.7, 0.4)), 
    n_ES_sim = n_ES_param(40, 3)
  )
  dat$Z1 <- sample(LETTERS[1:2], size = nrow(dat), replace = TRUE)
  
  
  m1_mle_A1 <-
    selection_model(
      data = subset(dat, Z1 == "A"),
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      vcov_type = "robust",
      estimator = "ML",
      optimizer = "Rvmmin"
    )
  m1_mle_A2 <-
    selection_model(
      data = dat,
      subset = Z1 == "A",
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      vcov_type = "robust",
      estimator = "ML",
      optimizer = "Rvmmin"
    )
  expect_identical(m1_mle_A1$est, m1_mle_A2$est)
  expect_identical(m1_mle_A1$vcov, m1_mle_A2$vcov)
  
  m1_mle_B1 <-
    selection_model(
      data = subset(dat, Z1 == "B"),
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      vcov_type = "robust",
      estimator = "ML",
      optimizer = "Rvmmin"
    )
  m1_mle_B2 <-
    selection_model(
      data = dat,
      subset = Z1 == "B",
      yi = d,
      vi = Va,
      steps = c(.025, .500),
      vcov_type = "robust",
      estimator = "ML",
      optimizer = "Rvmmin"
    )
  expect_identical(m1_mle_B1$est, m1_mle_B2$est)
  expect_identical(m1_mle_B1$vcov, m1_mle_B2$vcov)
  
  m1_mle_mod <-
    selection_model(
      data = dat,
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      mean_mods = ~ 0 + Z1,
      var_mods = ~ 0 + Z1,
      sel_mods = ~ 0 + Z1,
      vcov_type = "robust",
      estimator = "ML",
      optimizer = "Rvmmin"
    )

  expect_equal(
    subset(m1_mle_mod$est, substr(param, nchar(param)-3, nchar(param)) == "_Z1A", c(Est, SE)),
    subset(m1_mle_A1$est, select = c(Est, SE)),
    ignore_attr = TRUE,
    tolerance = 1e-6
  )
  
  expect_equal(
    subset(m1_mle_mod$est, substr(param, nchar(param)-3, nchar(param)) == "_Z1B", c(Est, SE)),
    subset(m1_mle_B1$est, select = c(Est, SE)),
    ignore_attr = TRUE,
    tolerance = 1e-6    
  )
  
  
  m1_hybrid_A1 <-
    selection_model(
      data = subset(dat, Z1 == "A"),
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      calc_vcov = TRUE,
      estimator = "hybrid"
    )
  m1_hybrid_A2 <-
    selection_model(
      data = dat,
      subset = Z1 == "A",
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      calc_vcov = TRUE,
      estimator = "hybrid"
    )

  expect_identical(m1_hybrid_A1$est, m1_hybrid_A2$est)
  expect_identical(m1_hybrid_A1$vcov, m1_hybrid_A2$vcov)
  expect_identical(m1_hybrid_A1$info, m1_hybrid_A2$info)
  
  
  m1_hybrid_B1 <-
    selection_model(
      data = subset(dat, Z1 == "B"),
      yi = d,
      vi = Va,
      steps = c(.025, .500),
      vcov_type = "robust",
      estimator = "hybrid"
    )
  m1_hybrid_B2 <-
    selection_model(
      data = dat,
      subset = Z1 == "B",
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      vcov_type = "robust",
      estimator = "hybrid"
    )
  
  expect_identical(m1_hybrid_B1$est, m1_hybrid_B2$est)
  expect_identical(m1_hybrid_B1$vcov, m1_hybrid_B2$vcov)
  expect_identical(m1_hybrid_B1$info, m1_hybrid_B2$info)
  
  
  set.seed(20240912)
  dat <- r_meta(
    mean_smd = 0.3, tau = 0.15, omega = 0,
    m = 800, cor_mu = 0.6, cor_sd = 0.001, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.7, 0.4)), 
    n_ES_sim = n_ES_param(40, 1)
  )
  dat$Z1 <- sample(LETTERS[1:2], size = nrow(dat), replace = TRUE)
  
  m1_hybrid_full_A1 <-
    selection_model(
      data = subset(dat, Z1 == "A"),
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      vcov_type = "robust",
      estimator = "hybrid-full"
    )
  
  m1_hybrid_full_B1 <-
    selection_model(
      data = subset(dat, Z1 == "B"),
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      vcov_type = "robust",
      estimator = "hybrid-full"
    )
  
  m1_hybrid_mod <-
    selection_model(
      data = dat,
      yi = d,
      sei = sda,
      steps = c(.025, .500),
      mean_mods = ~ 0 + Z1,
      var_mods = ~ 0 + Z1,
      sel_mods = ~ 0 + Z1,
      vcov_type = "robust",
      estimator = "hybrid-full"
    )
  
  expect_equal(
    subset(m1_hybrid_mod$est, substr(param, nchar(param)-3, nchar(param)) == "_Z1A", select = c(-param,-estimator,-p_value)),
    subset(m1_hybrid_full_A1$est, select = c(-param,-estimator,-p_value)),
    ignore_attr = TRUE,
    tolerance = 1e-6
  )
  
  expect_equal(
    subset(m1_hybrid_mod$est, substr(param, nchar(param)-3, nchar(param)) == "_Z1B", select = c(-param,-estimator,-p_value)),
    subset(m1_hybrid_full_B1$est, select = c(-param,-estimator,-p_value)),
    ignore_attr = TRUE,
    tolerance = 1e-6
  )
  
})
