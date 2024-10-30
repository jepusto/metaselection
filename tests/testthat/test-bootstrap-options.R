skip_if_not_installed("dplyr")

set.seed(20240806)
dat <- r_meta(
  mean_smd = 0,
  tau = .1,
  omega = .01,
  m = 50,
  cor_mu = .4,
  cor_sd = 0.001,
  censor_fun = step_fun(cut_vals = .025, weights = 0.3),
  n_ES_sim = n_ES_param(40, 3)
)

test_that("bootstrap_CI options for selection_model() are irrelevant when bootstrap = 'none'.", {

  expect_error(
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "none",
      CI_type = "percentile"
    )
  )  
  
  expect_error(
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "none",
      CI_type = "percentile",
      R = 999
    )
  )
  
  expect_error(
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial",
      R = 0
    )
  )

  expect_error(
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "exp",
      CI_type = "percentile",
      R = 0
    )
  )

})  

test_that("bootstrap_CI options for selection_model() work when bootstrap = 'multinomial'.", {
  
  aseed <- 20241030
  
  set.seed(aseed)
  
  step_large <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial", 
      CI_type = "large-sample",
      R = 19
    )
  
  expect_s3_class(step_large, "selmodel")
  expect_identical(names(step_large$est), c("estimator","param","Est","SE","CI_lo","CI_hi"))
  expect_false(is.null(step_large$bootstrap_reps))

  set.seed(aseed)
  
  step_perc <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial", 
      CI_type = "percentile",
      R = 19
    )
  
  expect_s3_class(step_perc, "boot.selmodel")
  expect_identical(names(step_perc$est), c("estimator","param","Est","SE","bootstraps","percentile_lower","percentile_upper"))
  expect_identical(table(step_perc$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 19L)))
  
  set.seed(aseed)
  
  step_t <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial", 
      CI_type = "student",
      R = 19
    )
  
  expect_s3_class(step_t, "boot.selmodel")
  expect_identical(names(step_t$est), c("estimator","param","Est","SE","bootstraps","student_lower","student_upper"))
  expect_identical(table(step_t$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 19L)))
  expect_identical(step_t$bootstrap_reps$Est, step_perc$bootstrap_reps$Est)
  
  set.seed(aseed)
  
  step_basic <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial", 
      CI_type = "basic",
      R = 19
    )
  
  expect_s3_class(step_basic, "boot.selmodel")
  expect_identical(names(step_basic$est), c("estimator","param","Est","SE","bootstraps","basic_lower","basic_upper"))
  expect_identical(table(step_basic$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 19L)))
  expect_identical(step_basic$bootstrap_reps$Est, step_perc$bootstrap_reps$Est)
  
  set.seed(aseed)
  
  step_all <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial", 
      CI_type = c("large-sample","student","percentile","basic"),
      R = 19
    )
  
  expect_s3_class(step_all, "boot.selmodel")
  expect_identical(names(step_all$est), c("estimator","param","Est","SE","CI_lo","CI_hi","bootstraps","basic_lower","basic_upper","student_lower","student_upper","percentile_lower","percentile_upper"))
  expect_identical(table(step_all$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 19L)))
  expect_identical(step_all$bootstrap_reps, step_t$bootstrap_reps)
  
  expect_equal(
    dplyr::bind_cols(
      step_large$est, 
      step_basic$est[,c("bootstraps","basic_lower","basic_upper")], 
      step_t$est[,c("student_lower","student_upper")],
      step_perc$est[,c("percentile_lower","percentile_upper")]
    ), 
    step_all$est
  )
  
  set.seed(aseed)
  
  step_long <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial", 
      CI_type = c("large-sample","student","percentile","basic"),
      R = 19,
      format = "long"
    )
  
  expect_s3_class(step_long, "boot.selmodel")
  expect_identical(names(step_long$est), c("estimator","param","Est","SE","CI_lo","CI_hi","boot_CIs"))
  expect_identical(table(step_long$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 19L)))
  expect_identical(step_all$bootstrap_reps, step_long$bootstrap_reps)
  
  expect_equal(
    step_all$est[,c("estimator","param","Est","SE","CI_lo","CI_hi")],
    step_long$est[,c("estimator","param","Est","SE","CI_lo","CI_hi")]
  )
  
  long_CIs <- 
    dplyr::bind_rows(step_long$est$boot_CIs, .id = "param") |>
    dplyr::arrange(type, param) |>
    dplyr::select(type, param, bootstraps, lower, upper)
  
  all_CIs <- 
    dplyr::bind_rows(
      basic = dplyr::select(step_basic$est, param, bootstraps, lower = basic_lower, upper = basic_upper),
      percentile = dplyr::select(step_perc$est, param, bootstraps, lower = percentile_lower, upper = percentile_upper),
      student = dplyr::select(step_t$est, param, bootstraps, lower = student_lower, upper = student_upper),
      .id = "type"
    )
  rownames(all_CIs) <- NULL
  expect_equal(long_CIs, all_CIs)
  
  set.seed(aseed)
  
  step_none <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial", 
      CI_type = "none",
      R = 29
    )
  
  expect_s3_class(step_none, "boot.selmodel")
  expect_identical(step_large$est[,c("estimator","param","Est","SE")], step_none$est)
  expect_identical(names(step_none$est), c("estimator","param","Est","SE"))
  expect_identical(table(step_none$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 29L)))
  expect_identical(step_none$bootstrap_reps[1:(3*19),], step_t$bootstrap_reps[,c("param","Est")])
  
})

test_that("bootstrap_CI options for selection_model() work when bootstrap = 'exp'.", {
  
  aseed <- 20241029
  set.seed(aseed)
  
  step_large <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "exp", 
      CI_type = "large-sample",
      R = 24L
    )
  
  expect_s3_class(step_large, "selmodel")
  expect_identical(names(step_large$est), c("estimator","param","Est","SE","CI_lo","CI_hi"))
  expect_false(is.null(step_large$bootstrap_reps))
  
  set.seed(aseed)
  
  step_perc <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      CI_type = "percentile",
      bootstrap = "exp", 
      R = 24
    )
  
  expect_s3_class(step_perc, "boot.selmodel")
  expect_identical(names(step_perc$est), c("estimator","param","Est","SE","bootstraps","percentile_lower","percentile_upper"))
  expect_identical(table(step_perc$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 24L)))
  
  set.seed(aseed)
  
  step_t <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "exp", 
      CI_type = "student",
      R = 24
    )
  
  expect_s3_class(step_t, "boot.selmodel")
  expect_identical(names(step_t$est), c("estimator","param","Est","SE","bootstraps","student_lower","student_upper"))
  expect_identical(table(step_t$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 24L)))
  expect_identical(step_t$bootstrap_reps$Est, step_perc$bootstrap_reps$Est)
  
  set.seed(aseed)
  
  step_basic <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "exp", 
      CI_type = "basic",
      R = 24
    )
  
  expect_s3_class(step_basic, "boot.selmodel")
  expect_identical(names(step_basic$est), c("estimator","param","Est","SE","bootstraps","basic_lower","basic_upper"))
  expect_identical(table(step_basic$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 24L)))
  expect_identical(step_basic$bootstrap_reps$Est, step_perc$bootstrap_reps$Est)
  
  set.seed(aseed)
  
  step_all <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "exp", 
      CI_type = c("large-sample","student","percentile","basic"),
      R = 24
    )
  
  expect_s3_class(step_all, "boot.selmodel")
  expect_identical(names(step_all$est), c("estimator","param","Est","SE","CI_lo","CI_hi","bootstraps","basic_lower","basic_upper","student_lower","student_upper","percentile_lower","percentile_upper"))
  expect_identical(table(step_all$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 24L)))
  expect_identical(step_all$bootstrap_reps, step_t$bootstrap_reps)
  
  expect_equal(
    dplyr::bind_cols(
      step_large$est, 
      step_basic$est[,c("bootstraps","basic_lower","basic_upper")], 
      step_t$est[,c("student_lower","student_upper")],
      step_perc$est[,c("percentile_lower","percentile_upper")]
    ), 
    step_all$est
  )
  
  set.seed(aseed)
  
  step_long <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "exp", 
      CI_type = c("large-sample","student","percentile","basic"),
      R = 24,
      format = "long"
    )
  
  expect_s3_class(step_long, "boot.selmodel")
  expect_identical(names(step_long$est), c("estimator","param","Est","SE","CI_lo","CI_hi","boot_CIs"))
  expect_identical(table(step_long$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 24L)))
  expect_identical(step_all$bootstrap_reps, step_long$bootstrap_reps)
  
  expect_equal(
    step_all$est[,c("estimator","param","Est","SE","CI_lo","CI_hi")],
    step_long$est[,c("estimator","param","Est","SE","CI_lo","CI_hi")]
  )
  
  long_CIs <- 
    dplyr::bind_rows(step_long$est$boot_CIs, .id = "param") |>
    dplyr::arrange(type, param) |>
    dplyr::select(type, param, bootstraps, lower, upper)
  
  all_CIs <- 
    dplyr::bind_rows(
      basic = dplyr::select(step_basic$est, param, bootstraps, lower = basic_lower, upper = basic_upper),
      percentile = dplyr::select(step_perc$est, param, bootstraps, lower = percentile_lower, upper = percentile_upper),
      student = dplyr::select(step_t$est, param, bootstraps, lower = student_lower, upper = student_upper),
      .id = "type"
    )
  rownames(all_CIs) <- NULL
  expect_equal(long_CIs, all_CIs)
  
  set.seed(aseed)
  
  step_none <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "exp", 
      CI_type = "none",
      R = 37
    )
  
  expect_s3_class(step_none, "boot.selmodel")
  expect_identical(step_large$est[,c("estimator","param","Est","SE")], step_none$est)
  expect_identical(names(step_none$est), c("estimator","param","Est","SE"))
  expect_identical(table(step_none$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 37L)))
  expect_identical(step_none$bootstrap_reps[1:(3*24),], step_t$bootstrap_reps[,c("param","Est")])
  
})

test_that("CI_type options agree with simhelpers::bootstrap_CIs.", {

  suppressWarnings(
    step_multinomial <- 
      selection_model(
        data = dat,
        yi = d,
        sei = sd_d,
        pi = p_onesided,
        cluster = studyid,
        steps = 0.025,
        estimator = "ML",
        bootstrap = "multinomial", 
        CI_type = c("student","percentile","basic"),
        R = 39,
        seed = 20241031,
        format = "long"
      )
  )
  
  multi_boot <- get_boot_CIs(step_multinomial, CI_type = c("percentile","student","basic"), R = 39, format = "long")
  
  expect_equal(step_multinomial$est$boot_CIs, multi_boot)
  
  suppressWarnings(
    step_multinomial_multiR <- 
      selection_model(
        data = dat,
        yi = d,
        sei = sd_d,
        pi = p_onesided,
        cluster = studyid,
        steps = 0.025,
        estimator = "ML",
        bootstrap = "multinomial", 
        CI_type = c("student","percentile","basic"),
        R = 59,
        seed = 20240819,
        format = "long"
      )
  )
  
  multi_boot_multiR <- get_boot_CIs(step_multinomial_multiR, CI_type = c("percentile","student","basic"), R = 59, seed = 20240819, format = "long")
  
  expect_equal(step_multinomial_multiR$est$boot_CIs, multi_boot_multiR)

  suppressWarnings(
    step_exponential <- 
      selection_model(
        data = dat,
        yi = d,
        sei = sd_d,
        pi = p_onesided,
        cluster = studyid,
        steps = 0.025,
        estimator = "ML",
        bootstrap = "exp", 
        CI_type = c("student","percentile"),
        R = 49,
        seed = 20241101,
        format = "long"
      )
  )
  
  exp_boot <- get_boot_CIs(step_exponential, CI_type = c("percentile","student"), R = 49, format = "long")

  expect_equal(step_exponential$est$boot_CIs, exp_boot)
  
  
})

test_that("bootstrapping works with parallel processing.", {
  
  skip_on_cran()
  skip_if_not_installed("future")

  library(future)
  aseed <- 20241028
  
  plan(sequential)
  set.seed(aseed)
  
  step_sequential <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial", 
      CI_type = c("large-sample","percentile"),
      R = 40
    )
  
  
  plan(multisession, workers = 2L)
  set.seed(aseed)
  
  step_parallel <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial", 
      CI_type = c("large-sample","percentile"),
      R = 40
    )
  
  expect_identical(class(step_sequential), class(step_parallel))
  expect_identical(step_sequential$est, step_parallel$est)
  expect_identical(step_sequential$bootstrap_reps, step_parallel$bootstrap_reps)
  
})  


test_that("bootstrap CIs appear in the right order in models with predictors.", {
  
  set.seed(20241025)
  dat <- r_meta_categories(
    mean_smd = c(0, 0.2, 0.4),
    tau = c(0.1, 0.05, 0.13),
    omega = 0,
    m = 50,
    cor_mu = 0.3,
    cor_sd = 0.01,
    censor_fun = step_fun(cut_vals = .1, weights = 0.6),
    n_ES_sim = n_ES_param(mean_N = 80, mean_ES = 1L)
  )
  
  levels(dat$X) <- c("D","B","A")
  
  dat$sig <- dat$p_onesided < .1
  with(dat, table(X, sig))
  
  step_multinomial <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.1,
      mean_mods = ~ X,
      estimator = "hybrid",
      bootstrap = "multinomial", 
      CI_type = c("student","percentile","basic"),
      R = 59
    )
  
  # Check that Est falls within CI bounds
  expect_true(all(with(step_multinomial$est, basic_lower < Est & Est < basic_upper)))  
  expect_true(all(with(step_multinomial$est, student_lower < Est & Est < student_upper)))
  expect_true(all(with(step_multinomial$est, percentile_lower < Est & Est < percentile_upper)))

  # Check against get_boot_CIs  
  multi_boot <- get_boot_CIs(step_multinomial, CI_type = c("percentile","student","basic"), R = 59)
  multi_boot <- do.call(rbind, multi_boot)
  multi_boot <- multi_boot[step_multinomial$est$param,]
  expect_equal(
    subset(step_multinomial$est, select = bootstraps:percentile_upper), 
    multi_boot
  )
  
  step_multinomial <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.1,
      mean_mods = ~ X,
      var_mods = ~ 0 + X, 
      estimator = "hybrid",
      bootstrap = "multinomial", 
      CI_type = c("percentile","basic"),
      R = 39
    )
  
  # Check that Est falls within CI bounds
  expect_true(all(with(step_multinomial$est, basic_lower < Est & Est < basic_upper)))  
  expect_true(all(with(step_multinomial$est, percentile_lower < Est & Est < percentile_upper)))
  
  # Check against get_boot_CIs  
  multi_boot <- get_boot_CIs(step_multinomial, CI_type = c("percentile","basic"), R = 39)
  multi_boot <- do.call(rbind, multi_boot)
  multi_boot <- multi_boot[step_multinomial$est$param,]
  expect_equal(
    subset(step_multinomial$est, select = bootstraps:percentile_upper), 
    multi_boot
  )
  
  step_multinomial <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.1,
      mean_mods = ~ X,
      sel_mods = ~ 0 + X, 
      estimator = "hybrid",
      bootstrap = "multinomial", 
      CI_type = "percentile",
      seed = 20241031,
      R = 49 
    )
  
  # Check that Est falls within CI bounds
  expect_true(all(with(step_multinomial$est, percentile_lower < Est & Est < percentile_upper)))
  
  # Check against get_boot_CIs  
  multi_boot <- get_boot_CIs(step_multinomial, CI_type = "percentile", R = 49)
  multi_boot <- do.call(rbind, multi_boot)
  multi_boot <- multi_boot[step_multinomial$est$param,]
  
  expect_equal(
    subset(step_multinomial$est, select = bootstraps:percentile_upper), 
    multi_boot
  )
  
})  
