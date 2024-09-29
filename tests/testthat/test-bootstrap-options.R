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
  step_RVE <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "none"
    )
  
  step_RVE_A <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "none",
      boot_CI = "percentile"
    )
  
  expect_identical(step_RVE_A$est, step_RVE$est)
  expect_identical(class(step_RVE_A), class(step_RVE))
  
  step_RVE_B <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "none",
      boot_CI = "percentile",
      R = 999
    )
  
  expect_identical(step_RVE_B$est, step_RVE$est)
  expect_identical(class(step_RVE_B), class(step_RVE))
  
  step_RVE_C <- 
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
  
  expect_identical(step_RVE_C$est, step_RVE$est)
  expect_identical(class(step_RVE_C), class(step_RVE))
  
  step_RVE_D <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "exp",
      boot_CI = "percentile",
      R = 0
    )
  
  expect_identical(step_RVE_D$est, step_RVE$est)
  expect_identical(class(step_RVE_D), class(step_RVE))
})  

test_that("bootstrap_CI options for selection_model() work when bootstrap = 'multinomial'.", {
  
  
  set.seed(20240807)
  
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
      boot_CI = "large-sample",
      R = 19
    )
  
  expect_s3_class(step_large, "selmodel")
  expect_identical(names(step_large$est), c("estimator","param","Est","SE","CI_lo","CI_hi"))
  expect_false(is.null(step_large$bootstrap_reps))

  set.seed(20240807)
  
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
      R = 19
    )
  
  expect_s3_class(step_perc, "boot.selmodel")
  expect_identical(names(step_perc$est), c("estimator","param","Est","SE","bootstraps","percentile_lower","percentile_upper"))
  expect_identical(table(step_perc$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 19L)))
  
  set.seed(20240807)
  
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
      boot_CI = "student",
      R = 19
    )
  
  expect_s3_class(step_t, "boot.selmodel")
  expect_identical(names(step_t$est), c("estimator","param","Est","SE","bootstraps","student_lower","student_upper"))
  expect_identical(table(step_t$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 19L)))
  expect_identical(step_t$bootstrap_reps$Est, step_perc$bootstrap_reps$Est)
  
  set.seed(20240807)
  
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
      boot_CI = "basic",
      R = 19
    )
  
  expect_s3_class(step_basic, "boot.selmodel")
  expect_identical(names(step_basic$est), c("estimator","param","Est","SE","bootstraps","basic_lower","basic_upper"))
  expect_identical(table(step_basic$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 19L)))
  expect_identical(step_basic$bootstrap_reps$Est, step_perc$bootstrap_reps$Est)
  
  set.seed(20240807)
  
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
      boot_CI = c("large-sample","student","percentile","basic"),
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
  
  set.seed(20240807)
  
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
      boot_CI = c("large-sample","student","percentile","basic"),
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
  expect_equal(long_CIs, all_CIs)
  
  set.seed(20240807)
  
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
      boot_CI = "none",
      R = 29
    )
  
  expect_s3_class(step_none, "boot.selmodel")
  expect_identical(step_large$est[,c("estimator","param","Est","SE")], step_none$est)
  expect_identical(names(step_none$est), c("estimator","param","Est","SE"))
  expect_identical(table(step_none$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 29L)))
  expect_identical(step_none$bootstrap_reps[1:(3*19),], step_t$bootstrap_reps[,c("param","Est")])
  
})

test_that("bootstrap_CI options for selection_model() work when bootstrap = 'exp'.", {
  
  
  set.seed(20240808)
  
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
      boot_CI = "large-sample",
      R = 24L
    )
  
  expect_s3_class(step_large, "selmodel")
  expect_identical(names(step_large$est), c("estimator","param","Est","SE","CI_lo","CI_hi"))
  expect_false(is.null(step_large$bootstrap_reps))
  
  set.seed(20240808)
  
  step_perc <- 
    selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      pi = p_onesided,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "exp", 
      R = 24
    )
  
  expect_s3_class(step_perc, "boot.selmodel")
  expect_identical(names(step_perc$est), c("estimator","param","Est","SE","bootstraps","percentile_lower","percentile_upper"))
  expect_identical(table(step_perc$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 24L)))
  
  set.seed(20240808)
  
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
      boot_CI = "student",
      R = 24
    )
  
  expect_s3_class(step_t, "boot.selmodel")
  expect_identical(names(step_t$est), c("estimator","param","Est","SE","bootstraps","student_lower","student_upper"))
  expect_identical(table(step_t$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 24L)))
  expect_identical(step_t$bootstrap_reps$Est, step_perc$bootstrap_reps$Est)
  
  set.seed(20240808)
  
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
      boot_CI = "basic",
      R = 24
    )
  
  expect_s3_class(step_basic, "boot.selmodel")
  expect_identical(names(step_basic$est), c("estimator","param","Est","SE","bootstraps","basic_lower","basic_upper"))
  expect_identical(table(step_basic$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 24L)))
  expect_identical(step_basic$bootstrap_reps$Est, step_perc$bootstrap_reps$Est)
  
  set.seed(20240808)
  
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
      boot_CI = c("large-sample","student","percentile","basic"),
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
  
  set.seed(20240808)
  
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
      boot_CI = c("large-sample","student","percentile","basic"),
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
  expect_equal(long_CIs, all_CIs)
  
  set.seed(20240808)
  
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
      boot_CI = "none",
      R = 37
    )
  
  expect_s3_class(step_none, "boot.selmodel")
  expect_identical(step_large$est[,c("estimator","param","Est","SE")], step_none$est)
  expect_identical(names(step_none$est), c("estimator","param","Est","SE"))
  expect_identical(table(step_none$bootstrap_reps$param), table(rep(c("beta","gamma","zeta1"), 37L)))
  expect_identical(step_none$bootstrap_reps[1:(3*24),], step_t$bootstrap_reps[,c("param","Est")])
  
})

test_that("boot_CI options agree with simhelpers::bootstrap_CIs.", {

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
        boot_CI = c("student","percentile","basic"),
        R = 39,
        format = "long"
      )
  )
  
  multi_boot <- get_boot_CIs(step_multinomial, type = c("percentile","student","basic"), format = "long")
  
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
        boot_CI = c("student","percentile","basic"),
        R = c(39,59,79,99),
        seed = 20240819,
        format = "long"
      )
  )
  
  multi_boot_multiR <- get_boot_CIs(step_multinomial_multiR, type = c("percentile","student","basic"), B_vals = c(39, 59, 79, 99), seed = 20240819, format = "long")
  
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
        boot_CI = c("student","percentile"),
        R = 49,
        format = "long"
      )
  )
  
  exp_boot <- get_boot_CIs(step_exponential, type = c("percentile","student"), format = "long")

  expect_equal(step_exponential$est$boot_CIs, exp_boot)
  
  
})

test_that("bootstrapping works with parallel processing.", {
  
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  
  library(future)
  
  plan(sequential)
  set.seed(20240916)
  
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
      boot_CI = c("large-sample","percentile"),
      R = 40
    )
  
  
  plan(multisession, workers = 2L)
  set.seed(20240916)
  
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
      boot_CI = c("large-sample","percentile"),
      R = 40
    )
  
  expect_identical(class(step_sequential), class(step_parallel))
  
})  
