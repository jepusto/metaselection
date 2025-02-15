test_that("calc_step_area() is correct.", {
  
  steps <- sort(runif(4))
  expect_equal(calc_step_area(rep(0,4), steps = steps), 1)
  
  zeta <- rnorm(4)  
  step_area <- calc_step_area(zeta = zeta, steps = steps)
  
  heights <- exp(zeta)
  widths <- diff(c(steps, 1))
  area_calc <- sum(heights * widths) / (1 - steps[1])
  
  expect_equal(step_area, area_calc)
  
})

test_that("calc_beta_area() is correct.", {
  
  steps <- c(runif(1, min = 0, max = 0.5), runif(1, min = 0.5, max = 1))
  expect_equal(calc_beta_area(c(0,0), steps = steps), 1)
  
  zeta1 <- rnorm(10)
  beta_areas <- sapply(zeta1, \(z) calc_beta_area(c(z, 0), steps = steps))
  calc_areas <- sapply(
    exp(zeta1), 
    \(l1) (diff(steps^l1) / l1 + (1 - steps[2]) * steps[2]^(l1 - 1)) / ((1 - steps[1]) * steps[1]^(l1 - 1))
  )
  expect_equal(beta_areas, calc_areas)
  
  zeta2 <- rnorm(10)
  beta_areas <- sapply(zeta2, \(z) calc_beta_area(c(0, z), steps = steps))
  calc_areas <- sapply(
    exp(zeta2), 
    \(l2) (diff(rev(1 - steps)^l2) / l2 + (1 - steps[2])^l2) / ((1 - steps[1])^l2)
  )
  expect_equal(beta_areas, calc_areas)
  
})

set.seed(20250215)

dat <- r_meta(
  mean_smd = 0, 
  tau = .1, omega = .01,
  m = 50, 
  cor_mu = .4, cor_sd = 0.001, 
  censor_fun = step_fun(cut_vals = c(.025,.50), weights = c(0.4,0.2)), 
  n_ES_sim = n_ES_param(40, 3)
)

test_that("p_area() works with 3-param step function models.", {
  
  step_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "CML",
    bootstrap = "none"
  )
  
  area_est <- p_area(step_noboot)
  expect_equal(area_est$Est, exp(step_noboot$est$Est[3]))
  expect_equal(dim(area_est), c(1L,2L))
  
  step_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "ARGL",
    bootstrap = "multinomial",
    CI_type = c("percentile","normal","basic","student","bias-corrected","BCa"),
    R = 19
  )
  
  area_est <- p_area(step_boot, CI_type = c("percentile","bias-corrected","BCa"))
  area_est$param <- NULL
  area_est$SE <- NULL
  area_est$bootstraps <- NULL
  expect_equal(area_est, exp(step_boot$est["zeta1",names(area_est)]), ignore_attr = TRUE)

  step_boot_nojack <- step_boot
  step_boot_nojack$jack_vals <- NULL
  
  expect_error(p_area(step_boot, CI_type = "student"))
  expect_error(p_area(step_boot_nojack, CI_type = "BCa"))
  
  expect_warning(
    all_boots <- p_area(step_boot)
  )
  
  expect_warning(
    perc_boots <- p_area(step_boot, CI_type = c("percentile","normal","student"))
  )

  expect_warning(
    BC_boots <- p_area(step_boot_nojack, CI_type = c("bias-corrected","BCa"))
  )
  
  BCa_boots <- p_area(step_boot, CI_type = c("bias-corrected","BCa"))

  all_boots_nowarn <- p_area(step_boot, warn = FALSE)
  expect_equal(all_boots, all_boots_nowarn)
  expect_equal(perc_boots, all_boots[,names(perc_boots)])
  expect_equal(BC_boots, all_boots[,names(BC_boots)])
  expect_equal(BCa_boots, all_boots[,names(BCa_boots)])
  
})

test_that("p_area() works with 4-param step function models.", {
  
  step_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,.500),
    estimator = "ARGL",
    bootstrap = "none"
  )
  area_calc <- sum(exp(step_noboot$est$Est[3:4]) * diff(c(step_noboot$steps, 1))) / (1 - step_noboot$steps[1])
  area_est <- p_area(step_noboot)
  expect_equal(area_est$Est, area_calc)
  expect_equal(dim(area_est), c(1L,2L))
  
  step_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,.500),
    estimator = "ARGL",
    bootstrap = "multinomial",
    CI_type = c("percentile","student","bias-corrected","BCa"),
    R = 12
  )
  step_boot_nojack <- step_boot
  step_boot_nojack$jack_vals <- NULL
  
  expect_error(p_area(step_boot, CI_type = "student"))
  expect_error(p_area(step_boot_nojack, CI_type = "BCa"))
  
  expect_warning(
    all_boots <- p_area(step_boot)
  )
  
  expect_warning(
    perc_boots <- p_area(step_boot, CI_type = c("percentile","student"))
  )
  
  expect_warning(
    BC_boots <- p_area(step_boot_nojack, CI_type = c("bias-corrected","BCa"))
  )
  
  BCa_boots <- p_area(step_boot, CI_type = c("bias-corrected","BCa"))
  
  all_boots_nowarn <- p_area(step_boot, warn = FALSE)
  expect_equal(all_boots, all_boots_nowarn)
  expect_equal(perc_boots, all_boots[,names(perc_boots)])
  expect_equal(BC_boots, all_boots[,names(BC_boots)])
  expect_equal(BCa_boots, all_boots[,names(BCa_boots)])
  
})

test_that("p_area() works with beta density models.", {
  
  beta_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    selection_type = "beta",
    steps = c(0.025,0.975),
    estimator = "CML",
    bootstrap = "none"
  )
  
  area_est <- p_area(beta_noboot)
  expect_equal(dim(area_est), c(1L,2L))
  
  suppressWarnings(
    beta_boot <- selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      cluster = studyid,
      steps = c(0.025,0.975),
      selection_type = "beta",
      estimator = "CML",
      bootstrap = "multinomial",
      CI_type = "percentile",
      R = 6
    )
  )
  
  expect_error(p_area(beta_boot, CI_type = "student"))
  expect_error(p_area(beta_boot, CI_type = "BCa"))
  
  perc_boots <- p_area(beta_boot)
  
  expect_warning(
    basic_boots <- p_area(beta_boot, CI_type = c("basic","student"))
  )
  
  expect_warning(
    BC_boots <- p_area(beta_boot, CI_type = c("bias-corrected","BCa"))
  )
  
  all_boots <- p_area(
    beta_boot, 
    CI_type = c("percentile","normal","basic","student","bias-corrected","BCa"), 
    warn = FALSE
  )
  
  expect_equal(perc_boots, all_boots[,names(perc_boots)])
  expect_equal(basic_boots, all_boots[,names(basic_boots)])
  expect_equal(BC_boots, all_boots[,names(BC_boots)])

})

