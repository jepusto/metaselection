
test_that("selection_plot() works for 3PSM.", {

  set.seed(20240911)
  
  dat <- r_meta(
    mean_smd = 0, 
    tau = .1, omega = .01,
    m = 50, 
    cor_mu = .4, cor_sd = 0.001, 
    censor_fun = step_fun(cut_vals = .025, weights = 0.4), 
    n_ES_sim = n_ES_param(40, 3)
  )
  
  ML_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "ML",
    bootstrap = "none"
  )
  
  wt_noboot <- selection_wts(ML_noboot, pts = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(ML_noboot)
  p_noboot
  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[2]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)
  
  ML_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "ML",
    bootstrap = "multinomial",
    boot_CI = "percentile",
    R = 29
  )
  
  expect_equal(ML_noboot$est$Est, ML_boot$est$Est)
  
  p_bootoff <- selection_plot(ML_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[2]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(ML_boot, pts = seq(0,1,length.out = 200))
  p_boot <- selection_plot(ML_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[2]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[2]]$data)
  

  hybrid_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "hybrid",
    bootstrap = "none"
  )
  
  wt_noboot <- selection_wts(hybrid_noboot, pts = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(hybrid_noboot)
  
  p_noboot
  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[2]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)
  
  hybrid_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "hybrid",
    bootstrap = "multinomial",
    boot_CI = "percentile",
    R = 199
  )
  
  expect_equal(hybrid_noboot$est$Est, hybrid_boot$est$Est)
  
  
  p_bootoff <- selection_plot(hybrid_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[2]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(hybrid_boot, pts = seq(0,1,length.out = 200))
  p_boot <- selection_plot(hybrid_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[2]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[2]]$data)
  
})

test_that("selection_plot() works for 4PSM.", {
  
  set.seed(20240913)
  
  dat <- r_meta(
    mean_smd = 0, 
    tau = .1, omega = .01,
    m = 60, 
    cor_mu = .4, cor_sd = 0.001, 
    censor_fun = step_fun(cut_vals = c(.025, .500), weights = c(0.3, 0.1)), 
    n_ES_sim = n_ES_param(40, 3)
  )
  
 
  ML_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "ML",
    bootstrap = "none"
  )
  
  wt_noboot <- selection_wts(ML_noboot, pts = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(ML_noboot)
  
  p_noboot
  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[2]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)
  
  ML_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "ML",
    bootstrap = "multinomial",
    boot_CI = "percentile",
    R = 29
  )
  
  expect_equal(ML_noboot$est$Est, ML_boot$est$Est)
  
  
  p_bootoff <- selection_plot(ML_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[2]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(ML_boot, pts = seq(0,1,length.out = 200))
  p_boot <- selection_plot(ML_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[2]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[2]]$data)
  
  
  hybrid_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "hybrid",
    bootstrap = "none"
  )
  
  wt_noboot <- selection_wts(hybrid_noboot, pts = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(hybrid_noboot)
  
  p_noboot
  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[2]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)
  
  hybrid_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "hybrid",
    bootstrap = "multinomial",
    boot_CI = "percentile",
    R = 199
  )
  
  expect_equal(hybrid_noboot$est$Est, hybrid_boot$est$Est)
  
  p_bootoff <- selection_plot(hybrid_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[2]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(hybrid_boot, pts = seq(0,1,length.out = 200))
  p_boot <- selection_plot(hybrid_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[2]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[2]]$data)
  
})


test_that("selection_plot() works for beta model", {
  
  skip_on_cran()
  
  set.seed(20240914)
  
  dat <- r_meta(
    mean_smd = 0, 
    tau = .1, omega = .01,
    m = 60, 
    cor_mu = .4, cor_sd = 0.001, 
    censor_fun = beta_fun(delta_1 = 0.2, delta_2 = 0.9), 
    n_ES_sim = n_ES_param(40, 3)
  )
  
  
  ML_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    selection_type = "beta",
    steps = c(0.025,0.975),
    estimator = "ML",
    bootstrap = "none"
  )
  
  wt_noboot <- selection_wts(ML_noboot, pts = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(ML_noboot)
  
  p_noboot
  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[2]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)
  
  suppressWarnings(
    ML_boot <- selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      cluster = studyid,
      steps = c(0.025,0.975),
      selection_type = "beta",
      estimator = "ML",
      bootstrap = "multinomial",
      boot_CI = "percentile",
      R = 19
    )
  )
  
  expect_equal(ML_noboot$est$Est, ML_boot$est$Est)
  
  p_bootoff <- selection_plot(ML_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[2]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(ML_boot, pts = seq(0,1,length.out = 200))
  p_boot <- selection_plot(ML_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[2]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[2]]$data)
  
})

