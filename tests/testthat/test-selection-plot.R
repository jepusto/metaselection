
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
  
  wt_noboot <- selection_wts(ML_noboot, pvals = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(ML_noboot)
  p_noboot
  

  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[3]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)
  
  wt_r_noboot <- selection_wts(ML_noboot, pvals = seq(0,1,length.out = 200), ref_pval = 0.9)
  p_r_noboot <- selection_plot(ML_noboot, ref_pval = 0.9)
  expect_identical(wt_r_noboot, p_r_noboot$data)
  
  p_noboot_sqrt <- selection_plot(ML_noboot, transform = "sqrt")
  expect_s3_class(p_noboot_sqrt, "ggplot")
  p_noboot_arcsin <- selection_plot(ML_noboot, transform = "asn")
  expect_s3_class(p_noboot_arcsin, "ggplot")
  
  suppressWarnings(
    ML_boot <- selection_model(
      data = dat,
      yi = d,
      sei = sd_d,
      cluster = studyid,
      steps = 0.025,
      estimator = "ML",
      bootstrap = "multinomial",
      CI_type = "percentile",
      R = 8
    )
  )
  
  expect_equal(ML_noboot$est$Est, ML_boot$est$Est)
  
  p_bootoff <- selection_plot(ML_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(ML_boot, pvals = seq(0,1,length.out = 200))
  p_boot <- selection_plot(ML_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[4]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[3]]$data)
  
  wt_r_boot <- selection_wts(ML_boot, pvals = seq(0,1,length.out = 200), ref_pval = 0.3)
  p_r_boot <- selection_plot(ML_boot, ref_pval = 0.3)
  expect_identical(wt_r_boot$wts, p_r_boot$data)
  expect_identical(wt_r_boot$boot_wts, p_r_boot$layers[[3]]$data)
  

  hybrid_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "hybrid",
    bootstrap = "none"
  )
  
  wt_noboot <- selection_wts(hybrid_noboot, pvals = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(hybrid_noboot)
  
  p_noboot
  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[3]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)
  
  hybrid_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "hybrid",
    bootstrap = "multinomial",
    CI_type = "percentile",
    R = 49
  )
  
  expect_equal(hybrid_noboot$est$Est, hybrid_boot$est$Est)
  
  
  p_bootoff <- selection_plot(hybrid_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(hybrid_boot, pvals = seq(0,1,length.out = 200))
  p_boot <- selection_plot(hybrid_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[4]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[3]]$data)
  

  hybrid_boot_default <- quick_boot_selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "hybrid",
    bootstrap = "multinomial",
    CI_type = "percentile"
  )
  
  p_boot_default <- selection_plot(hybrid_boot_default, draw_boots = FALSE)
  
  p_boot_default
  expect_s3_class(p_boot_default, "ggplot")
  expect_s3_class(p_boot_default$layers[[3]]$geom, "GeomLine")
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
  
  wt_noboot <- selection_wts(ML_noboot, pvals = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(ML_noboot)
  
  p_noboot
  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[3]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)
  
  wt_r_noboot <- selection_wts(ML_noboot, pvals = seq(0,1,length.out = 200), ref_pval = 0.47)
  p_r_noboot <- selection_plot(ML_noboot, ref_pval = 0.47)
  expect_identical(wt_r_noboot, p_r_noboot$data)
  
  p_noboot_sqrt <- selection_plot(ML_noboot, transform = "sqrt")
  expect_s3_class(p_noboot_sqrt, "ggplot")
  p_noboot_arcsin <- selection_plot(ML_noboot, transform = "asn")
  expect_s3_class(p_noboot_arcsin, "ggplot")
  
  
  ML_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "ML",
    bootstrap = "multinomial",
    CI_type = "percentile",
    R = 11
  )
  
  expect_equal(ML_noboot$est$Est, ML_boot$est$Est)
  
  
  p_bootoff <- selection_plot(ML_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(ML_boot, pvals = seq(0,1,length.out = 200))
  p_boot <- selection_plot(ML_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[4]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[3]]$data)
  
  wt_r_boot <- selection_wts(ML_boot, pvals = seq(0,1,length.out = 200), ref_pval = 0.34)
  p_r_boot <- selection_plot(ML_boot, ref_pval = 0.34)
  expect_identical(wt_r_boot$wts, p_r_boot$data)
  expect_identical(wt_r_boot$boot_wts, p_r_boot$layers[[3]]$data)
  
  
  ML_boot_default <- quick_boot_selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "ML",
    bootstrap = "multinomial",
    CI_type = "percentile"
  )
  
  p_boot_default <- selection_plot(ML_boot_default)
  
  p_boot_default
  expect_s3_class(p_boot_default, "ggplot")
  expect_s3_class(p_boot_default$layers[[3]]$geom, "GeomLine")
  
  
  hybrid_noboot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "hybrid",
    bootstrap = "none"
  )
  
  wt_noboot <- selection_wts(hybrid_noboot, pvals = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(hybrid_noboot)
  
  p_noboot
  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[3]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)
  
  hybrid_boot <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "hybrid",
    bootstrap = "multinomial",
    CI_type = "percentile",
    R = 39
  )
  
  expect_equal(hybrid_noboot$est$Est, hybrid_boot$est$Est)
  
  p_bootoff <- selection_plot(hybrid_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(hybrid_boot, pvals = seq(0,1,length.out = 200))
  p_boot <- selection_plot(hybrid_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[4]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[3]]$data)
  
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
  
  wt_noboot <- selection_wts(ML_noboot, pvals = seq(0,1,length.out = 200))
  p_noboot <- selection_plot(ML_noboot)
  
  p_noboot
  expect_s3_class(p_noboot, "ggplot")
  expect_s3_class(p_noboot$layers[[3]]$geom, "GeomArea")
  expect_identical(wt_noboot, p_noboot$data)

  wt_r_noboot <- selection_wts(ML_noboot, pvals = seq(0,1,length.out = 200), ref_pval = 0.6)
  p_r_noboot <- selection_plot(ML_noboot, ref_pval = 0.6)
  expect_identical(wt_r_noboot, p_r_noboot$data)
  
  
  p_noboot_sqrt <- selection_plot(ML_noboot, transform = "sqrt")
  expect_s3_class(p_noboot_sqrt, "ggplot")
  p_noboot_arcsin <- selection_plot(ML_noboot, transform = "asn")
  expect_s3_class(p_noboot_arcsin, "ggplot")
  
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
      CI_type = "percentile",
      R = 6
    )
  )
  
  expect_equal(ML_noboot$est$Est, ML_boot$est$Est)
  
  p_bootoff <- selection_plot(ML_boot, draw_boots = FALSE)
  
  p_bootoff
  expect_s3_class(p_bootoff, "ggplot")
  expect_s3_class(p_bootoff$layers[[3]]$geom, "GeomLine")
  expect_identical(wt_noboot, p_bootoff$data)
  
  wt_boot <- selection_wts(ML_boot, pvals = seq(0,1,length.out = 200))
  p_boot <- selection_plot(ML_boot)
  
  p_boot
  expect_s3_class(p_boot, "ggplot")
  expect_s3_class(p_boot$layers[[3]]$geom, "GeomLine")
  expect_s3_class(p_boot$layers[[4]]$geom, "GeomLine")
  expect_identical(wt_boot$wts, p_boot$data)
  expect_identical(wt_boot$boot_wts, p_boot$layers[[3]]$data)

  wt_r_boot <- selection_wts(ML_boot, pvals = seq(0,1,length.out = 200), ref_pval = 0.25)
  p_r_boot <- selection_plot(ML_boot, ref_pval = 0.25)

  expect_identical(wt_r_boot$wts, p_r_boot$data)
  expect_identical(wt_r_boot$boot_wts, p_r_boot$layers[[3]]$data)
  
})

