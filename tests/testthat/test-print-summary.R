set.seed(20240927)

dat <- r_meta_categories(
  mean_smd = c(0, 0.2, 0.4),
  tau = c(0.2, 0.05, 0.5),
  omega = 0,
  m = 20,
  cor_mu = 0.6,
  cor_sd = 0.01,
  censor_fun = step_fun(cut_vals = .025, weights = 0.3),
  n_ES_sim = n_ES_param(mean_N = 80, mean_ES = 3)
)



test_that("print() and summary() work for selmodel objects with no predictors.", {
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "CML",
    bootstrap = "none"
  )
  
  expect_output(print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","p_value","CI_lo","CI_hi"))
  expect_identical(mod_print[-1,1], c("beta","tau2","lambda1"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta","gamma","lambda1"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta","tau2","zeta1"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta","gamma","zeta1"))
  expect_identical(raw_params, row.names(mod$est))
  
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    steps = 0.025,
    estimator = "CML",
    vcov_type = "model-based",
    bootstrap = "none"
  )
  
  expect_output(print(mod))
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, row.names(mod$est))
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    selection_type = "beta",
    estimator = "CML",
    bootstrap = "none"
  )
  
  expect_output(print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","p_value","CI_lo","CI_hi"))
  expect_identical(mod_print[-1,1], c("beta","tau2","lambda1","lambda2"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta","gamma","lambda1","lambda2"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta","tau2","zeta1","zeta2"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta","gamma","zeta1","zeta2"))
  expect_identical(raw_params, row.names(mod$est))
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "ARGL",
    bootstrap = "multinomial",
    CI_type = "percentile",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","percentile_lower","percentile_upper"))
  expect_identical(mod_print[-1,1], c("beta","tau2","lambda1","lambda2"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta","gamma","lambda1","lambda2"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta","tau2","zeta1","zeta2"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta","gamma","zeta1","zeta2"))
  expect_identical(raw_params, row.names(mod$est))
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  check_selmodel_summary(mod, transf_gamma = FALSE)
  check_selmodel_summary(mod, transf_zeta = FALSE)
  check_selmodel_summary(mod, transf_gamma = FALSE, transf_zeta = FALSE)
  

})


test_that("print() and summary() work for selmodel objects with mean predictors.", {
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    mean_mods = ~ 0 + X, 
    estimator = "ARGL",
    bootstrap = "multinomial",
    CI_type = c("percentile","student"),
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","student_lower","student_upper","percentile_lower","percentile_upper"))
  expect_identical(mod_print[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","lambda1","lambda2"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta_XA","beta_XB","beta_XC","gamma","lambda1","lambda2"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","zeta1","zeta2"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta_XA","beta_XB","beta_XC","gamma","zeta1","zeta2"))
  expect_identical(raw_params, row.names(mod$est))
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    mean_mods = ~ 0 + X, 
    selection_type = "beta",
    estimator = "CML"
  )
  
  expect_output(print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","p_value","CI_lo","CI_hi"))
  expect_identical(mod_print[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","lambda1","lambda2"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta_XA","beta_XB","beta_XC","gamma","lambda1","lambda2"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","zeta1","zeta2"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta_XA","beta_XB","beta_XC","gamma","zeta1","zeta2"))
  expect_identical(raw_params, row.names(mod$est))
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  check_selmodel_summary(mod, transf_gamma = FALSE)
  check_selmodel_summary(mod, transf_zeta = FALSE)
  check_selmodel_summary(mod, transf_gamma = FALSE, transf_zeta = FALSE)
  
})

test_that("print() works for selmodel objects with variance predictors.", {
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    var_mods = ~ 0 + X, 
    estimator = "ARGL",
    bootstrap = "multinomial",
    CI_type = c("large-sample","basic"),
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","p_value","CI_lo","CI_hi","basic_lower","basic_upper"))
  expect_identical(mod_print[-1,1], c("beta","tau2_XA","tau2_XB","tau2_XC","lambda1"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta","gamma_XA","gamma_XB","gamma_XC","lambda1"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta","tau2_XA","tau2_XB","tau2_XC","zeta1"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta","gamma_XA","gamma_XB","gamma_XC","zeta1"))
  expect_identical(raw_params, row.names(mod$est))
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    mean_mods = ~ 0 + X,
    var_mods = ~ 0 + X, 
    estimator = "ARGL",
    bootstrap = "multinomial",
    CI_type = "student",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","student_lower","student_upper"))
  expect_identical(mod_print[-1,1], c("beta_XA","beta_XB","beta_XC","tau2_XA","tau2_XB","tau2_XC","lambda1"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta_XA","beta_XB","beta_XC","gamma_XA","gamma_XB","gamma_XC","lambda1"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2_XA","tau2_XB","tau2_XC","zeta1"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta_XA","beta_XB","beta_XC","gamma_XA","gamma_XB","gamma_XC","zeta1"))
  expect_identical(raw_params, row.names(mod$est))
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  check_selmodel_summary(mod, transf_gamma = FALSE)
  check_selmodel_summary(mod, transf_zeta = FALSE)
  check_selmodel_summary(mod, transf_gamma = FALSE, transf_zeta = FALSE)
  
})

test_that("print() works for selmodel objects with selection predictors.", {
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    sel_mods = ~ 0 + X, 
    estimator = "ARGL",
    bootstrap = "exponential",
    CI_type = "percentile",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","percentile_lower","percentile_upper"))
  expect_identical(mod_print[-1,1], c("beta","tau2","lambda1_XA","lambda1_XB","lambda1_XC"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta","gamma","lambda1_XA","lambda1_XB","lambda1_XC"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta","tau2","zeta1_XA","zeta1_XB","zeta1_XC"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta","gamma","zeta1_XA","zeta1_XB","zeta1_XC"))
  expect_identical(raw_params, row.names(mod$est))
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    mean_mods = ~ 0 + X,
    sel_mods = ~ 0 + X, 
    estimator = "ARGL",
    bootstrap = "exponential",
    CI_type = "student",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","student_lower","student_upper"))
  expect_identical(mod_print[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","lambda1_XA","lambda1_XB","lambda1_XC"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta_XA","beta_XB","beta_XC","gamma","lambda1_XA","lambda1_XB","lambda1_XC"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta_XA","beta_XB","beta_XC","gamma","zeta1_XA","zeta1_XB","zeta1_XC"))
  expect_identical(raw_params, row.names(mod$est))
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  check_selmodel_summary(mod, transf_gamma = FALSE)
  check_selmodel_summary(mod, transf_zeta = FALSE)
  check_selmodel_summary(mod, transf_gamma = FALSE, transf_zeta = FALSE)
  
})


test_that("print() works for selmodel objects with sel_zero predictors.", {
  
  dat$XisA <- as.integer(dat$X=="A")
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    sel_zero_mods = ~ XisA,
    sel_mods = ~ 0 + X, 
    estimator = "CML",
    bootstrap = "none"
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","p_value","CI_lo","CI_hi"))
  expect_identical(mod_print[-1,1], c("beta","tau2","lambda0_XisA","lambda1_XA","lambda1_XB","lambda1_XC"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta","gamma","lambda0_XisA","lambda1_XA","lambda1_XB","lambda1_XC"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta","tau2","zeta0_XisA","zeta1_XA","zeta1_XB","zeta1_XC"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta","gamma","zeta0_XisA","zeta1_XA","zeta1_XB","zeta1_XC"))
  expect_identical(raw_params, row.names(mod$est))

  expect_output(summary(mod))
  check_selmodel_summary(mod)
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    mean_mods = ~ 0 + X,
    sel_zero_mods = ~ XisA,
    sel_mods = ~ 0 + X, 
    estimator = "ARGL",
    bootstrap = "exponential",
    CI_type = "student",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","student_lower","student_upper"))
  expect_identical(mod_print[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","lambda0_XisA","lambda1_XA","lambda1_XB","lambda1_XC"))
  expect_identical(print_and_parse(mod, transf_gamma = FALSE)[-1,1], c("beta_XA","beta_XB","beta_XC","gamma","lambda0_XisA","lambda1_XA","lambda1_XB","lambda1_XC"))
  expect_identical(print_and_parse(mod, transf_zeta = FALSE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","zeta0_XisA","zeta1_XA","zeta1_XB","zeta1_XC"))
  raw_params <- print_and_parse(mod, transf_gamma = FALSE, transf_zeta = FALSE)[-1,1]
  expect_identical(raw_params, c("beta_XA","beta_XB","beta_XC","gamma","zeta0_XisA","zeta1_XA","zeta1_XB","zeta1_XC"))
  expect_identical(raw_params, row.names(mod$est))
  
  expect_output(summary(mod))
  check_selmodel_summary(mod)
  check_selmodel_summary(mod, transf_gamma = FALSE)
  check_selmodel_summary(mod, transf_zeta = FALSE)
  check_selmodel_summary(mod, transf_gamma = FALSE, transf_zeta = FALSE)
  
})

