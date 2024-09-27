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


print_and_parse <- function(mod,...) {
  p <- capture_output_lines(print(mod,...))
  d <- do.call(rbind, strsplit(p, " +"))
  d[,-1]
}

pull_argument <- function(s, str) {
  x <- s[grepl(str, s)]
  substr(x, nchar(str) + 2L, nchar(x))
}

summary_and_parse <- function(mod, ...) {
  s <- capture_output_lines(summary(mod, ...))
  mod <- s[1]
  steps <- pull_argument(s, "Steps:")
  estimator <- pull_argument(s,"Estimator:",) 
  boot_type <- pull_argument(s, "Bootstrap type:")
  R <- pull_argument(s, "Number of replications:")
  CI_type <- pull_argument(s, "CI type:")
  
  headers <- which(grepl("estimates:", p))
  
}

test_that("print() and summary() work for selmodel objects with no predictors.", {
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = 0.025,
    estimator = "ML",
    bootstrap = "none"
  )
  
  expect_output(print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","CI_lo","CI_hi"))
  expect_identical(mod_print[-1,1], c("beta","gamma","zeta1"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE)[-1,1], c("beta","tau2","zeta1"))
  expect_identical(print_and_parse(mod, transf_zeta = TRUE)[-1,1], c("beta","gamma","lambda_1"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE, transf_zeta = TRUE)[-1,1], c("beta","tau2","lambda_1"))
  
  expect_output(summary(mod))
  
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    estimator = "hybrid",
    bootstrap = "multinomial",
    boot_CI = "percentile",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","percentile_lower","percentile_upper"))
  expect_identical(mod_print[-1,1], c("beta","gamma","zeta1","zeta2"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE)[-1,1], c("beta","tau2","zeta1","zeta2"))
  expect_identical(print_and_parse(mod, transf_zeta = TRUE)[-1,1], c("beta","gamma","lambda_1","lambda_2"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE, transf_zeta = TRUE)[-1,1], c("beta","tau2","lambda_1","lambda_2"))
  

})


test_that("print() works for selmodel objects with mean predictors.", {
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.025,0.500),
    mean_mods = ~ 0 + X, 
    estimator = "hybrid",
    bootstrap = "multinomial",
    boot_CI = "percentile",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","percentile_lower","percentile_upper"))
  expect_identical(mod_print[-1,1], c("beta_XA","beta_XB","beta_XC","gamma","zeta1","zeta2"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","zeta1","zeta2"))
  expect_identical(print_and_parse(mod, transf_zeta = TRUE)[-1,1], c("beta_XA","beta_XB","beta_XC","gamma","lambda_1","lambda_2"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE, transf_zeta = TRUE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","lambda_1","lambda_2"))
  
  
})

test_that("print() works for selmodel objects with variance predictors.", {
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    var_mods = ~ 0 + X, 
    estimator = "hybrid",
    bootstrap = "multinomial",
    boot_CI = "basic",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","basic_lower","basic_upper"))
  expect_identical(mod_print[-1,1], c("beta","gamma_XA","gamma_XB","gamma_XC","zeta1"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE)[-1,1], c("beta","tau2_XA","tau2_XB","tau2_XC","zeta1"))
  expect_identical(print_and_parse(mod, transf_zeta = TRUE)[-1,1], c("beta","gamma_XA","gamma_XB","gamma_XC","lambda_1"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE, transf_zeta = TRUE)[-1,1], c("beta","tau2_XA","tau2_XB","tau2_XC","lambda_1"))
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    mean_mods = ~ 0 + X,
    var_mods = ~ 0 + X, 
    estimator = "hybrid",
    bootstrap = "multinomial",
    boot_CI = "student",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","student_lower","student_upper"))
  expect_identical(mod_print[-1,1], c("beta_XA","beta_XB","beta_XC","gamma_XA","gamma_XB","gamma_XC","zeta1"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2_XA","tau2_XB","tau2_XC","zeta1"))
  expect_identical(print_and_parse(mod, transf_zeta = TRUE)[-1,1], c("beta_XA","beta_XB","beta_XC","gamma_XA","gamma_XB","gamma_XC","lambda_1"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE, transf_zeta = TRUE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2_XA","tau2_XB","tau2_XC","lambda_1"))
  
})

test_that("print() works for selmodel objects with selection predictors.", {
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    sel_mods = ~ 0 + X, 
    estimator = "hybrid",
    bootstrap = "exponential",
    boot_CI = "percentile",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","percentile_lower","percentile_upper"))
  expect_identical(mod_print[-1,1], c("beta","gamma","zeta1_XA","zeta1_XB","zeta1_XC"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE)[-1,1], c("beta","tau2","zeta1_XA","zeta1_XB","zeta1_XC"))
  expect_identical(print_and_parse(mod, transf_zeta = TRUE)[-1,1], c("beta","gamma","lambda_1_XA","lambda_1_XB","lambda_1_XC"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE, transf_zeta = TRUE)[-1,1], c("beta","tau2","lambda_1_XA","lambda_1_XB","lambda_1_XC"))
  
  mod <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    cluster = studyid,
    steps = c(0.05),
    mean_mods = ~ 0 + X,
    sel_mods = ~ 0 + X, 
    estimator = "hybrid",
    bootstrap = "exponential",
    boot_CI = "student",
    R = 9
  )
  
  expect_output(mod_print <- print(mod))
  mod_print <- print_and_parse(mod)
  expect_identical(mod_print[1,], c("param","Est","SE","student_lower","student_upper"))
  expect_identical(mod_print[-1,1], c("beta_XA","beta_XB","beta_XC","gamma","zeta1_XA","zeta1_XB","zeta1_XC"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","zeta1_XA","zeta1_XB","zeta1_XC"))
  expect_identical(print_and_parse(mod, transf_zeta = TRUE)[-1,1], c("beta_XA","beta_XB","beta_XC","gamma","lambda_1_XA","lambda_1_XB","lambda_1_XC"))
  expect_identical(print_and_parse(mod, transf_gamma = TRUE, transf_zeta = TRUE)[-1,1], c("beta_XA","beta_XB","beta_XC","tau2","lambda_1_XA","lambda_1_XB","lambda_1_XC"))
  
})