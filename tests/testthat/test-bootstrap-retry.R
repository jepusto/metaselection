skip_on_cran()

data("wwc_es")
wwc_es$n <- round(wwc_es$n / 3)
n_ES_sim <- n_ES_empirical(wwc_es)

test_that("retry_bootstrap argument works correctly.",{
  
  set.seed(20251115L)
  
  dat <- r_meta(mean_smd = 0.8,
                tau = 0.15,
                omega = 0,
                m = 60,
                cor_mu = 0.8,
                cor_sd = 0.05,
                censor_fun =  step_fun(cut_vals = .025, weights = 0.05),  
                n_ES_sim = n_ES_sim)
  
  selmod0 <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    pi = p_onesided,
    cluster = studyid,
    steps = 0.025,
    priors = NULL,
    estimator = "CML",
    bootstrap = "multinomial",
    CI_type = c("percentile","student"),
    R = 199L,
    retry_bootstrap = 0L
  )
  
  expect_lt(unique(selmod0$est$bootstraps), 199L)
  
  set.seed(20251115L)
  
  dat <- r_meta(mean_smd = 0.8,
                tau = 0.15,
                omega = 0,
                m = 60,
                cor_mu = 0.8,
                cor_sd = 0.05,
                censor_fun =  step_fun(cut_vals = .025, weights = 0.05),  
                n_ES_sim = n_ES_sim)
  
  selmod3 <- selection_model(
    data = dat,
    yi = d,
    sei = sd_d,
    pi = p_onesided,
    cluster = studyid,
    steps = 0.025,
    priors = NULL,
    estimator = "CML",
    bootstrap = "multinomial",
    CI_type = c("percentile","student"),
    R = 199L,
    retry_bootstrap = 3L
  )
  
  expect_equal(unique(selmod3$est$bootstraps), 199L)
  
})
