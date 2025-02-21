
skip()

set.seed(20250221)

dat <- r_meta(
  mean_smd = 0.8,
  tau = .05,
  omega = .00,
  m = 40,
  cor_mu = .8,
  cor_sd = 0.001,
  censor_fun = step_fun(cut_vals = .025, weights = 0.05),
  n_ES_sim = n_ES_param(40, 3)
)

selmod <- selection_model(
  data = dat,
  yi = d,
  sei = sd_d,
  pi = p_onesided,
  cluster = studyid,
  steps = 0.025,
  estimator = "CML",
  bootstrap = "multinomial",
  CI_type = "percentile",
  R = 399L,
  retry_bootstrap = 3L
)

selmod$est

