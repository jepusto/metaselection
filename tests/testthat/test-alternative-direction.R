skip_if_not_installed("metadat")
data("dat.hackshaw1998", package = "metadat")
dat <- dat.hackshaw1998
dat$yi_neg <- -1 * dat$yi
dat$sei <- sqrt(dat$vi)
dat$year_c <- (dat$year - 1990L) / 10

test_that("CML step models are consistent when alternative = 'less'.", {

  steps <- c(.05, .10, .50)
  test_mu_gamma <- c(0, log(0.1) / 2)
  test_zeta <- c(-0.05, 0.3, -0.6)
  test_chi <- c(rev(test_zeta[-length(test_zeta)]), 0) - test_zeta[length(test_zeta)]
  test_param <- c(test_mu_gamma, test_zeta)
  test_param_trans <- c(test_mu_gamma, test_chi)
  
  # equivalence of log-likelihoods
  ll_gt <- step_loglik(
    test_param, yi = dat$yi, sei = dat$sei,
    steps = steps
  )
  ll_ls <- step_loglik(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    steps = rev(1 - steps),
    Hsgn = -1L
  )
  expect_equal(ll_gt, ll_ls)

  # equivalence of scores
  score_gt <- step_score(
    test_param, yi = dat$yi, sei = dat$sei,
    steps = steps
  )
  score_ls <- step_score(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    steps = rev(1 - steps),
    Hsgn = -1L
  )
  expect_equal(score_gt[1:2], score_ls[1:2])
  expect_equal(
    as.vector(score_gt[2+1:length(steps)]), 
    as.vector(c(rev(score_ls[2+1:(length(steps) - 1L)]), -sum(score_ls[2+1:length(steps)])))
  )
  

  # equivalence of Hessians
  Hess_gt <- step_hessian(
    test_param, yi = dat$yi, sei = dat$sei,
    steps = steps
  )
  Hess_ls <- step_hessian(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    steps = rev(1 - steps),
    Hsgn = -1L
  )
  dimnames(Hess_gt) <- NULL
  dimnames(Hess_ls) <- NULL
  expect_equal(Hess_gt[1:2,1:2], Hess_ls[1:2,1:2])
  expect_equal(
    Hess_gt[1:2,2+1:length(steps)], 
    cbind(Hess_ls[1:2,2+(length(steps) - 1L):1], -rowSums(Hess_ls[1:2,2+1:length(steps)])),
    ignore_attr = TRUE
  )

  
  # equivalence of parameter estimates
  check_valence_equivalence(
    data = dat, 
    yi = yi, yi_neg = yi_neg, sei = sei,
    steps = c(.05, .10, .50),
    priors = NULL
  )
  
  check_valence_equivalence(
    data = dat, 
    yi = yi, yi_neg = yi_neg, sei = sei,
    mean_mods = ~ design + year_c,
    steps = c(.025, .50),
    priors = NULL,
    vcov_type = "model-based"
  )
  
  check_valence_equivalence(
    data = dat, 
    yi = yi, yi_neg = yi_neg, sei = sei,
    steps = c(.025),
    priors = NULL,
    bootstrap = "multinomial",
    CI_type = c("large-sample","percentile","normal","basic","bias-corrected","student"),
    R = 19L
  )
  
  
})


test_that("CML beta models are consistent when alternative = 'less' and steps are symmetric.", {
  
  steps <- c(.025, .975)
  test_mu_gamma <- c(0, log(0.1) / 2)
  test_zeta <- c(-0.9, -0.1)
  test_param <- c(test_mu_gamma, test_zeta)
  test_param_trans <- c(test_mu_gamma, rev(test_zeta))
  
  # equivalence of calculated parameters
  params_gt <- parse_beta_params(
    test_param, yi = dat$yi, sei = dat$sei,
    alpha = steps
  )
  params_ls <- parse_beta_params(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    alpha = steps,
    Hsgn = -1L
  )
  expect_equal(
    params_gt[c("k","x_dim","u_dim","beta","gamma","mu","tausq","eta","alpha","weight_vec","H_names")],
    params_ls[c("k","x_dim","u_dim","beta","gamma","mu","tausq","eta","alpha","weight_vec","H_names")]
  )
  expect_equal(params_gt$zeta, rev(params_ls$zeta))
  expect_equal(params_gt$lambda, rev(params_ls$lambda))
  expect_equal(params_gt$alpha_lambda, rev(params_ls$alpha_lambda))
  expect_equal(params_gt$pi_tilde, 1 - params_ls$pi_tilde)
  
  # equivalence of log-likelihoods
  ll_gt <- beta_loglik(
    test_param, yi = dat$yi, sei = dat$sei,
    steps = steps
  )
  ll_ls <- beta_loglik(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    steps = steps,
    Hsgn = -1L
  )
  expect_equal(ll_gt, ll_ls)
  
  # equivalence of scores
  score_gt <- beta_score(
    test_param, yi = dat$yi, sei = dat$sei,
    steps = steps,
  )
  score_ls <- beta_score(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    steps = rev(1 - steps),
    Hsgn = -1L
  )
  expect_equal(score_gt, score_ls[c(1,2,4,3)], ignore_attr = TRUE)
  
  
  # equivalence of Hessians
  Hess_gt <- beta_hessian(
    test_param, yi = dat$yi, sei = dat$sei,
    steps = steps
  )
  Hess_ls <- beta_hessian(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    steps = rev(1 - steps),
    Hsgn = -1L
  )
  expect_equal(Hess_gt, Hess_ls[c(1,2,4,3),c(1,2,4,3)], ignore_attr = TRUE)
  
  
  # equivalence of parameter estimates
  check_valence_equivalence(
    data = dat, 
    yi = yi, yi_neg = yi_neg, sei = sei,
    selection_type = "beta",
    steps = steps,
    priors = NULL
  )
  
  check_valence_equivalence(
    data = dat, 
    yi = yi, yi_neg = yi_neg, sei = sei,
    mean_mods = ~ design + year_c,
    selection_type = "beta",
    steps = steps,
    priors = NULL,
    vcov_type = "model-based"
  )
  

})


test_that("CML beta models are consistent when alternative = 'less' and steps are symmetric.", {

  steps <- c(.05, .50)
  test_mu_gamma <- c(0.2, log(0.02^2))
  test_zeta <- c(-0.2, 0.1)
  test_param <- c(test_mu_gamma, test_zeta)
  test_param_trans <- c(test_mu_gamma, rev(test_zeta))
  
  # equivalence of calculated parameters
  params_gt <- parse_beta_params(
    test_param, yi = dat$yi, sei = dat$sei,
    alpha = steps,
    Hsgn = 1L,
    calc_Ai = TRUE,
    calc_Ai_deriv = TRUE
  )
  params_ls <- parse_beta_params(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    alpha = rev(1 - steps),
    Hsgn = -1L,
    calc_Ai = TRUE,
    calc_Ai_deriv = TRUE
  )
  
  
  expect_equal(
    params_gt[c("k","x_dim","u_dim","beta","gamma","mu","tausq","eta","weight_vec","H_names","Ai")],
    params_ls[c("k","x_dim","u_dim","beta","gamma","mu","tausq","eta","weight_vec","H_names","Ai")]
  )
  
  expect_equal(params_gt$zeta, rev(params_ls$zeta))
  expect_equal(params_gt$lambda, rev(params_ls$lambda))
  expect_equal(params_gt$alpha_lambda, rev(params_ls$alpha_lambda))
  expect_equal(params_gt$g_dot, rev(params_ls$g_dot))
  expect_equal(params_gt$pi_tilde, 1 - params_ls$pi_tilde)
  
  expect_equal(params_gt$c_1ij, -params_ls$c_2ij)
  expect_equal(params_gt$c_2ij, -params_ls$c_1ij)
  expect_equal(params_gt$B_0ij, params_ls$B_2ij)
  expect_equal(params_gt$B_2ij, params_ls$B_0ij)
  
  expect_equal(params_gt$dA_dmu, params_ls$dA_dmu)
  expect_equal(params_gt$dA_deta, params_ls$dA_deta)
  expect_equal(params_gt$dA_dlambda1, params_ls$dA_dlambda2)
  expect_equal(params_gt$dA_dlambda1, params_ls$dA_dlambda2)
  
  
  # equivalence of log-likelihoods
  ll_gt <- beta_loglik(
    test_param, yi = dat$yi, sei = dat$sei,
    steps = steps
  )
  ll_ls <- beta_loglik(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    steps = 1 - rev(steps),
    Hsgn = -1L
  )
  expect_equal(ll_gt, ll_ls)
  
  # equivalence of scores
  score_gt <- beta_score(
    test_param, yi = dat$yi, sei = dat$sei,
    steps = steps,
  )
  score_ls <- beta_score(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    steps = rev(1 - steps),
    Hsgn = -1L
  )
  expect_equal(score_gt, score_ls[c(1,2,4,3)], ignore_attr = TRUE)
  
  
  # equivalence of Hessians
  Hess_gt <- beta_hessian(
    test_param, yi = dat$yi, sei = dat$sei,
    steps = steps
  )
  Hess_ls <- beta_hessian(
    test_param_trans, yi = dat$yi, sei = dat$sei,
    steps = rev(1 - steps),
    Hsgn = -1L
  )
  expect_equal(Hess_gt, Hess_ls[c(1,2,4,3),c(1,2,4,3)], ignore_attr = TRUE)
  
  
  # equivalence of parameter estimates
  check_valence_equivalence(
    data = dat, 
    yi = yi, yi_neg = yi_neg, sei = sei,
    selection_type = "beta",
    steps = steps,
    priors = NULL
  )
  
  check_valence_equivalence(
    data = dat, 
    yi = yi, yi_neg = yi_neg, sei = sei,
    mean_mods = ~ design + year_c,
    selection_type = "beta",
    steps = steps,
    priors = NULL,
    vcov_type = "model-based"
  )
  
})

