#-------------------------------------------------------------------------------
# Functions for comparing scores to numerical derivatives

diffentiate_loglik <- function(
    param,
    f_ll, f_score, f_hess,
    steps,
    theta,
    from, to, N,
    yi,
    sei
) {
  
  param_seq <- seq(from, to, length.out = N + 2)
  
  f_ll_ <- function(x) {
    theta[param] <- x
    f_ll(theta = theta, yi = yi, sei = sei, steps = steps)
  }
  ll_seq <- sapply(param_seq, f_ll_)
  
  f_score_ <- function(x) {
    theta[param] <- x
    f_score(theta = theta, yi = yi, sei = sei, steps = steps)
  }
  score_seq <- sapply(param_seq, f_score_)

  f_hess_ <- function(x) {
    theta[param] <- x
    f_hess(theta = theta, yi = yi, sei = sei, steps = steps)[param,]
  }
  hess_seq <- sapply(param_seq, f_hess_)
  delta <- param_seq[3:(N + 2)] - param_seq[1:N]
  num_score <- (ll_seq[3:(N + 2)] - ll_seq[1:N]) / delta
  
  num_hess <- t(score_seq[,3:(N + 2)] - score_seq[,1:N]) / delta 
  
  data.frame(
    val = param_seq[2:(N + 1)],
    d_num = num_score,
    score = score_seq[param, 2:(N + 1)],
    h_num = num_hess,
    hess = t(hess_seq[,2:(N + 1)])
  )
}

check_all_derivatives <- function(
    steps = .025,
    selection_type = "step",
    estimator = "ML",
    ...,
    N = 100,
    crit = qnorm(0.975),
    params = "all"
) {
  
  cl <- match.call(expand.dots = TRUE)
  cl$N <- NULL
  cl$crit <- NULL
  cl[[1]] <- quote(selection_model)
  selmod_fit <- eval(cl, parent.frame())
  
  m <- match(c("data","yi", "sei"), names(cl), 0L)
  mf <- cl[c(1L, m)]
  mf[[1L]] <- quote(build_model_frame)
  mf <- eval(mf, parent.frame())
  
  # Evaluate yi, sei from model frame
  yi <- eval(cl$yi, envir = mf)
  sei <- eval(cl$sei, envir = mf)
  
  theta <- selmod_fit$est$Est
  
  if (identical(params,"all")) params <- 1:length(theta)
  
  from <- selmod_fit$est$Est - crit * selmod_fit$est$SE
  to <- selmod_fit$est$Est + crit * selmod_fit$est$SE
  
  if (selection_type == "step") {
    if (estimator == "ML") {
      f_ll <- step_loglik
      f_score <- step_score
      f_hess <- step_hessian
    } else if (estimator == "hybrid") {
      f_ll <- step_weighted_logpartlik
      f_score <- step_hybrid_score
      f_hess <- step_hybrid_jacobian
    }
  } else if (selection_type == "beta") {
    f_ll <- beta_loglik
    f_score <- beta_score
    f_hess <- beta_hessian
  }
  
  derivs <- mapply(
    diffentiate_loglik, 
    param = params, from = from[params], to = to[params],
    MoreArgs = list(
      f_ll = f_ll,
      f_score = f_score,
      f_hess = f_hess,
      steps = steps,
      theta = theta, 
      yi = yi,
      sei = sei,
      N = N
    ),
    SIMPLIFY = FALSE
  )
  
  names(derivs) <- selmod_fit$est$param[params]
  score_max_diff <- sapply(derivs, \(x) max(abs(x$d_num - x$score), na.rm = TRUE))
  score_range <- sapply(derivs, \(x) diff(range(x$score, na.rm = TRUE)))
  score_diff_over_range <- score_max_diff / pmax(1, score_range)
  p_dim <- length(theta)
  hess_max_diff <- sapply(derivs, \(x) apply(abs(x[,3 + 1:p_dim] - x[, 3 + p_dim + 1:p_dim]), 2, max, na.rm = TRUE))
  hess_range <- sapply(derivs, \(x) apply(x[,3 + p_dim + 1:p_dim], 2, \(x) diff(range(x, na.rm = TRUE)))) 
  hess_diff_over_range <- hess_max_diff / pmax(1, hess_range)
  deriv_dat <- do.call(rbind, derivs)
  
  deriv_dat$param <- rep(names(derivs), each = N)
  
  return(list(
    selmod_fit = selmod_fit, 
    score_diff_over_range = score_diff_over_range, 
    hess_diff_over_range = hess_diff_over_range,
    data = deriv_dat
  ))
  
}

#-------------------------------------------------------------------------------
# Functions for checking selection_model()


get_RE_params <- function(mod) {
  beta <- as.numeric(mod$beta)
  gamma <- log(mod$tau2)
  list(beta = beta, gamma = gamma)
}

get_selmodel_params <- function(sel_mod, theta = TRUE) {
  beta <- as.numeric(sel_mod$beta)
  gamma <- log(sel_mod$tau2)
  zeta <- if (sel_mod$type == "stepfun") log(sel_mod$delta[-1]) else log(sel_mod$delta)
  if (theta) {
    c(beta, gamma, zeta)
  } else {
    list(beta = beta, gamma = gamma, zeta = zeta)
  }
}

check_against_metafor_selmodel <- function(
    mod, type = "stepfun", steps = .025, 
    tol_LRT = 1e-8, tol_score = 5e-5, tol_param = 1e-4,
    ...,
    verbose = FALSE
) {

  suppressWarnings(
    sel_mod <- metafor::selmodel(
      mod, type = type, steps = steps, 
      control=list(optimizer = "nlminb", rel.tol = 1e-10, ...)
    )
  )

  theta <- get_selmodel_params(sel_mod)
  
  dat <- mod$data
  yi <- mod$call$yi
  sei <- mod$call$sei
  mods <- mod$call$mods
  
  if (type == "stepfun") {
    
    ll_full <- step_loglik(theta = theta, 
                           yi = sel_mod$yi, sei = sqrt(sel_mod$vi), 
                           X = sel_mod$X,
                           steps = steps)
    ll_null <- step_loglik(yi = sel_mod$yi, sei = sqrt(sel_mod$vi), 
                           X = sel_mod$X,
                           steps = steps,
                           beta = as.numeric(mod$beta), 
                           gamma = log(mod$tau2),
                           zeta = rep(0, length(steps)))
    score_full <- step_score(theta = theta, 
                             yi = sel_mod$yi, sei = sqrt(sel_mod$vi), 
                             X = sel_mod$X,
                             steps = steps)
    
    if (mod$method == "FE") {
      score_full <- score_full[-(length(mod$beta)+1)]
    } else {
      
      metafor_est <- c(
        sel_mod$beta,
        log(sel_mod$tau2),
        log(sel_mod$delta[-1])
      )
      
      my_est <- selection_model(
        data = dat,
        yi = yi,
        sei = sei,
        steps = steps,
        mean_mods = mods,
        selection_type = "step",
        estimator = "ML",
        theta = metafor_est
      )
      
    }
    
  } else if (type == "beta") {
    
    steps <- c(1e-05, 1 - 1e-05)
    ll_full <- beta_loglik(theta = theta, 
                           yi = sel_mod$yi, sei = sqrt(sel_mod$vi), 
                           X = sel_mod$X,
                           steps = steps)
    ll_null <- beta_loglik(yi = sel_mod$yi, sei = sqrt(sel_mod$vi), 
                           X = sel_mod$X,
                           steps = steps,
                           beta = as.numeric(mod$beta), 
                           gamma = log(mod$tau2),
                           zeta = rep(0, 2))
    score_full <- beta_score(theta = theta, 
                             yi = sel_mod$yi, sei = sqrt(sel_mod$vi), 
                             X = sel_mod$X,
                             steps = steps)
    
    if (mod$method == "FE") {
      score_full <- score_full[-(length(mod$beta)+1)]
    } else {
      
      metafor_est <- c(
        sel_mod$beta,
        log(sel_mod$tau2),
        log(sel_mod$delta)
      )
      
      my_est <- selection_model(
        data = dat,
        yi = yi,
        sei = sei,
        steps = steps,
        mean_mods = mods,
        selection_type = "beta",
        estimator = "ML",
        theta = metafor_est
      )
      
    }
  }

  # Likelihood ratio test statistic
  testthat::expect_equal(2 * (ll_full - ll_null), sel_mod$LRT, tolerance = tol_LRT)
  
  # Score vector
  if (max(abs(score_full)) > tol_score) print(score_full)
  testthat::expect_lt(max(abs(score_full)), tol_score)
  
  if (mod$method != "FE") {
    # parameter estimates
    est_compare <- data.frame(
      param = my_est$est$param,
      metafor = metafor_est,
      package = my_est$est$Est,
      SE = my_est$est$SE,
      diff = metafor_est - my_est$est$Est
    )
    
    if (verbose) print(est_compare)
    
    testthat::expect_equal(my_est$est$Est, metafor_est, tolerance = tol_param)
    
  }
}

check_dims <- function(mf, rows, cols) {
  testthat::expect_identical(nrow(mf), rows)
  testthat::expect_identical(ncol(mf), cols)
}


#-------------------------------------------------------------------------------
# Functions for checking r_meta()

check_target_m <- function(
    mean_smd, tau, omega, m, 
    cor_mu = 0.6, cor_sd = 0, 
    censor_fun = step_fun(), 
    n_ES_sim = n_ES_param(40, 1), 
    m_multiplier = 2
) {
  
  max_param_length <- max(lengths(list(mean_smd, tau, omega, m, cor_mu, cor_sd, censor_fun)))
  
  dat <- r_meta_categories(
    mean_smd = mean_smd, 
    tau = tau, 
    omega = omega, 
    m = m, 
    cor_mu = cor_mu, 
    cor_sd = cor_sd, 
    censor_fun = censor_fun, 
    n_ES_sim = n_ES_sim, 
    m_multiplier = m_multiplier
  )
  
  m_per_cat <- rep(m, length.out = length(table(dat$X)))
  studies_per_cat <- tapply(dat$studyid, dat$X, \(x) nlevels(droplevels(x)))
  
  # Number of levels of X equal to maximum number of parameter values
  testthat::expect_identical(max_param_length, nlevels(dat$X))

  # Number of studies per level of X equal to specified m's
  testthat::expect_identical(as.integer(studies_per_cat), m_per_cat)
  testthat::expect_identical(sum(m_per_cat), nlevels(dat$studyid))
  
  # Number of effect sizes per level of X gte specified m's
  testthat::expect_true(all(table(dat$X) >= m_per_cat))
  
}


check_model_structure <- function(
    mean_smd, tau, omega, m, 
    cor_mu = 0.6, cor_sd = 0, 
    censor_fun = step_fun(), 
    n_ES_sim = n_ES_param(40, 1), 
    m_multiplier = 1,
    check_varcomp = TRUE,
    verbose = FALSE
) {
  
  # Generate data with specified parameters
  
  dat <- r_meta_categories(
    mean_smd = mean_smd, 
    tau = tau, 
    omega = omega, 
    m = m, 
    cor_mu = cor_mu, 
    cor_sd = cor_sd, 
    censor_fun = censor_fun, 
    n_ES_sim = n_ES_sim, 
    m_multiplier = m_multiplier
  )

  
  # Determine moderator specification
  
  if (length(mean_smd) > 1L) {
    mods <- ~ 0 + X
  } else {
    mods <- ~ 1
  }
  
  
  # Determine random effects specification
  
  if (length(tau) > 1L) {
    RE_study <- ~ X | studyid
    struct_study <- "DIAG"
  } else {
    RE_study <- ~ 1 | studyid
    struct_study <- NULL
  }
  
  if (length(omega) > 1L) {
    RE_effect <- ~ X | esid
    struct_effect <- "DIAG"
  } else {
    RE_effect <- ~ 1 | esid
    struct_effect <- NULL
  } 
  
  random_spec <- list(RE_study, RE_effect)
  
  
  # Determine random effects structure and order of variance components in output
  
  if (is.null(struct_study) & is.null(struct_effect)) {
    struct <- "CS"
    var_params <- c(sigma.1 = tau, sigma.2 = omega)
  } else {
    struct <- c(struct_study, struct_effect)
    if (is.null(struct_study)) {
      var_params <- c(sigma = tau, tau. = omega)
    } else if (is.null(struct_effect)) {
      var_params <- c(sigma = omega, tau. = tau)
    } else {
      var_params <- c(tau. = tau, gamma. = omega)
    }
  }
  
  
  # Fit model in metafor
  
  if (cor_mu > 0) {
    Vmat <- metafor::vcalc(
      vi = dat$var_d, cluster = dat$studyid, obs = dat$esid,  
      rho = cor_mu, checkpd = FALSE
    )
    
    suppressWarnings(
      m_fit <- metafor::rma.mv(
        yi = dat$d, V = Vmat, data = dat,
        mods = mods,
        random = random_spec, struct = struct,
        sparse = TRUE
      )
    )
    
  } else {
    
    suppressWarnings(
      m_fit <- metafor::rma.mv(
        yi = dat$d, V = dat$var_d, data = dat,
        mods = mods,
        random = random_spec, struct = struct,
        sparse = TRUE
      )
    )
    
  }
  
  # Check CIs for fixed effects
  
  if (verbose) {
    print(m_fit)
    cat("mean_smd: ", mean_smd, "\n")
  }

  testthat::expect_true(all(mean_smd > m_fit$ci.lb))
  testthat::expect_true(all(mean_smd < m_fit$ci.ub))
  
  
  # Check CIs for variance components
  
  if (check_varcomp) {
    
    varcomp_CIs <- metafor::confint.rma.mv(m_fit, verbose = verbose)
    
    for (i in seq_along(var_params)) {
      var_name <- names(var_params)[i]
      if (verbose) {
        print(varcomp_CIs[[i]])
        cat(var_name, ":", var_params[i], "\n")
      }
      testthat::expect_gte(var_params[i], varcomp_CIs[[i]]$random[var_name,"ci.lb"])
      testthat::expect_lte(var_params[i], varcomp_CIs[[i]]$random[var_name,"ci.ub"])
    }
    
  }

}

check_step_score_hessian_bias <- function(
    mean_smd, 
    tau, 
    m, 
    cor_mu = 0.6, cor_sd = 0, 
    steps = .025, 
    weights = 1,
    mean_N = 60,
    mean_kj = 1,
    m_multiplier = 1,
    check_Hotelling = TRUE,
    hessian_threshold = 1,
    verbose = TRUE,
    seed = NULL,
    just_pvals = FALSE,
    score_type = "score"
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # create censoring functions given steps and weights
  if (is.list(weights)) {
    censor_fun <- lapply(weights, \(w) step_fun(cut_vals = steps, weights = w))
  } else {
    censor_fun <- step_fun(cut_vals = steps, weights = weights)
  }
  
  # create sample size sampler
  n_ES_sim <- n_ES_param(mean_N = mean_N, mean_ES = mean_kj)
  
  # Generate data with specified parameters
  
  dat <- r_meta_categories(
    mean_smd = mean_smd, 
    tau = tau, 
    omega = 0, 
    m = m, 
    cor_mu = cor_mu, 
    cor_sd = cor_sd, 
    censor_fun = censor_fun, 
    n_ES_sim = n_ES_sim, 
    m_multiplier = m_multiplier
  )
  dat$sig_bucket <- cut(dat$p_onesided, c(0, steps, 1))
  
  if (verbose) {
    print(table(dat$X, dat$sig_bucket))
  }
  
  # Determine parameter values
  beta <- mean_smd
  gamma <- log(tau^2)
  zeta <- if (is.list(weights)) as.vector(t(sapply(weights, log))) else log(weights)
  
  if (nlevels(dat$X) > 1L) {
    Xmat <- model.matrix(~ 0 + X, data = dat)
  }
  
  # Determine X specification
  X <- if (length(mean_smd) > 1L) Xmat else NULL
  
  # Determine U specification
  U <- if (length(tau) > 1L) Xmat else NULL

  # Determine Z specification
  if (is.list(weights)) {
    Z <- Xmat
    if (length(steps) > 1L) Z <- rep(list(Z), length(steps))
  } else {
    Z <- NULL
  }
  
 
  # Calculate score contributions
  
  if (score_type == "score") {
    
    score_i <- step_score(
      yi = dat$d,
      sei = dat$sd_d,
      beta = beta,
      gamma = gamma,
      zeta = zeta,
      steps = steps,
      X = X,
      U = U,
      Z = Z,
      contributions = TRUE
    )
  } else if (score_type == "hybrid") {
    
    score_i <- step_hybrid_score(
      yi = dat$d,
      sei = dat$sd_d,
      beta = beta,
      gamma = gamma,
      zeta = zeta,
      steps = steps,
      X = X,
      U = U,
      Z = Z,
      contributions = TRUE
    )    
    
    
  }

  one_var_ps <- apply(score_i, 2, \(x) t.test(x)$p.value)
  Hotelling_T2 <- DescTools::HotellingsT2Test(score_i, mu = rep(0, ncol(score_i)))
  
  if (just_pvals) {
    pval_dat <- data.frame(
      component = c(names(one_var_ps), "Hotelling"), 
      pval = c(one_var_ps, Hotelling_T2$p.value)
    )
    return(pval_dat)
  }
  
  if (verbose) {
    summary_stats <- data.frame(
      M = colMeans(score_i), 
      SD = apply(score_i, 2, sd),
      pval = one_var_ps
    )
    print(summary_stats)
  }

  testthat::expect_true(all(one_var_ps > .1))
  
  if (check_Hotelling) {
    if (verbose) {
      print(Hotelling_T2)
    }
    
    testthat::expect_gt(Hotelling_T2$p.value, .1)
  }

  if (!is.null(hessian_threshold)) {
    hess_mat <- step_hessian(
      yi = dat$d,
      sei = dat$sd_d,
      beta = beta,
      gamma = gamma,
      zeta = zeta,
      steps = steps,
      X = X,
      U = U,
      Z = Z
    )
    
    var_score <- crossprod(score_i) / m
    Bartlett_diff <- var_score + hess_mat / m
    trace_stat <- sum(diag(var_score %*% solve(-hess_mat / m))) - nrow(hess_mat)
    if (verbose) {
      cat("-Hessian: \n")
      print(-hess_mat / m)
      cat("Var(score_i):\n")
      print(var_score)
      cat("Difference:\n")
      print(Bartlett_diff)
      cat("Maximum difference:", max(abs(Bartlett_diff)), "\n")
      cat("Trace statistic:", trace_stat, "\n")
    }
    
    testthat::expect_true(all(abs(Bartlett_diff) < hessian_threshold))
  }

}


check_beta_score_hessian_bias <- function(
    mean_smd, 
    tau, 
    m, 
    cor_mu = 0.6, cor_sd = 0, 
    steps = c(.025, .975),
    lambdas = c(1, 1),
    mean_N = 60,
    mean_kj = 1,
    m_multiplier = 1,
    check_Hotelling = TRUE,
    hessian_threshold = 1,
    verbose = TRUE,
    seed = NULL,
    just_pvals = FALSE,
    fit_mod = FALSE
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  censor_fun <- beta_fun(delta_1 = lambdas[1], delta_2 = lambdas[2],
                             trunc_1 = steps[1], trunc_2 = steps[2])
  
  # create sample size sampler
  n_ES_sim <- n_ES_param(mean_N = mean_N, mean_ES = mean_kj)
  
  # Generate data with specified parameters
  
  dat <- r_meta_categories(
    mean_smd = mean_smd, 
    tau = tau, 
    omega = 0, 
    m = m, 
    cor_mu = cor_mu, 
    cor_sd = cor_sd, 
    censor_fun = censor_fun, 
    n_ES_sim = n_ES_sim, 
    m_multiplier = m_multiplier
  )
  
  dat$sig_bucket <- cut(dat$p_onesided, c(0, steps, 1))
  
  if (verbose) {
    print(table(dat$X, dat$sig_bucket))
  }
  
  
  # Determine parameter values
  beta <- mean_smd
  gamma <- log(tau^2)
  zeta <- log(lambdas)
  
  if (nlevels(dat$X) > 1L) {
    Xmat <- model.matrix(~ 0 + X, data = dat)
  }
  
  # Determine X specification
  if (length(mean_smd) > 1L) {
    X <- Xmat
    X_frm <- ~ X
  } else {
    X <- NULL
    X_frm <- NULL
  }
  
  # Determine U specification
  if (length(tau) > 1L) {
    U <- Xmat
    U_frm <- ~ X
  } else {
    U <- NULL
    U_frm <- NULL
  }

  # Fit model to the data
  if (fit_mod) {
    selmod_fit <- selection_model(
      data = dat, 
      yi = quote(d), sei = quote(sd_d), 
      mean_mods = X_frm,
      var_mods = U_frm,
      steps = steps, 
      selection_type = "beta"
    )
    selmod_fit$est$param_val <- c(beta, gamma, zeta)
    print(selmod_fit$est)
  }
  
  # Calculate score contributions
  
  score_i <- beta_score(
    yi = dat$d,
    sei = dat$sd_d,
    beta = beta,
    gamma = gamma,
    zeta = zeta,
    steps = steps,
    X = X,
    U = U,
    contributions = TRUE
  )
  
  one_var_ps <- apply(score_i, 2, \(x) t.test(x)$p.value)
  Hotelling_T2 <- DescTools::HotellingsT2Test(score_i, mu = rep(0, ncol(score_i)))

  if (just_pvals) {
    pval_dat <- data.frame(
      component = c(names(one_var_ps), "Hotelling"), 
      pval = c(one_var_ps, Hotelling_T2$p.value)
    )
    return(pval_dat)
  }
  
  if (verbose) {
    summary_stats <- data.frame(
      M = colMeans(score_i), 
      SD = apply(score_i, 2, sd),
      pval = one_var_ps
    )
    print(summary_stats)
  }
  
  testthat::expect_true(all(one_var_ps > .1))
  
  if (check_Hotelling) {
    if (verbose) {
      print(Hotelling_T2)
    }
    
    testthat::expect_gt(Hotelling_T2$p.value, .1)
  }

  if (!is.null(hessian_threshold)) {
    hess_mat <- beta_hessian(
      yi = dat$d,
      sei = dat$sd_d,
      beta = beta,
      gamma = gamma,
      zeta = zeta,
      steps = steps,
      X = X,
      U = U
    )
    
    var_score <- crossprod(score_i) / m
    Bartlett_diff <- var_score + hess_mat / m
    trace_stat <- sum(diag(var_score %*% solve(-hess_mat / m))) - nrow(hess_mat)
    if (verbose) {
      cat("-Hessian: \n")
      print(-hess_mat / m)
      cat("Var(score_i):\n")
      print(var_score)
      cat("Difference:\n")
      print(Bartlett_diff)
      cat("Maximum difference:", max(abs(Bartlett_diff)), "\n")
      cat("Trace statistic:", trace_stat, "\n")
    }
    
    testthat::expect_true(all(abs(Bartlett_diff) < hessian_threshold))
  }
}


eval_beta_range <- function(delta, steps, eps = 1e-3, plot = FALSE) {
  p <- seq(0,1,eps)
  b <- beta_fun(delta_1 = delta[1], delta_2 = delta[2], 
                    trunc_1 = steps[1], trunc_2 = steps[2])
  wts <- b(p)
  if (plot) plot(p, wts)
  max_val <- (delta[1] - 1) / (sum(delta) - 2)
  
  list(range = range(wts), max_val = max_val)
}


check_profiling_equivalence <- function(
    yi, sei, 
    steps, 
    ai = NULL, 
    cluster = NULL, 
    X = NULL, 
    U = NULL, 
    Z0 = NULL, 
    Z = NULL, 
    subset = NULL,
    vcov_type = "robust",
    optimizer_control = list(),
    tol = 1e-8,
    score_tol = 1e-8,
    jac_tol = 1e-8,
    verbose = FALSE
) {
  
  
  hybrid_fit <- fit_selection_model(
    yi = yi, sei = sei, pi = pnorm(yi / sei, lower.tail = FALSE),
    steps = steps, 
    ai = ai, 
    cluster = cluster, 
    X = X, 
    U = U, 
    Z0 = Z0, 
    Z = Z, 
    subset = subset,
    vcov_type = vcov_type,
    selection_type = "step",
    estimator = "hybrid-full",
    optimizer_control = list(),
  )

  theta <- hybrid_fit$est
  x_index <- if (is.null(X)) 1L else 1:ncol(X)
  
  parse_params_full <- 
    parse_step_params(
      theta = theta,
      yi = yi,
      sei = sei,
      ai = ai,
      steps = steps,
      X = X,
      U = U,
      Z0 = Z0,
      Z = Z,
      calc_Ai = TRUE
    )
  
  scores_full <- step_hybrid_score(
    theta = theta,
    yi = yi,
    sei = sei,
    ai = ai,
    steps = steps,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z
  )

  jac_full <- step_hybrid_jacobian(
    theta = theta,
    yi = yi,
    sei = sei,
    ai = ai,
    steps = steps,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z
  )
  
  theta_na <- theta
  theta_na[x_index] <- NA_real_
  
  profile_fit <- fit_selection_model(
    yi = yi, sei = sei, pi = pnorm(yi / sei, lower.tail = FALSE),
    steps = steps, 
    ai = ai, 
    cluster = cluster, 
    X = X, 
    U = U, 
    Z0 = Z0, 
    Z = Z, 
    subset = subset,
    vcov_type = vcov_type,
    selection_type = "step",
    estimator = "hybrid",
    optimizer_control = list(),
  )
  
  parse_params_prof <- 
    parse_step_params(
      theta = theta_na,
      yi = yi,
      sei = sei,
      ai = ai,
      steps = steps,
      X = X,
      U = U,
      Z0 = Z0,
      Z = Z,
      calc_Ai = TRUE
    )
  
  scores_prof <- step_hybrid_profile_score(
    theta = theta_na[-x_index],
    yi = yi,
    sei = sei,
    ai = ai,
    steps = steps,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z
  )
  
  jac_prof <- step_hybrid_profile_jacobian(
    theta = theta_na[-x_index],
    yi = yi,
    sei = sei,
    ai = ai,
    steps = steps,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z
  )

  if (verbose) {
    cat("hybrid estimates:", hybrid_fit$est, "\n")
    cat("profiled estimates:", profile_fit$est, "\n")
  }
  
  
  if (length(tol) == 1) tol <- rep(tol, 3L)
  testthat::expect_equal(parse_params_full$mu, parse_params_prof$mu, ignore_attr = FALSE, tolerance = tol[1])
  testthat::expect_equal(parse_params_full$eta, parse_params_prof$eta, ignore_attr = FALSE, tolerance = tol[2])
  testthat::expect_equal(parse_params_full$lambda_full, parse_params_prof$lambda_full, ignore_attr = FALSE, tolerance = tol[3])

  if (verbose) {
    cat("full scores:", scores_full, "\n")
    cat("profiled scores:", scores_prof, "\n")
  }
  
  if (!is.infinite(score_tol)) {
    testthat::expect_lt(max(abs(scores_full)), score_tol)
    testthat::expect_lt(max(abs(scores_prof)), score_tol)
  }

  if (verbose) {
    cat("full jacobian:\n", jac_full, "\n")
    cat("profiled jacobian:\n", jac_prof, "\n")
  }  
  
  testthat::expect_equal(jac_full[-x_index, -x_index],jac_prof, tolerance = jac_tol)
}


#-------------------------------------------------------------------------------
# Functions for checking print() and summary()

print_and_parse <- function(mod,...) {
  p <- testthat::capture_output_lines(print(mod,...), width = 1000)
  d <- do.call(rbind, strsplit(p, " +"))
  d[,-1]
}

pull_argument <- function(s, str) {
  x <- s[grepl(str, s)]
  trimws(substr(x, nchar(str) + 2L, nchar(x)))
}

check_selmodel_summary <- function(mod, ...) {
  s <- testthat::capture_output_lines(summary(mod, ...), width = 1000)
  p <- testthat::capture_output_lines(print(mod, ...), width = 1000)
  mod_str <- trimws(s[1])
  steps <- pull_argument(s, "Steps:") |> strsplit(",") |> unlist() |> as.numeric()
  estimator <- pull_argument(s,"Estimator:") 
  boot_type <- pull_argument(s, "Bootstrap type:")
  R <- pull_argument(s, "Number of bootstrap replications:") |> as.integer()
  
  headers <- which(grepl("estimates:", s))
  mean_effects <- s[(headers[1] + 3L):(headers[2] - 2)]
  var_effects <- s[(headers[2] + 3L):(headers[3] - 2)]
  sel_effects <- s[(headers[3] + 3L):length(s)]
  sum_effects <- strsplit(c(mean_effects, var_effects, sel_effects), " +")
  sum_effects <- sapply(sum_effects, \(x) x[2])
  p_effects <- sapply(strsplit(p[-1], " +"), \(x) x[2])

  mod_type <- if (mod$selection_type == "step") "Step Function Model" else "Beta Density Model"
  if (inherits(mod, "boot.selmodel")) {
    mod_type <- paste(mod_type, "with Cluster Bootstrapping")
    
    testthat::expect_identical(boot_type, unique(mod$bootstrap_type))
    testthat::expect_identical(R, unique(mod$est$bootstraps))
  } 
  
  testthat::expect_identical(mod_str, mod_type)
  testthat::expect_identical(steps, mod$steps)
  
  estimator_type <- if (mod$estimator == "ML") "maximum likelihood" else "hybrid estimating equations"
  testthat::expect_identical(estimator, estimator_type)
 
  testthat::expect_identical(sum_effects, p_effects) 
}
