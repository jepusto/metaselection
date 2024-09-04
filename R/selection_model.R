#' 
#' @importFrom Formula as.Formula


build_model_frame <- function(
    data,
    yi,
    sei,
    pi,
    ai,
    cluster,
    subset = NULL,
    mean_mods = NULL,
    var_mods = NULL,
    sel_mods = NULL,
    sel_zero_mods = NULL
) {
  if (missing(yi) || missing(sei)) stop("You must specify an effect size variable yi and a standard error variable sei.")
  yi_str <- deparse(substitute(yi))
  sei_str <- deparse(substitute(sei))
  pi_str <- if (missing(pi)) NULL else deparse(substitute(pi))
  ai_str <- if (missing(ai)) NULL else deparse(substitute(ai))
  cl_str <- if (missing(cluster)) NULL else deparse(substitute(cluster))
  yi_sei <- reformulate(c(sei_str, pi_str, ai_str, cl_str), response = yi_str)
  subset <- eval(substitute(subset), envir = data)
  
  formula_list <- list(yi_sei = yi_sei)
  if (!is.null(mean_mods)) formula_list$mean_mods <- mean_mods
  if (!is.null(var_mods)) formula_list$var_mods <- var_mods
  if (!is.null(sel_zero_mods)) formula_list$sel_zero_mods <- sel_zero_mods
  if (!is.null(sel_mods)) {
    if (is.list(sel_mods)) formula_list <- c(formula_list, sel_mods) else formula_list$sel_modss <- sel_mods
  } 

  names(formula_list) <- NULL
  combined_formula <- do.call(Formula::as.Formula, formula_list)
  model.frame(combined_formula, data = data, subset = subset)
}

find_starting_values <- function(
  yi, sei, 
  selection_type,
  steps, 
  X = NULL, 
  U = NULL, 
  Z0 = NULL, 
  Z = NULL
) {
  # Use weighted lm() for betas
  
  if (is.null(X)) {
    beta_start <- weighted.mean(x = yi, w = 1 / sei^2)
    ri_sq <- (yi - beta_start)^2
  } else {
    beta_fit <- stats::lm.wfit(x = X, y = yi, w = 1 / sei^2)
    beta_start <- unname(beta_fit$coefficients)
    ri_sq <- stats::residuals(beta_fit)^2
  }
  
  # Compute squared residuals
  if (is.null(U)) {
    gamma_start <- log(mean(ri_sq))
  } else {
    gamma_fit <- stats::glm.fit(x = U, y = ri_sq, family = stats::quasipoisson())
    gamma_start <- unname(gamma_fit$coefficients)
  }
  
  # Set selection parameters to zero throughout
  if (selection_type == "step") {
    zeta0_start <- if (is.null(Z0)) NULL else rep(0, ncol(Z0))
    zeta_start <- if (is.null(Z)) rep(0, length(steps)) else lapply(Z, \(x) rep(0, ncol(x)))
  } else {
    zeta0_start <- NULL
    zeta_start <- rep(0, 2L)
  }
  
  theta <- unlist(c(beta = beta_start, gamma = gamma_start, zeta0 = zeta0_start, zeta = zeta_start), use.names = FALSE)
  
  return(theta)
  
}

fit_selection_model <- function(
  yi, sei, pi, steps, 
  ai = NULL, 
  cluster = NULL, 
  X = NULL, 
  U = NULL, 
  Z0 = NULL, 
  Z = NULL, 
  subset = NULL,
  make_sandwich = TRUE,
  selection_type = "step",
  estimator = "ML",
  theta = NULL,
  optimizer = if (estimator == "ML") "BFGS" else "nleqslv",
  optimizer_control = list(),
  use_jac = TRUE
) {

  if (!is.null(subset)) {
    yi <- yi[subset]
    sei <- sei[subset]
    pi <- pi[subset]
    if (!is.null(ai)) ai <- ai[subset]
    if (!is.null(cluster)) cluster <- cluster[subset]
    if (!is.null(X)) X <- X[subset,,drop=FALSE]
    if (!is.null(U)) U <- U[subset,,drop=FALSE]
    if (!is.null(Z0)) Z0 <- Z0[subset,,drop=FALSE]
    if (!is.null(Z)) Z <- Z[subset,,drop=FALSE]
  }
  
  if (is.null(theta)) {
    # Compute starting values for parameters if not provided
    theta <- find_starting_values(
      yi = yi, sei = sei, 
      selection_type = selection_type, steps = steps,
      X = X, U = U, Z0 = Z0, Z = Z
    )
  }
  
  params <- parse_step_params(
    theta = theta,
    yi = yi, sei = sei,
    pi = pi, ai = ai,
    steps = steps,
    X = X, U = U, Z0 = Z0, Z = Z, 
    calc_Ai = FALSE
  )
  
  if (estimator == "ML") {
    
    # Optimize the log likelihood using optimx()
    
    optimizer_control$maximize <- TRUE
    
    utils::capture.output(
      if (selection_type == "step") {
        
        mle_est <- optimx::optimx(
          par = theta, 
          fn = step_loglik, 
          gr = step_score,
          yi = yi, sei = sei, pi = pi, ai = ai,
          steps = steps,
          X = X, U = U, Z0 = Z0, Z = Z,
          method = optimizer,
          control = optimizer_control
        )
        
      } else if (selection_type == "beta") {
        
        mle_est <- optimx::optimx(
          par = theta, 
          fn = beta_loglik, 
          gr = beta_score,
          yi = yi, sei = sei, pi = pi, ai = ai,
          steps = steps,
          X = X, U = U, 
          method = optimizer,
          control = optimizer_control
        )
        
      }
    )
    
    if (make_sandwich == "raw") {
      return(mle_est)
    } 
    
    max_method <- row.names(mle_est)[which.max(mle_est$value)]
    theta_names <- 1:length(theta)
    theta <- as.numeric(mle_est[max_method, theta_names])
    info <- mle_est[max_method, -theta_names]
    names(theta) <- params$H_names
    
    if (isFALSE(make_sandwich)) {
      return(theta)
    }
    
  } else if (estimator == "hybrid-full") {
    
    optimizer <- match.arg(optimizer, choices = c("nleqslv","rootSolve"), several.ok = TRUE)
    
    info <- list()

    jac <- if (use_jac) step_hybrid_jacobian else NULL
    
    if ("nleqslv" %in% optimizer) {
      
      nleqslv_args <- list(
        x = theta, 
        fn = step_hybrid_score,
        jac = jac,
        yi = yi, sei = sei, pi = pi, ai = ai,
        steps = steps,
        X = X, U = U, Z0 = Z0, Z = Z,
        jacobian = TRUE
      )
      if ("method" %in% names(optimizer_control)) {
        nleqslv_args$method <- optimizer_control$method
        optimizer_control$method <- NULL
      }
      if ("global" %in% names(optimizer_control)) {
        nleqslv_args$global <- optimizer_control$global
        optimizer_control$global <- NULL
      }
      if ("xscalm" %in% names(optimizer_control)) {
        nleqslv_args$xscalm <- optimizer_control$xscalm
        optimizer_control$xscalm <- NULL
      }
      nleqslv_args$control <- optimizer_control
      
      nleqslv_res <- do.call(nleqslv::nleqslv, args = nleqslv_args)
      
      nleqslv_norm <- sum(nleqslv_res$fvec^2)
      info$nleqslv <- nleqslv_res[-1]  
      info$nleqslv$f_norm <- nleqslv_norm
      
    }
    
    if ("rootSolve" %in% optimizer) {
      rootsolve <- list(
        f = step_hybrid_score,
        start = theta, 
        jacfun = jac,
        yi = yi, sei = sei, pi = pi, ai = ai,
        steps = steps,
        X = X, U = U, Z0 = Z0, Z = Z
      )
      rootsolve_args <- c(rootsolve, optimizer_control)
      
      rootsolve_res <- do.call(rootSolve::multiroot, args = rootsolve_args)
      
      rootsolve_norm <- sum(rootsolve_res$f.root^2)
      info$rootsolve <- rootsolve_res[-1]
      info$rootsolve$f_norm <- rootsolve_norm
      
    }

    if (identical(optimizer, "nleqslv") || (all(c("nleqslv","rootSolve") %in% optimizer) && nleqslv_norm < rootsolve_norm)) {
      max_method <- if ("method" %in% names(nleqslv_args)) nleqslv_args$method else "default"
      max_method <- paste("nleqslv", max_method)
      theta <- nleqslv_res$x
    } else  {
      max_method <- "rootSolve"
      theta <- rootsolve_res$root
    }
    
    names(theta) <- params$H_names
    
    if (make_sandwich == "raw") {
      return(list(est = theta, max_method = max_method, info = info))
    }
    
    if (isFALSE(make_sandwich)) {
      return(theta)
    }
    
  }  else if (estimator == "hybrid") {
    
    optimizer <- match.arg(optimizer, choices = c("nleqslv","rootSolve"), several.ok = TRUE)
    
    info <- list()
    
    x_index <- if (is.null(X)) 1L else 1:ncol(X)
    theta <- theta[-x_index]
     
    jac <- if (use_jac) step_hybrid_profile_jacobian else NULL
    
    if ("nleqslv" %in% optimizer) {
      
      nleqslv_args <- list(
        x = theta, 
        fn = step_hybrid_profile_score,
        jac = jac,
        yi = yi, sei = sei, pi = pi, ai = ai,
        steps = steps,
        X = X, U = U, Z0 = Z0, Z = Z,
        jacobian = TRUE
      )
      if ("method" %in% names(optimizer_control)) {
        nleqslv_args$method <- optimizer_control$method
        optimizer_control$method <- NULL
      }
      if ("global" %in% names(optimizer_control)) {
        nleqslv_args$global <- optimizer_control$global
        optimizer_control$global <- NULL
      }
      if ("xscalm" %in% names(optimizer_control)) {
        nleqslv_args$xscalm <- optimizer_control$xscalm
        optimizer_control$xscalm <- NULL
      }
      nleqslv_args$control <- optimizer_control
      
      nleqslv_res <- do.call(nleqslv::nleqslv, args = nleqslv_args)
      
      nleqslv_norm <- sum(nleqslv_res$fvec^2)
      info$nleqslv <- nleqslv_res[-1]  
      info$nleqslv$f_norm <- nleqslv_norm
      
    }
    
    if ("rootSolve" %in% optimizer) {
      rootsolve <- list(
        f = step_hybrid_profile_score,
        start = theta, 
        jacfun = jac,
        yi = yi, sei = sei, pi = pi, ai = ai,
        steps = steps,
        X = X, U = U, Z0 = Z0, Z = Z
      )
      rootsolve_args <- c(rootsolve, optimizer_control)
      
      rootsolve_res <- do.call(rootSolve::multiroot, args = rootsolve_args)
      
      rootsolve_norm <- sum(rootsolve_res$f.root^2)
      info$rootsolve <- rootsolve_res[-1]
      info$rootsolve$f_norm <- rootsolve_norm
      
    }
    
    if (identical(optimizer, "nleqslv") || (all(c("nleqslv","rootSolve") %in% optimizer) && nleqslv_norm < rootsolve_norm)) {
      max_method <- if ("method" %in% names(nleqslv_args)) nleqslv_args$method else "default"
      max_method <- paste("nleqslv", max_method)
      theta <- nleqslv_res$x
    } else  {
      max_method <- "rootSolve"
      theta <- rootsolve_res$root
    }
    
    theta <- c(rep(NA_real_, max(x_index)), theta)
    
    params <- parse_step_params(
      theta = theta,
      yi = yi, sei = sei,
      pi = pi, ai = ai,
      steps = steps,
      X = X, U = U, Z0 = Z0, Z = Z, 
      calc_Ai = FALSE
    )
    
    theta[x_index] <- params$beta

    names(theta) <- params$H_names
    
    if (make_sandwich == "raw") {
      return(list(est = theta, max_method = max_method, info = info))
    }

    if (isFALSE(make_sandwich)) {
      return(theta)
    }
  }

  
  # Compute score contributions using MLEs
  
  if (estimator == "ML") {
    
    if (selection_type == "step") {
      
      scores <- step_score(theta = theta, 
                           yi = yi, sei = sei, pi = pi, ai = ai,
                           steps = steps,
                           X = X, U = U, Z0 = Z0, Z = Z, contributions = TRUE)
      
    } else if (selection_type == "beta") {
      
      scores <- beta_score(theta = theta, 
                           yi = yi, sei = sei, pi = pi, ai = ai,
                           steps = steps,
                           X = X, U = U, contributions = TRUE)
      
    }
    
  
  } else if (estimator %in% c("hybrid-full","hybrid")) {
    
    scores <- step_hybrid_score(theta = theta, 
                                yi = yi, sei = sei, pi = pi, ai = ai,
                                steps = steps,
                                X = X, U = U, Z0 = Z0, Z = Z, contributions = TRUE)
    
    
  }

  
  # Compute Hessian using MLEs
  
  if (estimator == "ML") {
    if (selection_type == "step") {
      
      hess <- step_hessian(theta = theta, 
                           yi = yi, sei = sei, pi = pi, ai = ai,
                           steps = steps,
                           X = X, U = U, Z0 = Z0, Z = Z)
      
    } else if (selection_type == "beta") {
      
      hess <- beta_hessian(theta = theta, 
                           yi = yi, sei = sei, pi = pi, ai = ai,
                           steps = steps,
                           X = X, U = U)
      
    }
    
  
  } else if (estimator %in% c("hybrid-full","hybrid")) {
    
    hess <- step_hybrid_jacobian(theta = theta, 
                                 yi = yi, sei = sei, pi = pi, ai = ai,
                                 steps = steps,
                                 X = X, U = U, Z0 = Z0, Z = Z)
  }
  
  hess_inv <- MASS::ginv(hess)
  
  # Compute sandwich variance estimator using score and Hessian (Eq. 52-54)
  
  if (is.null(cluster)) {
    meat <- crossprod(scores)
  } else {
    score_j <- rowsum(scores, group = cluster)
    meat <- crossprod(score_j)
  }

  sandwich <- hess_inv %*% meat %*% t(hess_inv)
  
  rownames(hess) <- names(theta)
  colnames(hess) <- names(theta)
  
  list(
    est = theta, 
    vcov = sandwich, 
    method = max_method,
    info = info
  )
  
}

bootstrap_selmodel <- function(
    yi, 
    sei, 
    pi, 
    steps, 
    ai = NULL,
    cluster = NULL, 
    X = NULL, 
    U = NULL, 
    Z0 = NULL, 
    Z = NULL, 
    make_sandwich = TRUE,
    selection_type = "step",
    estimator = "ML",
    theta = NULL,
    optimizer = if (estimator == "ML") "BFGS" else "nleqslv",
    optimizer_control = list(),
    use_jac = TRUE,
    wtype = c("multinomial", "exponential")
) {
  
  if (!is.null(cluster)) {
    cluster <- factor(cluster)
    m <- nlevels(cluster)
    cluster_numeric <- as.integer(cluster)
  } else {
    m <- length(yi)
    cluster_numeric <- 1:m
  }
  
  if (wtype == "multinomial") {
    
    i <- sample.int(n = m, replace = TRUE)
    dim(i) <- c(1, m)
    cluster_w <- apply(i, 1, tabulate, m)
    
  } else if (wtype == "exponential") {
    
    cluster_w <- stats::rexp(n = m)
    
  }
  
  if (is.null(ai)) {
    wi <- cluster_w[cluster_numeric]
  } else {
    wi <- ai * cluster_w[cluster_numeric]
  }
  
  non_zero_cl <- cluster_w > 0
  cl_subset <- if (all(non_zero_cl)) NULL else non_zero_cl[cluster_numeric]
  res <- fit_selection_model(
    yi = yi, 
    sei = sei, 
    pi = pi, 
    steps = steps,
    ai = wi, 
    cluster = cluster, 
    subset = cl_subset,
    X = X, U = U, Z0 = Z0, Z = Z,
    make_sandwich = make_sandwich, 
    selection_type = selection_type,
    estimator = estimator, 
    theta = theta,
    optimizer = optimizer, 
    optimizer_control = optimizer_control,
    use_jac = use_jac
  )
  
  if (make_sandwich) {
    est <- data.frame(
      param = names(res$est),
      Est = as.numeric(res$est),
      SE = sqrt(diag(res$vcov))
    )
  } else {
    est <- data.frame(
      param = names(res),
      Est = as.numeric(res)
    )
  }
  
  return(est)
}

#' @title Estimate step or beta selection model
#' 
#' @description Estimate step or beta selection model
#' 
#' 
#' @param data a data frame or tibble containing the meta-analytic data
#' @param yi vector of effect sizes
#' @param sei vector of sampling standard error
#' @param pi vector of one-sided p-values
#' @param ai vector of analytic weight
#' @param cluster vector indicating cluster membership 
#' @param selection_type character string specifying the type selection model to estimate, with possible options \code{"step"} or \code{"beta"}.
#' @param steps numeric vector of one or more values specifying the thresholds of p-values where the likelihood of effects being selected changes 
#' @param mean_mods optional model formula indicating moderators related to the average effect size magnitude
#' @param var_mods optional model formula indicating moderators related to effect size heterogeneity
#' @param sel_mods optional model formula indicating moderators related to the probability of selection 
#' @param sel_zero_mods optional model formula indicating moderators related to the probability of selection 
#' @param subset optional logical expression indicating the subset of the data to keep in the analysis
#' @param make_sandwich logical with \code{TRUE} (the default) indicating to calculate standard errors using sandwich estimators
#' @param conf_level desired coverage level for confidence intervals, with the default value set to \code{.95}
#' @param estimator vector indicating whether to use the maximum likelihood or the hybrid estimator, with possible options \code{"ML"}, \code{"hybrid"}, and \code{"hybrid-full"}. If running the beta selection model, only the maximum likelihood estimator, \code{"ML"}, is available. For step function models, both maximum likelihood and hybrid estimators are available.  
#' @param theta optional numeric vector of containing starting values of the parameters to be estimated 
#' @param optimizer character string indicating the optimizer to use
#' @param optimizer_control an optional list of control parameters to be used for optimization
#' @param use_jac logical with \code{TRUE} (the default) indicating that the Jacobian ...???
#' @param bootstrap character string specifying the type of bootstrap to run, with possible options \code{"none"} (the default), \code{"exponential"}, \code{"multinomial"}.
#' @param boot_CI character string specifying the type of bootstrap confidence interval to calculate, with possible options \code{"percentile"} (the default), \code{"student"}, \code{"large-sample"}, or \code{"none"}. 
#' @param R number of bootstrap replications, with the default set to \code{1999}
#' @param ... further arguments passed to \code{simhelpers::bootstrap_CIs}.
#' 
#' @returns A numeric vector.
#' 
#' @export
#' 
#' 
#' @examples
#' 
#' 
#' selection_model(data = dat,
#'                 yi = d,
#'                 sei = sd_d,
#'                 pi = p_onesided,
#'                 cluster = studyid,
#'                 steps = 0.025,
#'                 estimator = "ML",
#'                 bootstrap = "multinomial",
#'                 boot_CI = "percentile",
#'                 R = 49)
#' 


selection_model <- function(
    data,
    yi,
    sei,
    pi,
    ai,
    cluster,
    selection_type = c("step","beta"),
    steps = if (selection_type == "step") .025 else c(.025, .975),
    mean_mods = NULL,
    var_mods = NULL,
    sel_mods = NULL,
    sel_zero_mods = NULL,
    subset = NULL,
    make_sandwich = TRUE,
    conf_level = .95,
    estimator = c("ML","hybrid","hybrid-full"),
    theta = NULL,
    optimizer = NULL,
    optimizer_control = list(),
    use_jac = TRUE,
    bootstrap = "none",
    boot_CI = "percentile",
    R = 1999,
    ...
) {
  
  selection_type <- match.arg(selection_type)
  estimator <- match.arg(estimator)
  bootstrap <- match.arg(bootstrap, c("none","exponential","multinomial"))
  
  if (identical(R, 0)) {
    bootstrap <- "none"
    boot_CI <- "large-sample"
  } else if (bootstrap == "none") {
    boot_CI <- "large-sample"
  } else {
    boot_CI <- match.arg(boot_CI, c("large-sample","percentile","student","basic", "none"), several.ok = TRUE)  
  }
  
  if (is.null(optimizer)) {
    if (selection_type == "step") {
      optimizer <- if (estimator == "ML") "Rvmmin" else "nleqslv"
    } else {
      optimizer <- "nlminb"
    }
  }
  
  if (missing(steps)) stop("Must supply steps argument.")
  if (!inherits(steps, "numeric") || min(steps) <= 0 || max(steps) >= 1) stop("steps must be a numeric vector with all entries in the interval (0,1).")
  if (selection_type == "beta") {
    if (length(steps) != 2L) stop("steps must be a numeric vector of length 2 when selection_type = 'beta'.")
    if (!is.null(sel_mods)) stop("sel_mods must be NULL when selection_type = 'beta'.")
    if (!is.null(sel_zero_mods)) stop("sel_zero_mods must be NULL when selection_type = 'beta'.")
    if (estimator != "ML") stop("estimator must be equal to 'ML' when selection_type = 'beta'.")
  }
  
  # Create common model frame
  
  cl <- match.call()
  m <- match(c("data","yi", "sei", "pi", "ai", "cluster","subset", "mean_mods", 
               "var_mods", "sel_mods", "sel_zero_mods"), names(cl), 0L)
  mf <- cl[c(1L, m)]
  mf[[1L]] <- quote(build_model_frame)
  mf <- eval(mf, parent.frame())
  
  # Evaluate yi, sei, pi, ai from model frame
  
  yi <- eval(cl$yi, envir = mf)
  sei <- eval(cl$sei, envir = mf)
  pi <- if (missing(pi)) pnorm(yi / sei, lower.tail = FALSE) else eval(cl$pi, envir = mf)
  ai <- if (missing(ai)) NULL else eval(cl$ai, envir = mf)
  cluster <- if (missing(cluster)) NULL else eval(cl$cluster, envir = mf)
  
  # Create matrices X, U, Z0, Z_1,... from data formulas
  
  X <- if (is.null(mean_mods)) NULL else do.call(model.matrix, list(object = mean_mods, data = mf))
  U <- if (is.null(var_mods)) NULL else do.call(model.matrix, list(object = var_mods, data = mf))
  Z0 <- if (is.null(sel_zero_mods)) NULL else do.call(model.matrix, list(object = sel_zero_mods, data = mf))
  if (is.null(sel_mods)) {
    Z <- NULL
  } else {
    if (is.list(sel_mods)) {
      if (length(sel_mods) != length(steps)) stop("sel_mods must be a list with length equal to the number of steps.")
      Z <- lapply(sel_mods, model.matrix, data = mf)
    } else {
      Z1 <- do.call(model.matrix, list(object = sel_mods, data = mf))
      Z <- rep(list(Z1), times = length(steps))
    }
  }
  
  res <- fit_selection_model(
    yi = yi, sei = sei, pi = pi, ai = ai, cluster = cluster, 
    X = X, U = U, Z0 = Z0, Z = Z,
    steps = steps,
    make_sandwich = make_sandwich, 
    selection_type = selection_type,
    estimator = estimator, 
    theta = theta,
    optimizer = optimizer, 
    optimizer_control = optimizer_control,
    use_jac = use_jac
  ) 
  
  if ((make_sandwich == "raw") && (bootstrap=="none")) {
    return(res)
  }
  
  # build selmodel object
  
  res$est <- data.frame(
    estimator = estimator,
    param = names(res$est), 
    Est = as.numeric(res$est), 
    SE = sqrt(diag(res$vcov))
  )
  
  if ("large-sample" %in% boot_CI) {
    qz <- qnorm(1 - (1 - conf_level) / 2)
    res$est$CI_lo <- res$est$Est - qz * res$est$SE
    res$est$CI_hi <- res$est$Est + qz * res$est$SE
  } 

  if (selection_type == "step") {
    class(res) <- c("step.selmodel","selmodel")
  } else if (selection_type == "beta") {
    class(res) <- c("beta.selmodel","selmodel")
  }
  
  if (bootstrap != "none") {
    
    boot_sandwich <- any("student" %in% boot_CI)
    
    booties_df <- pbapply::pbreplicate(n = max(R), {
      bootstrap_selmodel(
        yi = yi, sei = sei, pi = pi, ai = ai, cluster = cluster, 
        X = X, U = U, Z0 = Z0, Z = Z,
        steps = steps,
        make_sandwich = boot_sandwich, 
        selection_type = selection_type,
        estimator = estimator, 
        theta = theta,
        optimizer = optimizer, 
        optimizer_control = optimizer_control,
        use_jac = use_jac,
        wtype = bootstrap
      )
    }, simplify = FALSE)
    
    res$bootstrap_reps <- do.call(rbind, booties_df)
    res$est$bootstrap <- bootstrap
    
    if (any(c("percentile","basic","student") %in% boot_CI)) {
      
      est <- split(res$est$Est, res$est$param)
      se <- split(res$est$SE, res$est$param)
      boots <- by(res$bootstrap_reps, res$bootstrap_reps$param, identity)
      
      boot_CIs <- mapply(
        \(e, s, b) simhelpers::bootstrap_CIs(
          boot_est = b$Est, boot_se = b$SE, 
          est = e, se = s, 
          CI_type = boot_CI, level = conf_level, B_vals = R, ...
        ),
        e = est,
        s = se,
        b = boots,
        SIMPLIFY = FALSE
      )
      
      boot_lengths <- sapply(boot_CIs, nrow)
      if (all(boot_lengths == 1L)) {
        boot_CIs <- do.call(rbind, c(boot_CIs, make.row.names = FALSE))
        res$est <- cbind(res$est, boot_CIs)
      } else {
        res$est$boot_CIs <- boot_CIs
      }
    }
    
    class(res) <- c("boot.selmodel", class(res))
    
  }
  
  res$cl <- cl
  res$mf <- mf
  
  return(res)
}
