#' 
#' @importFrom Formula as.Formula


build_model_frame <- function(
    data,
    yi,
    vi,
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
  if (missing(yi) || (missing(sei) & missing(vi))) stop("You must specify an effect size variable yi and a variance variable vi or a standard error variable sei.")
  if (!missing(sei) & !missing(vi)) stop("You must specify either a variance variable vi or a standard error variable sei, but not both.")
  yi_str <- deparse(substitute(yi))
  vi_str <- if (missing(vi)) NULL else deparse(substitute(vi))
  sei_str <- if(missing(sei)) NULL else deparse(substitute(sei))
  pi_str <- if (missing(pi)) NULL else deparse(substitute(pi))
  ai_str <- if (missing(ai)) NULL else deparse(substitute(ai))
  cl_str <- if (missing(cluster)) NULL else deparse(substitute(cluster))
  yi_sei <- reformulate(c(vi_str, sei_str, pi_str, ai_str, cl_str), response = yi_str)
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
  yi, 
  sei, 
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
    gamma_start <- log(max(1e-4, mean(ri_sq) - mean(sei^2)))
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
  subset = NULL,
  vcov_type = "robust",
  selection_type = "step",
  estimator = "CML",
  theta = NULL,
  optimizer = "BFGS",
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
    if (!is.null(Z)) Z <- lapply(Z, \(z) z[subset,,drop=FALSE])
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
  
  if (estimator %in% c("ML","CML")) {
    
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
    
    if (vcov_type == "raw") {
      return(mle_est)
    } 
    
    max_method <- row.names(mle_est)[which.max(mle_est$value)]
    theta_names <- 1:length(theta)
    theta <- as.numeric(mle_est[max_method, theta_names])
    info <- mle_est[max_method, -theta_names]
    names(theta) <- params$H_names
    
    if (any(is.na(theta))) stop("Could not obtain parameter estimates. Perhaps try a different optimizer?")
    
    if (vcov_type == "none") {
      return(theta)
    }
    
  } else if (estimator %in% c("ARGL-full","hybrid-full")) {
    
    info <- list()

    jac <- if (use_jac) step_hybrid_jacobian else NULL
    
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
    
    max_method <- if ("method" %in% names(nleqslv_args)) nleqslv_args$method else "default"
    max_method <- paste("nleqslv", max_method)
    theta <- nleqslv_res$x
    
    names(theta) <- params$H_names
    
    if (vcov_type == "raw") {
      return(list(est = theta, max_method = max_method, info = info))
    }
    
    if (any(is.na(theta))) stop("Could not obtain parameter estimates. Perhaps try a different optimizer?")
    
    if (vcov_type == "none") {
      return(theta)
    }
    
  }  else if (estimator %in% c("ARGL","hybrid")) {
    
    info <- list()
    
    x_index <- if (is.null(X)) 1L else 1:ncol(X)
    theta <- theta[-x_index]
     
    jac <- if (use_jac) step_hybrid_profile_jacobian else NULL
    
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
    
    max_method <- if ("method" %in% names(nleqslv_args)) nleqslv_args$method else "default"
    max_method <- paste("nleqslv", max_method)
    theta <- nleqslv_res$x
    
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
    
    if (vcov_type == "raw") {
      return(list(est = theta, max_method = max_method, info = info))
    }

    if (any(is.na(theta))) stop("Could not obtain parameter estimates. Perhaps try a different optimizer?")
    
    if (vcov_type == "none") {
      return(theta)
    }
  }

  # Compute composite log likelihood
  
  norm_const <-  log(2 * base::pi) * length(yi) / 2
  
  if (selection_type == "step") {
    
    log_lik <- step_loglik(
      theta = theta,  
      yi = yi, sei = sei, pi = pi, ai = ai,
      steps = steps,
      X = X, U = U, Z0 = Z0, Z = Z
    ) - norm_const
    
    wt_partial_log_lik <- step_weighted_logpartlik(
      theta = theta,  
      yi = yi, sei = sei, pi = pi, ai = ai,
      steps = steps,
      X = X, U = U, Z0 = Z0, Z = Z
    )
    
  } else if (selection_type == "beta") {
    
    log_lik <- beta_loglik(
      theta, 
      yi = yi, sei = sei, pi = pi, ai = ai,
      steps = steps,
      X = X, U = U
    ) - norm_const
    
    wt_partial_log_lik <- NULL
  }

  # Compute score contributions using MLEs
  
  if (estimator %in% c("ML","CML")) {
    
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
    
  
  } else if (estimator %in% c("ARGL-full","ARGL","hybrid-full","hybrid")) {
    
    scores <- step_hybrid_score(theta = theta, 
                                yi = yi, sei = sei, pi = pi, ai = ai,
                                steps = steps,
                                X = X, U = U, Z0 = Z0, Z = Z, contributions = TRUE)
    
    
  }

  
  # Compute Hessian using MLEs
  
  if (estimator %in% c("ML","CML")) {
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
    
  
  } else if (estimator %in% c("ARGL-full","ARGL","hybrid-full","hybrid")) {
    
    hess <- step_hybrid_jacobian(theta = theta, 
                                 yi = yi, sei = sei, pi = pi, ai = ai,
                                 steps = steps,
                                 X = X, U = U, Z0 = Z0, Z = Z)
  }
  
  hess_inv <- tryCatch(MASS::ginv(hess), error = function(e) e)
  if (inherits(hess_inv, "error")) {
    warning("Hessian matrix is not invertible. Cluster-robust standard errors cannot be computed.")
    hess_inv <- matrix(NA_real_, nrow = nrow(hess), ncol = ncol(hess))
  }
  
  if (vcov_type == "robust") {
    
    # Compute sandwich variance estimator using score and Hessian (Eq. 52-54)
    if (is.null(cluster)) {
      meat <- crossprod(scores)
    } else {
      score_j <- rowsum(scores, group = cluster)
      meat <- crossprod(score_j)
    }
    
    vcov_mat <- hess_inv %*% meat %*% t(hess_inv)

  } else {
    vcov_mat = -hess_inv
  }
  
  
  rownames(vcov_mat) <- names(theta)
  colnames(vcov_mat) <- names(theta)
  
  list(
    est = theta, 
    vcov = vcov_mat, 
    method = max_method,
    info = info,
    ll = log_lik,
    wpll = wt_partial_log_lik
  )
  
}

two_stage_sample <- function(cluster) {
  m <- nlevels(cluster)
  i <- 1L
  while (length(unique(i)) == 1) {
    i <- sample.int(n = m, replace = TRUE)
  }
  dim(i) <- c(1, m)
  bi <- apply(i, 1, tabulate, m)
  
  ki <- table(cluster)
  bi_list <- mapply(
    function(ki,a) tabulate(sample.int(ki, a * ki, replace = TRUE), nbins = ki),
    ki = ki, 
    a = bi[,1],
    SIMPLIFY = FALSE
  )
  
  unsplit(bi_list,cluster)
  
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
    vcov_type = "robust",
    selection_type = "step",
    estimator = "CML",
    theta = NULL,
    optimizer = "BFGS",
    optimizer_control = list(),
    use_jac = TRUE,
    wtype = c("two-stage","multinomial", "exponential"),
    retry = 0L
) {
  
  if (!is.null(cluster)) {
    cluster <- factor(cluster)
    m <- nlevels(cluster)
    cluster_numeric <- as.integer(cluster)
  } else {
    m <- length(yi)
    cluster_numeric <- 1:m
  }
  
  if (wtype == "two-stage") {
    
    wi <- two_stage_sample(cluster)
    
    if (!is.null(ai)) {
      wi <- ai * wi
    }
    
    non_zero_cl <- wi > 0
    cl_subset <- if (all(non_zero_cl)) NULL else non_zero_cl
    
  } else {
    
    if (wtype == "multinomial") {
      i <- 1L
      while (length(unique(i)) == 1) {
        i <- sample.int(n = m, replace = TRUE)
      }
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
    
  }
  
  res <- suppressWarnings(tryCatch(
    fit_selection_model(
      yi = yi, 
      sei = sei, 
      pi = pi, 
      steps = steps,
      ai = wi, 
      cluster = cluster, 
      subset = cl_subset,
      X = X, U = U, Z0 = Z0, Z = Z,
      vcov_type = vcov_type, 
      selection_type = selection_type,
      estimator = estimator, 
      theta = theta,
      optimizer = optimizer, 
      optimizer_control = optimizer_control,
      use_jac = use_jac
    ),
    error = \(e) e
  ))
  
  if (inherits(res, "error") || (vcov_type != "none" && any(is.na(diag(res$vcov))))) {
    if (retry == 0L) return(NULL)
    
    res <- bootstrap_selmodel(
      yi = yi, 
      sei = sei, 
      pi = pi, 
      steps = steps, 
      ai = ai,
      cluster = cluster, 
      X = X, 
      U = U, 
      Z0 = Z0, 
      Z = Z, 
      vcov_type = vcov_type,
      selection_type = selection_type,
      estimator = estimator,
      theta = theta,
      optimizer = optimizer,
      optimizer_control = optimizer_control,
      use_jac = use_jac,
      wtype = wtype,
      retry = retry - 1L
    )
    return(res)
  }
  
  if (vcov_type != "none") {
    est <- data.frame(
      param = names(res$est),
      Est = as.numeric(res$est),
      SE = as.numeric(sqrt(diag(res$vcov)))
    )
  } else {
    est <- data.frame(
      param = names(res),
      Est = as.numeric(res)
    )
  }
  
  return(est)
}



jackknife_selmodel <- function(
    est,
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
    selection_type = "step",
    estimator = "CML",
    optimizer = "BFGS",
    optimizer_control = list(),
    use_jac = TRUE
) {
  
  if (is.null(cluster)) {
    cluster_jk <- 1:length(yi)
    index_jk <- cluster_jk
  } else {
    cluster_jk <- cluster
    index_jk <- unique(cluster)
  }
  
  q <- progressr::progressor(length(index_jk))
  
  jack_vals <- future.apply::future_lapply(
    index_jk, \(i) {
      q()
      res_i <- suppressWarnings(tryCatch(
        fit_selection_model(
          yi = yi, sei = sei, pi = pi, ai = ai, cluster = cluster, 
          X = X, U = U, Z0 = Z0, Z = Z,
          subset = cluster_jk != i,
          steps = steps,
          vcov_type = "none", 
          selection_type = selection_type,
          estimator = estimator, 
          theta = est,
          optimizer = optimizer, 
          optimizer_control = optimizer_control,
          use_jac = use_jac
        ),
        error = \(e) NULL
      ))
      res_i
    })
  
  do.call(rbind, jack_vals)
}

#' @title Estimate step or beta selection model
#'
#' @description Estimate step or beta selection model, with standard errors and
#'   confidence intervals based on either cluster-robust variance estimators
#'   (i.e., sandwich estimators) or cluster-level bootstrapping to handle
#'   dependent effect size estimates.
#'
#'
#' @param data \code{data.frame} or \code{tibble} containing the meta-analytic
#'   data
#' @param yi vector of effect sizes estimates.
#' @param vi vector of sample variances. If \code{vi} is specified, then the
#'   \code{sei} argument must be omitted.
#' @param sei vector of sampling standard errors. If \code{sei} is specified,
#'   then the \code{vi} argument must be omitted.
#' @param pi optional vector of one-sided p-values. If not specified, p-values
#'   will be computed from \code{yi} and \code{sei}.
#' @param ai optional vector of analytic weights.
#' @param cluster vector indicating which observations belong to the same
#'   cluster.
#' @param selection_type character string specifying the type selection model to
#'   estimate, with possible options \code{"step"} or \code{"beta"}.
#' @param steps If \code{selection_type = "step"}, a numeric vector of one or
#'   more values specifying the thresholds (or steps) where the selection
#'   probability changes, with a default of \code{steps = .025}. If
#'   \code{selection_type = "beta"}, then a numeric vector of two values
#'   specifying the thresholds beyond which the selection function is truncated,
#'   with a default of \code{steps = c(.025, .975)}.
#' @param mean_mods optional model formula for moderators related to average
#'   effect size magnitude.
#' @param var_mods optional model formula for moderators related to effect size
#'   heterogeneity.
#' @param sel_mods optional model formula for moderators related to the
#'   probability of selection. Only relevant for \code{selection_type = "step"}.
#' @param sel_zero_mods optional model formula for moderators related to the
#'   probability of selection for p-values below the lowest threshold value of
#'   \code{steps}. Only relevant for \code{selection_type = "step"}.
#' @param subset optional logical expression indicating a subset of observations
#'   to use for estimation.
#' @param estimator vector indicating whether to use the composite marginal
#'   likelihood estimator (option \code{"CML"}) or the augmented and reweighted
#'   Gaussian likelihood estimator (option \code{"ARGL"} or \code{"ARGL-full"}).
#'   If \code{selection_type = "beta"}, only the composite marginal likelihood
#'   estimator, \code{"CML"}, is available. For step function models, both
#'   estimators are available.
#' @param vcov_type character string specifying the type of variance-covariance
#'   matrix to calculate, with possible options \code{"robust"} for robust or
#'   cluster-robust standard errors, \code{"model-based"} for model-based
#'   standard errors, or \code{"none"}.
#' @param CI_type character string specifying the type of confidence interval to
#'   calculate, with possible options \code{"large-sample"} for large-sample
#'   normal interval (the default), \code{"percentile"} for a percentile
#'   interval, \code{"BCa"} for a bias-corrected-and-accelerated interval,
#'   \code{"bias-corrected"} for a bias-corrected percentile interval (without
#'   acceleration correction), \code{"normal"} for a standard normal interval,
#'   \code{"basic"} for a basic interval, \code{"student"} for a studentized
#'   interval, or \code{"none"}.
#' @param conf_level desired coverage level for confidence intervals, with the
#'   default value set to \code{.95}.
#' @param theta optional numeric vector of starting values to use in
#'   optimization routines.
#' @param optimizer character string indicating the optimizer to use. Ignored if
#'   \code{estimator = "ARGL"} or \code{"ARGL-full"}.
#' @param optimizer_control an optional list of control parameters to be used
#'   for optimization
#' @param use_jac logical with \code{TRUE} (the default) indicating to use the
#'   Jacobian of the estimating equations for optimization.
#' @param bootstrap character string specifying the type of bootstrap to run,
#'   with possible options \code{"none"} (the default), \code{"exponential"} for
#'   the fractionally re-weighted cluster bootstrap, \code{"multinomial"} for
#'   a conventional clustered bootstrap, or , \code{"two-stage"} for
#'   a two-stage clustered bootstrap.
#' @param R number of bootstrap replications, with a default of \code{1999}.
#' @param retry_bootstrap number of times to re-draw a bootstrap sample in the event of non-convergence, with a default of \code{0}.
#' @param ... further arguments passed to \code{simhelpers::bootstrap_CIs}.
#'
#' @returns An object of class \code{"selmodel"} containing the following
#'   components:
#' \describe{
#'   \item{\code{est}}{A data frame with parameter estimates, standard errors, and
#'   confidence intervals. Note that the results do not include p-values so
#'   as to focus interpretation on the parameter estimates, rather than on
#'   the statistical significance of any given parameter.}
#'   \item{\code{vcov}}{A matrix containing the estimated variance-covariance matrix
#'   of the parameter estiamtes}
#'   \item{\code{method}}{Character string indicating the optimization method used to solve for parameter estimates.}
#'   \item{\code{info}}{Further informaton about the optimization results.}
#'   \item{\code{ll}}{Log likelihood of the model evaluated at the reported parameter estimates.}
#'   \item{\code{wpll}}{Weighted partial log likelihood of the random effects model, with weights corresponding to inverse selection probabilities}
#'   \item{\code{n_clusters}}{Number of independent clusters of effect sizes.}
#'   \item{\code{n_effects}}{Number of effect size estimates in the data.}
#'   \item{\code{...}}{Some additional elements containing information about the methods used to estimate the model.}
#' }
#'
#' @export
#'
#' @examples
#' res_ML <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   steps = 0.025,
#'   estimator = "CML",
#'   bootstrap = "none"
#' )
#'
#' res_ML
#' summary(res_ML)
#'
#' # configure progress bar
#' progressr::handlers(global = TRUE)
#'
#' res_hybrid <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   steps = 0.025,
#'   estimator = "ARGL",
#'   bootstrap = "multinomial",
#'   CI_type = "percentile",
#'   R = 19
#' )
#'
#' res_hybrid
#' summary(res_hybrid)
#' 


selection_model <- function(
    data,
    yi,
    vi,
    sei,
    pi,
    ai,
    cluster,
    selection_type = c("step","beta"),
    steps = NULL,
    mean_mods = NULL,
    var_mods = NULL,
    sel_mods = NULL,
    sel_zero_mods = NULL,
    subset = NULL,
    estimator = "CML",
    vcov_type = "robust",
    CI_type = "large-sample",
    conf_level = .95,
    theta = NULL,
    optimizer = NULL,
    optimizer_control = list(),
    use_jac = TRUE,
    bootstrap = "none",
    R = 1999,
    retry_bootstrap = 0L,
    ...
) {
  
  selection_type <- match.arg(selection_type)
  estimator <- match.arg(estimator, c("CML","ML","ARGL","ARGL-full","hybrid","hybrid-full"))
  vcov_type <- match.arg(vcov_type, c("model-based","robust","none","raw"))
  bootstrap <- match.arg(bootstrap, c("none","exponential","multinomial","two-stage"))
  CI_type <- match.arg(CI_type, c("large-sample","percentile","normal","student","basic","bias-corrected","BCa", "none"), several.ok = TRUE)  

  if (vcov_type == "model-based") {
    if (!(estimator %in% c("ML","CML"))) stop("vcov_type = 'model-based' is only allowed for estimator = 'CML'.")
    if (!missing(cluster)) stop("vcov_type = 'model-based' does not allow the use of a clustering variable.")
  }
  if (bootstrap %in% c("exponential","multinomial","two-stage")) {
    if (identical(as.integer(R), 0L)) stop("Bootstrap methods require setting R > 0.") 
  }
  if (any(c("percentile","normal","basic","student","bias-corrected","BCa") %in% CI_type)) {
    if (identical(as.integer(R), 0L)) stop("Bootstrap confidence intervals require setting R > 0.")
    if (bootstrap == "none") stop("Bootstrap confidence intervals require setting bootstrap to 'multinomial' or 'exponential'.")
  }

  if (is.null(optimizer)) {
    if (selection_type == "step") {
      optimizer <- if (estimator %in% c("ML","CML")) "Rvmmin" else "nleqslv"
    } else {
      optimizer <- "nlminb"
    }
  }
  
  if (missing(steps)) {
    steps <- if (selection_type == "step") .025 else c(.025, .975)
  }
  
  if (!inherits(steps, "numeric") || min(steps) <= 0 || max(steps) >= 1) stop("steps must be a numeric vector with all entries in the interval (0,1).")
  if (selection_type == "beta") {
    if (length(steps) != 2L) stop("steps must be a numeric vector of length 2 when selection_type = 'beta'.")
    if (!is.null(sel_mods)) stop("sel_mods must be NULL when selection_type = 'beta'.")
    if (!is.null(sel_zero_mods)) stop("sel_zero_mods must be NULL when selection_type = 'beta'.")
    if (!(estimator %in% c("ML","CML"))) stop("estimator must be equal to 'CML' when selection_type = 'beta'.")
  }
  
  # Create common model frame
  
  cl <- match.call()
  m <- match(c("data","yi", "vi", "sei", "pi", "ai", "cluster","subset", "mean_mods", 
               "var_mods", "sel_mods", "sel_zero_mods"), names(cl), 0L)
  mf <- cl[c(1L, m)]
  mf[[1L]] <- str2lang("metaselection:::build_model_frame")
  mf <- eval(mf, parent.frame())
  
  # Evaluate yi, sei, pi, ai from model frame
  
  yi <- eval(cl$yi, envir = mf)
  
  vi <- if (missing(vi)) NULL else eval(cl$vi, envir = mf)
  sei <- if (missing(sei)) sqrt(vi) else eval(cl$sei, envir = mf)
  
  pi <- if (missing(pi)) pnorm(yi / sei, lower.tail = FALSE) else eval(cl$pi, envir = mf)
  ai <- if (missing(ai)) NULL else eval(cl$ai, envir = mf)
  cluster <- if (missing(cluster)) NULL else eval(cl$cluster, envir = mf)
  
  if (!is.null(cluster)) {
    n_clusters <- length(unique(cluster))
    if (n_clusters < 3L) stop("The data contain fewer than three unique clusters. This method will not provide valid results.")
  }
  
  # Create matrices X, U, Z0, Z_1,... from data formulas
  
  X <- if (is.null(mean_mods)) NULL else do.call(model.matrix, list(object = mean_mods, data = mf))
  U <- if (is.null(var_mods)) NULL else do.call(model.matrix, list(object = var_mods, data = mf))
  Z0 <- if (is.null(sel_zero_mods)) NULL else {
    sel_zero_terms <- stats::terms(sel_zero_mods)
    attr(sel_zero_terms, "intercept") <- 0L
    do.call(model.matrix, list(object = sel_zero_terms, data = mf))
  }
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
    vcov_type = vcov_type, 
    selection_type = selection_type,
    estimator = estimator, 
    theta = theta,
    optimizer = optimizer, 
    optimizer_control = optimizer_control,
    use_jac = use_jac
  ) 
  
  if ((vcov_type == "raw") && (bootstrap=="none")) {
    return(res)
  }
  
  # build selmodel object
  
  res$est <- data.frame(
    estimator = estimator,
    param = names(res$est), 
    Est = as.numeric(res$est), 
    SE = sqrt(diag(res$vcov))
  )
  
  if ("large-sample" %in% CI_type) {
    qz <- qnorm(1 - (1 - conf_level) / 2)
    
    res$est$p_value <- 2 * pnorm(abs(res$est$Est / res$est$SE), lower.tail = FALSE)
    
    if (is.null(var_mods)) {
      
      res$est["gamma", "p_value"] <- NA
      
    } else {
      
      gamma_intercept <- grepl("gamma_\\(Intercept\\)", res$est$param)
      res$est[gamma_intercept, "p_value"] <- NA
      
    }
    

    res$est$CI_lo <- res$est$Est - qz * res$est$SE
    res$est$CI_hi <- res$est$Est + qz * res$est$SE
  } 

  if (selection_type == "step") {
    class(res) <- c("step.selmodel","selmodel")
  } else if (selection_type == "beta") {
    class(res) <- c("beta.selmodel","selmodel")
  }
  
  # bootstrap calculations
  
  if (bootstrap != "none") {
    
    if ("student" %in% CI_type) {
      boot_sandwich <- "robust"
    } else {
      boot_sandwich <- "none"
    }

    reps <- max(R)
    p <- progressr::progressor(reps)
    booties_df <- future.apply::future_replicate(reps, {
      p()
      bootstrap_selmodel(
        yi = yi, sei = sei, pi = pi, ai = ai, cluster = cluster, 
        X = X, U = U, Z0 = Z0, Z = Z,
        steps = steps,
        vcov_type = boot_sandwich, 
        selection_type = selection_type,
        estimator = estimator, 
        theta = theta,
        optimizer = optimizer, 
        optimizer_control = optimizer_control,
        use_jac = use_jac,
        wtype = bootstrap,
        retry = retry_bootstrap
      )
    }, simplify = FALSE, future.seed = TRUE)
    
    res$R <- sum(!sapply(booties_df, is.null))
    res$bootstrap_reps <- do.call(rbind, booties_df)
    res$bootstrap_type <- bootstrap
    
    # Jackknife estimate of empirical influence
    if ("BCa" %in% CI_type) {
      res$jack_vals <- jackknife_selmodel(
        est = res$est$Est,
        yi = yi, sei = sei, pi = pi, ai = ai, cluster = cluster, 
        X = X, U = U, Z0 = Z0, Z = Z,
        steps = steps,
        selection_type = selection_type,
        estimator = estimator, 
        optimizer = optimizer, 
        optimizer_control = optimizer_control,
        use_jac = use_jac
      )
    }
    
    if (any(c("percentile","normal","basic","student","bias-corrected","BCa") %in% CI_type)) {
      
      boot_CIs <- get_boot_CIs(
        bmod = res, 
        CI_type = CI_type, conf_level = conf_level, 
        R = pmin(R, res$R), ...
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
  
  # Evaluate predictions
  if (selection_type == "step") {
    predictions <- parse_step_params(
      theta = res$est$Est,
      yi = yi, sei = sei, pi = pi, ai = ai, 
      steps = steps, 
      X = X, U = U, Z0 = Z0, Z = Z, calc_Ai = TRUE
    )
    sel_params <- sum(predictions$z_dim)
    sel_zero_params <- predictions$z0_dim
    
  } else if (selection_type == "beta") {
    predictions <- parse_beta_params(
      theta = res$est$Est,
      yi = yi, sei = sei, pi = pi, 
      alpha = steps, 
      X = X, U = U, calc_Ai = TRUE
    )
    sel_params <- 2L
    sel_zero_params <- 0L
    
  }
  
  # Finish building selmodel object
  
  res$cl <- cl
  res$mf <- mf
  res$selection_type <- selection_type
  res$steps <- steps
  res$estimator <- estimator
  res$vcov_type <- vcov_type
  res$conf_level <- conf_level
  res$param_dim <- c(
    mean = predictions$x_dim, 
    var = predictions$u_dim, 
    sel = sel_params, 
    sel_zero = sel_zero_params
  )
  res$predictions <- predictions[c("mu","tausq","eta","lambda_full","weight_vec","cats","B_mat","Ai")]
  
  if (!is.null(cluster)) res$n_clusters <- n_clusters else res$n_clusters <- NULL
  res$n_effects <- predictions$k
  
  res$ptable <- create_ptable(
    pvals = pi,
    cluster = cluster,
    steps = steps
  )
  
  res$selmods <- Z
  
  return(res)
}


get_boot_CIs <- function(bmod, CI_type, R, conf_level = 0.95, ...) {
  
  param_f <- factor(bmod$est$param, levels = bmod$est$param)
  est <- split(bmod$est$Est, param_f)
  se <- split(bmod$est$SE, param_f)
  boots <- by(
    bmod$bootstrap_reps, 
    factor(bmod$bootstrap_reps$param, levels = levels(param_f)), 
    identity
  )
  if (is.null(bmod$jack_vals)) {
    
    boot_CIs <- future.apply::future_mapply(
      \(e, s, b) simhelpers::bootstrap_CIs(
        boot_est = b$Est, boot_se = b$SE, 
        est = e, se = s, 
        CI_type = CI_type, level = conf_level, B_vals = R, ...
      ),
      e = est,
      s = se,
      b = boots,
      SIMPLIFY = FALSE,
      future.seed = TRUE
    )
  } else {
    inf_vals <- 
      bmod$jack_vals |>
      apply(1, \(x) bmod$est$Est - x) |>
      apply(1, identity, simplify = FALSE)
    
    boot_CIs <- future.apply::future_mapply(
      \(e, s, b, i) simhelpers::bootstrap_CIs(
        boot_est = b$Est, boot_se = b$SE, 
        est = e, se = s, influence = i,
        CI_type = CI_type, level = conf_level, B_vals = R, ...
      ),
      e = est,
      s = se,
      b = boots,
      i = inf_vals,
      SIMPLIFY = FALSE,
      future.seed = TRUE
    )
  }
  
  boot_CIs
  
}



create_ptable <- function(
    pvals,
    cluster,
    steps
) {
  
  steps <- c(0, steps, 1)
  
  psteps_l <- as.character(steps[-length(steps)])
  psteps_r <- as.character(steps[-1])
  len_l    <- nchar(psteps_l)
  psteps   <- paste0(psteps_l, " < p <= ", psteps_r)
  
  effects_group <- cut(pvals, breaks = steps, include.lowest = TRUE, labels = psteps)
  
  ptable   <- table(effects_group)
  ptable   <- data.frame(step = names(ptable), k = as.vector(ptable))
  
  if (!is.null(cluster)) {
  
    m <- tapply(cluster, effects_group, function(x) length(unique(x)))
    m[is.na(m)] <- 0L
    ptable$m <- m
    
    ptable <- ptable[, c("step", "m", "k")]
  
  }
  
  return(ptable)
  
}

