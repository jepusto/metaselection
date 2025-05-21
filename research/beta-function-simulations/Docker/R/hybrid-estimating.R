#-------------------------------------------------------------------------------
# Log likelihood function and selection constraint function

step_weighted_logpartlik <- function(
    theta,                                      # full parameter vector
    yi,                                         # outcome vector
    sei,                                        # sampling standard errors
    pi = pnorm(yi / sei, lower.tail = FALSE),   # one-sided p-values
    ai = NULL,                                  # analytic weight 
    beta = NULL,                                # mean parameter coefficients
    gamma = NULL,                               # variance component coefficients
    omega0 = NULL,                              # selection model coefficients
    omega = NULL,                               # selection model coefficients
    steps = .025,                               # steps / cut-points
    X = NULL,                                   # mean parameter design matrix
    U = NULL,                                   # variance component design matrix
    Z0 = NULL,                                  # selection model design matrix for highest step
    Z = NULL,                                   # selection model design matrices for each cut-point
    contributions = FALSE,                      # not used
    negate = FALSE                              # not used
) {
  
  params <- parse_step_params(
    theta = theta,
    yi = yi,
    sei = sei,
    pi = pi,
    beta = beta,
    gamma = gamma,
    omega0 = omega0,
    omega = omega,
    steps = steps,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z
  )
  
  # likelihood contributions
  log_lik_i <- ((yi - params$mu)^2 / (-2 * params$eta) - log(params$eta) / 2 ) / params$weight_vec

  # weighted log likelihood (Eq. 51)
  log_lik <- if (is.null(ai)) sum(log_lik_i) else sum(ai * log_lik_i)
  
  return(log_lik)
  
}


step_selection_constraint <- function(
    theta,                                      # full parameter vector
    yi,                                         # outcome vector
    sei,                                        # sampling standard errors
    pi = pnorm(yi / sei, lower.tail = FALSE),   # one-sided p-values
    ai = NULL,                                  # analytic weight
    beta = NULL,                                # mean parameter coefficients
    gamma = NULL,                               # variance component coefficients
    omega0 = NULL,                              # selection model coefficients
    omega = NULL,                               # selection model coefficients
    steps = .025,                               # steps / cut-points
    X = NULL,                                   # mean parameter design matrix
    U = NULL,                                   # variance component design matrix
    Z0 = NULL,                                  # selection model design matrix for highest step
    Z = NULL,                                   # selection model design matrices for each cut-point
    contributions = FALSE,                      # not used
    negate = FALSE                              # not used
) {
  
  params <- parse_step_params(
    theta = theta,
    yi = yi,
    sei = sei,
    pi = pi,
    beta = beta,
    gamma = gamma,
    omega0 = omega0,
    omega = omega,
    steps = steps,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z,
    calc_Ai = TRUE
  )

  scores <- calculate_S_omega_ij(
    B_mat = params$B_mat,
    Ai = params$Ai,
    weight_vec = params$weight_vec,
    cats = params$cats,
    Z0 = Z0,
    Z = Z,
    lambda0 = params$lambda0,
    lambda = params$lambda
  )
  
  S_omega <- if (is.null(Z0)) scores$S_omega_ij else cbind(scores$S_omega0_ij, scores$S_omega_ij)
  
  colSums(S_omega)
  
}

#-------------------------------------------------------------------------------
# Hybrid score function, profiled in beta

step_hybrid_profile_score <- function(
    theta,                                      # parameter vector excluding beta
    yi,                                         # outcome vector
    sei,                                        # sampling standard errors
    pi = pnorm(yi / sei, lower.tail = FALSE),   # one-sided p-values
    ai = NULL,                                  # analytic weight
    steps = .025,                               # steps / cut-points
    X = NULL,                                   # mean parameter design matrix
    U = NULL,                                   # variance component design matrix
    Z0 = NULL,                                  # selection model design matrix for highest step
    Z = NULL                                   # selection model design matrices for each cut-point
) {
  
  if (is.null(X)) {
    x_dim <- 1L
    theta_na <- c(NA_real_, theta)
  } else {
    x_dim <- ncol(X)
    theta_na <- c(rep(NA_real_, x_dim), theta)
  }
  
  scores <- step_hybrid_score(
    theta = theta_na,
    yi = yi, sei = sei, pi = pi, ai = ai,
    steps = steps, 
    X = X, U = U, Z0 = Z0, Z = Z
  )  
  scores[-(1:x_dim)]
}

#-------------------------------------------------------------------------------
# Hybrid score function

step_hybrid_score <- function(
  theta,                                      # full parameter vector
  yi,                                         # outcome vector
  sei,                                        # sampling standard errors
  pi = pnorm(yi / sei, lower.tail = FALSE),   # one-sided p-values
  ai = NULL,                                  # analytic weight
  beta = NULL,                                # mean parameter coefficients
  gamma = NULL,                               # variance component coefficients
  omega0 = NULL,                              # selection model coefficients
  omega = NULL,                               # selection model coefficients
  steps = .025,                               # steps / cut-points
  X = NULL,                                   # mean parameter design matrix
  U = NULL,                                   # variance component design matrix
  Z0 = NULL,                                  # selection model design matrix for highest step
  Z = NULL,                                   # selection model design matrices for each cut-point
  contributions = FALSE,                      # whether to return matrix of score contributions,
  negate = FALSE                              # whether to return the negative of the scores
) {

  params <- parse_step_params(
    theta = theta,
    yi = yi,
    sei = sei,
    pi = pi,
    beta = beta,
    gamma = gamma,
    omega0 = omega0,
    omega = omega,
    steps = steps,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z,
    calc_Ai = TRUE
  )


  #----------------------------------------------------------
  # score contributions (Eq. 12, just for observation ij)
  zeta_i <- (yi - params$mu) / params$eta

  # score contribution for mu
  S_beta_ij  <-  zeta_i / params$weight_vec
  if (!is.null(X)) S_beta_ij <- X * S_beta_ij

  # score contribution for gamma
  S_gamma_ij <- params$tausq * (zeta_i^2 - 1 / params$eta) / (2 * params$weight_vec)
  if (!is.null(U)) S_gamma_ij <- U * S_gamma_ij

  scores <- calculate_S_omega_ij(
    B_mat = params$B_mat,
    Ai = params$Ai,
    weight_vec = params$weight_vec,
    cats = params$cats,
    Z0 = Z0,
    Z = Z,
    lambda0 = params$lambda0,
    lambda = params$lambda
  )
  
  S_omega_ij <- scores$S_omega_ij
  S_omega0_ij <- if (is.null(Z0)) NULL else scores$S_omega0_ij

  score_contributions <- get_score_contributions(
    S_beta_ij = S_beta_ij,
    S_gamma_ij = S_gamma_ij,
    S_omega_ij = S_omega_ij,
    S_omega0_ij = S_omega0_ij,
    ai = ai,
    contributions = contributions
  )
  
  if (negate) return(-score_contributions) else return(score_contributions)
  
}



#-------------------------------------------------------------------------------
# Jacobian function

step_hybrid_profile_jacobian <- function(
    theta,                                      # full parameter vector
    yi,                                         # outcome vector
    sei,                                        # sampling standard errors
    pi = pnorm(yi / sei, lower.tail = FALSE),   # one-sided p-values  
    ai = NULL,                                  # analytic weight
    steps = .025,                               # steps / cut-points
    X = NULL,                                   # mean parameter design matrix
    U = NULL,                                   # variance component design matrix
    Z0 = NULL,                                  # selection model design matrix for highest step
    Z = NULL                                    # selection model design matrices for each cut-point
) {
  
  if (is.null(X)) {
    x_dim <- 1L
    theta_na <- c(NA_real_, theta)
  } else {
    x_dim <- ncol(X)
    theta_na <- c(rep(NA_real_, x_dim), theta)
  }
  
  jac <- step_hybrid_jacobian(
    theta = theta_na,
    yi = yi, sei = sei, pi = pi, ai = ai,
    steps = steps, 
    X = X, U = U, Z0 = Z0, Z = Z
  )
  
  jac[-(1:x_dim),-(1:x_dim)]
  
}

step_hybrid_jacobian <- function(
    theta,                                      # full parameter vector
    yi,                                         # outcome vector
    sei,                                        # sampling standard errors
    pi = pnorm(yi / sei, lower.tail = FALSE),   # one-sided p-values  
    ai = NULL,                                  # analytic weight
    beta = NULL,                                # mean parameter coefficients
    gamma = NULL,                               # variance component coefficients
    omega0 = NULL,                              # selection model coefficients
    omega = NULL,                               # selection model coefficients
    steps = .025,                               # steps / cut-points
    X = NULL,                                   # mean parameter design matrix
    U = NULL,                                   # variance component design matrix
    Z0 = NULL,                                  # selection model design matrix for highest step
    Z = NULL                                    # selection model design matrices for each cut-point
) {
  
  params <- parse_step_params(
    theta = theta,
    yi = yi,
    sei = sei,
    pi = pi,
    beta = beta,
    gamma = gamma,
    omega0 = omega0,
    omega = omega,
    steps = steps,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z,
    calc_Ai = TRUE
  )
  
  # analytic weights
  if (is.null(ai)) ai <- rep(1, params$k)
  
  
  #----------------------------------------------------------
  
  zeta_i <- (yi - params$mu) / params$eta
  
  # J_beta_beta is x_dim * x_dim
  J_beta_beta_right <-  1 / (params$eta * params$weight_vec)
  J_beta_beta <- -1 * matrix_diag_crossprod(d = ai * J_beta_beta_right, A = X, B = X)
  
  # J_beta_gamma is x_dim * u_dim
  J_beta_gamma_right <- params$tausq * zeta_i / (params$eta * params$weight_vec)
  J_beta_gamma <- -1 * matrix_diag_crossprod(d = ai * J_beta_gamma_right, A = X, B = U)

  # J_gamma_gamma is u_dim * u_dim
  J_gamma_gamma_right <- (params$tausq^2 * (zeta_i^2 / params$eta - 1 / (2 * params$eta^2)) - 
                            params$tausq * (zeta_i^2 / 2 - 1 / (2 * params$eta))) / params$weight_vec 
  J_gamma_gamma <- -1 * matrix_diag_crossprod(d = ai * J_gamma_gamma_right, A = U, B = U)
  
  if (is.null(Z0)) {
    J_beta_omega_right <- params$lambda * zeta_i / params$weight_vec^2
    J_gamma_omega_right <-  params$tausq * params$lambda * (zeta_i^2 - 1 / params$eta) / (2 * params$weight_vec^2)
    
    Z_full <- Z
    z_index <- c(0,cumsum(params$z_dim))
    
  } else {
    J_beta_omega_right <- params$lambda_full * zeta_i / params$weight_vec^2
    J_gamma_omega_right <-  params$tausq * params$lambda_full * (zeta_i^2 - 1 / params$eta) * (2 * params$weight_vec^2)
    
    Z_full <- if (is.list(Z)) c(list(Z0),Z) else list(Z0, Z)
    z_index <- c(0,cumsum(c(params$z0_dim, params$z_dim)))
    
  }
  
  if (is.null(Z)) {
    if (is.null(Z0)) {
      J_beta_omega <- -1 * matrix_diag_crossprod(A = X, B = ai * J_beta_omega_right)
      J_gamma_omega <- -1 * matrix_diag_crossprod(A = U, B = ai * J_gamma_omega_right)
      
    } else {
      
      J_beta_omega_right <- cbind(J_beta_omega_right[,1] * Z0, J_beta_omega_right[,-1])
      J_beta_omega <- -1 * matrix_diag_crossprod(A = X, B = ai * J_beta_omega_right)
      
      J_gamma_omega_right <- cbind(J_gamma_omega_right[,1] * Z0, J_gamma_omega_right[,-1])
      J_gamma_omega <- -1 * matrix_diag_crossprod(A = U, B = ai * J_gamma_omega_right)
      
    }
    
  } else if (is.matrix(Z_full)) {
    
    J_beta_omega <- -1 * matrix_diag_crossprod(A = X, d = ai * J_beta_omega_right, B = Z_full)
    J_gamma_omega <- -1 * matrix_diag_crossprod(A = U, d = ai * J_gamma_omega_right, B = Z_full)
    
  } else {

    J_beta_omega_right_h <-  apply(ai * J_beta_omega_right, 2, identity, simplify = FALSE)
    J_beta_omega_h <- mapply(matrix_diag_crossprod, d = J_beta_omega_right_h, B = Z_full, MoreArgs = list(A = X), SIMPLIFY = FALSE)
    J_beta_omega <- -1 * do.call(cbind, J_beta_omega_h)
    
    J_gamma_omega_right_h <- apply(ai * J_gamma_omega_right, 2, identity, simplify = FALSE)
    J_gamma_omega_h <- mapply(matrix_diag_crossprod, d = J_gamma_omega_right_h, B = Z_full, MoreArgs = list(A = U), SIMPLIFY = FALSE)
    J_gamma_omega <- -1 * do.call(cbind, J_gamma_omega_h)
    
  }
  
  
  H_omega <- calculate_H_omega(
    B_mat = params$B_mat,
    Ai = params$Ai,
    weight_vec = params$weight_vec,
    cats = params$cats,
    eta = params$eta,
    c_mat = params$c_mat,
    k = params$k,
    H = params$H,
    tausq = params$tausq,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z,
    ai = ai,
    z0_dim = params$z0_dim,
    z_dim = params$z_dim,
    lambda0 = params$lambda0,
    lambda = params$lambda,
    lambda_full = params$lambda_full
  )      

  H_omega_ <- cbind(t(H_omega$H_beta_omega), t(H_omega$H_gamma_omega), H_omega$H_omega_omega)
                         

  # Put the whole Hessian together
  
  J_matrix <- rbind(
    cbind(J_beta_beta, J_beta_gamma, J_beta_omega),
    cbind(t(J_beta_gamma), J_gamma_gamma, J_gamma_omega),
    H_omega_
  )
  colnames(J_matrix) <- rownames(J_matrix) <- params$H_names
  
  return(J_matrix)
  
}


step_selection_constraint_grad <- function(
    theta,                                      # full parameter vector
    yi,                                         # outcome vector
    sei,                                        # sampling standard errors
    pi = pnorm(yi / sei, lower.tail = FALSE),   # one-sided p-values
    ai = NULL,                                  # analytic weight
    beta = NULL,                                # mean parameter coefficients
    gamma = NULL,                               # variance component coefficients
    omega0 = NULL,                              # selection model coefficients
    omega = NULL,                               # selection model coefficients
    steps = .025,                               # steps / cut-points
    X = NULL,                                   # mean parameter design matrix
    U = NULL,                                   # variance component design matrix
    Z0 = NULL,                                  # selection model design matrix for highest step
    Z = NULL,                                   # selection model design matrices for each cut-point
    contributions = FALSE                       # not used
) {
 
  params <- parse_step_params(
    theta = theta,
    yi = yi,
    sei = sei,
    pi = pi,
    beta = beta,
    gamma = gamma,
    omega0 = omega0,
    omega = omega,
    steps = steps,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z,
    calc_Ai = TRUE
  )
  
  # analytic weights
  if (is.null(ai)) ai <- rep(1, params$k)
  
  H_omega <- calculate_H_omega(
    B_mat = params$B_mat,
    Ai = params$Ai,
    weight_vec = params$weight_vec,
    cats = params$cats,
    eta = params$eta,
    c_mat = params$c_mat,
    k = params$k,
    H = params$H,
    tausq = params$tausq,
    X = X,
    U = U,
    Z0 = Z0,
    Z = Z,
    ai = ai,
    z0_dim = params$z0_dim,
    z_dim = params$z_dim,
    lambda0 = params$lambda0,
    lambda = params$lambda,
    lambda_full = params$lambda_full
  )
  
  rbind(H_omega$H_beta_omega, H_omega$H_gamma_omega, H_omega$H_omega_omega)
}