
#-------------------------------------------------------------------------------
# Log likelihood function

#' Estimate the log likelihood
#' 
#' @param theta full parameter vector
#' @param yi outcome vector
#' @param sei sampling standard error
#' @param pi one-sided p-values
#' @param ai analytic weight
#' @param beta mean parameter coefficients
#' @param gamma variance component coefficients
#' @param omega0 selection model coefficients
#' @param omega selection model coefficients
#' @param steps steps or cut-points
#' @param X mean parameter design matrix
#' @param U variance component design matrix
#' @param Z0 selection model design matrix for highest step
#' @param Z selection model design matrix for each cut-point
#' @returns A numeric vector.
# 

step_loglik <- function(
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
  
  # likelihood contributions
  log_lik_i <- log(params$weight_vec) - (yi - params$mu)^2 / (2 * params$eta) - log(params$eta) / 2 - log(params$Ai)
  
  # weighted log likelihood (Eq. 51)
  log_lik <- if (is.null(ai)) sum(log_lik_i) else sum(ai * log_lik_i)

  return(log_lik)
  
}

#-------------------------------------------------------------------------------
# Score function

step_score <- function(
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
    contributions = FALSE                       # whether to return matrix of score contributions
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
  
  k <- params$k
  H <- params$H
  
  #----------------------------------------------------------
  # Derivatives of A w/r/t mu, eta, lambda_h (Eq. 21-23)
  
  dB <- dB_dmu_eta(k = k, H = H, c_mat = params$c_mat)
  dA_dmu <- rowSums(dB$dmu * params$lambda_full) / sqrt(params$eta)
  dA_deta <- rowSums(dB$deta * params$lambda_full) / (2 * params$eta)
  
  dA_dlambda <- params$B_mat

  # Derivatives of weights w/r/t lambda_h (Eq. 20)
  dw_dlambda <- model.matrix(yi ~ 0 + params$cats)
  
  # Derivative of g w/r/t (z_hij * omega_h)
  # (For VHSM this is the same as g(z_hij * omega_h) = lambda_h.)
  
  #----------------------------------------------------------
  # score contributions (Eq. 12, just for observation ij)
  zeta_i <- (yi - params$mu) / params$eta
  
  # score contribution for mu (Eq. 15) should be x_dim X k
  S_beta_ij <- zeta_i  - dA_dmu / params$Ai
  if (!is.null(X)) S_beta_ij <- X * S_beta_ij

  # score contribution for gamma (Eq. 16) should be u_dim X k
  S_gamma_ij <- params$tausq * (zeta_i^2 / 2 - 1 / (2 * params$eta)  - dA_deta / params$Ai)
  if (!is.null(U)) S_gamma_ij <- U * S_gamma_ij
  
  # score contribution for lambda_h (Eq. 17) should be z_dim[h] X k
  scores <- calculate_S_omega_ij(B_mat = params$B_mat,
                                  Ai = params$Ai,
                                  weight_vec = params$weight_vec,
                                  cats = params$cats,
                                  Z0 = Z0,
                                  Z = Z,
                                  lambda0 = params$lambda0,
                                  lambda = params$lambda)
  
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
  
  return(score_contributions)
  
}

#-------------------------------------------------------------------------------
# Hessian function

step_hessian <- function(
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
  
  k <- params$k
  H <- params$H

  # analytic weights
  if (is.null(ai)) ai <- rep(1, k)
  
  d1_mat <- cbind(rep(0,k), dnorm(params$c_mat), rep(0,k))
  d2_mat <- cbind(rep(0,k), params$c_mat * d1_mat[,2:H,drop=FALSE], rep(0,k))
  d3_mat <- cbind(rep(0,k), (params$c_mat^2 - 1) * d1_mat[,2:H,drop=FALSE], rep(0,k))
  d4_mat <- cbind(rep(0,k), (params$c_mat^3 - 3 * params$c_mat) * d1_mat[,2:H,drop=FALSE], rep(0,k))
  
  
  # First derivatives of A w/r/t mu, eta, lambda_h (Eq. 21-23)
  
  dB_dmu <- d1_mat[,2:(H+1L),drop=FALSE] - d1_mat[,1:H,drop=FALSE]
  dA_dmu <- rowSums(dB_dmu * params$lambda_full) / sqrt(params$eta)
  
  dB_deta <- d2_mat[,2:(H+1L),drop=FALSE] - d2_mat[,1:H,drop=FALSE]
  dA_deta <- rowSums(dB_deta * params$lambda_full) / (2 * params$eta)
  
  dA_dlambda <- if (is.null(Z0)) {
    params$B_mat[,-1,drop=FALSE]
  } else {
    params$B_mat
  }
  
  # Second derivatives of A w/r/t mu, eta, lambda_h (Eq. 25-30)
  
  # dA_dmu_dmu is k * 1
  dA_dmu_dmu <- rowSums(dB_deta * params$lambda_full) / params$eta
  
  # dA_dmu_deta is k * 1
  dB_d3 <- d3_mat[,2:(H+1L)] - d3_mat[,1:H]
  dA_dmu_deta <- rowSums(dB_d3 * params$lambda_full) / (2 * params$eta^1.5)
  
  # dA_deta_deta is k * 1
  dB_d4 <- d4_mat[,2:(H+1L)] - d4_mat[,1:H]
  dA_deta_deta <- rowSums(dB_d4 * params$lambda_full) / (4 * params$eta^2)
  
  # dA_dmu_dlambda is k * (H - 1) or k * H if !is.null(Z0)
  # dA_deta_dlambda is k * (H - 1) or k * H if  !is.null(Z0)
  
  if (is.null(Z0)) {
    dA_dmu_dlambda <- dB_dmu[,-1,drop=FALSE] / sqrt(params$eta)
    dA_deta_dlambda <- dB_deta[,-1,drop=FALSE] / (2 * params$eta)
    
    # First derivatives of weights w/r/t lambda_h (Eq. 20)
    dw_dlambda <- model.matrix(~ 0 + params$cats)[,-1,drop=FALSE]
    
  } else {
    dA_dmu_dlambda <- dB_dmu / sqrt(params$eta)
    dA_deta_dlambda <- dB_deta / (2 * params$eta)
    
    # First derivatives of weights w/r/t lambda_h (Eq. 20)
    dw_dlambda <- model.matrix(~ 0 + params$cats)
    
  }
  
  # dA_dlambda_dlambda is all zeros (Eq. 30)
  
  # Second derivatives of weights are all zero (Eq. 24)
  
  
  # First and second derivative of g w/r/t (z_hij * omega_h)
  # (For VHSM this is the same as g(z_hij * omega_h) = lambda_h.)
  
  #----------------------------------------------------------
  # calculate weighted Hessian (Eq. 53 on p. 16)
  # See p. 10 for component matrices
  # Calculate the components by first calculating the pieces inside parentheses,
  # then multiplying by the analytic weights (ai),
  # then taking cross-product with matrix terms (X, U, Z).
  
  # MJ: do i need to multiply by ai somewhere here?
  
  #----------------------------------------------------------
  
  zeta_i <- (yi - params$mu) / params$eta
  
  # H_beta_beta is x_dim * x_dim
  H_beta_beta_right <- dA_dmu^2 / params$Ai^2 - dA_dmu_dmu / params$Ai - 1/params$eta 
  H_beta_beta <- matrix_diag_crossprod(d = ai * H_beta_beta_right, A = X, B = X)

  # H_beta_gamma is x_dim * u_dim
  H_beta_gamma_right <- params$tausq * (dA_dmu * dA_deta / params$Ai^2 - dA_dmu_deta / params$Ai - (yi - params$mu) / params$eta^2)
  H_beta_gamma <- matrix_diag_crossprod(d = ai * H_beta_gamma_right, A = X, B = U)

  # H_gamma_gamma is u_dim * u_dim
  H_gamma_gamma_right <- params$tausq * (zeta_i^2 / 2 - 1 / (2 * params$eta) - dA_deta / params$Ai) + 
                          params$tausq^2 * (dA_deta^2 / params$Ai^2 - dA_deta_deta / params$Ai -  zeta_i^2 / params$eta + 1 / (2 * params$eta^2))
  H_gamma_gamma <- matrix_diag_crossprod(d = ai * H_gamma_gamma_right, A = U, B = U)
  
  # H_beta_omega is x_dim * sum(z_dim) or x_dim * (z0_dim + sum(z_dim)) if !is.null(Z0)
  # H_gamma_omega is u_dim * sum(z_dim) or x_dim * (z0_dim + sum(z_dim)) if !is.null(Z0)
  
  H_omega <- calculate_H_omega(B_mat = params$B_mat,
                               Ai = params$Ai,
                               weight_vec = params$weight_vec,
                               cats = params$cats,
                               eta = params$eta,
                               c_mat = params$c_mat,
                               k= params$k,
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
                               lambda_full = params$lambda_full,
                               dB = list(dmu = dB_dmu, deta = dB_deta))    
  
  H_omega_ <- cbind(t(H_omega$H_beta_omega), t(H_omega$H_gamma_omega), H_omega$H_omega_omega)

   
  # Put the whole Hessian together

  H_matrix <- rbind(
    cbind(H_beta_beta, H_beta_gamma, H_omega$H_beta_omega),
    cbind(t(H_beta_gamma), H_gamma_gamma, H_omega$H_gamma_omega),
    H_omega_
  )
  colnames(H_matrix) <- rownames(H_matrix) <- params$H_names
  
  return(H_matrix)
  
}

