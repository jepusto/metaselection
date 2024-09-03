
#-------------------------------------------------------------------------------
# Log likelihood function

beta_loglik <- function(
    theta,                                     # full parameter vector
    yi,                                        # outcome vector
    sei,                                       # sampling standard errors
    pi = pnorm(yi / sei, lower.tail = FALSE),  # one-sided p-values
    ai = NULL,                                 # analytic weight
    beta = NULL,                               # mean parameter coefficients
    gamma = NULL,                              # variance component coefficients
    zeta = NULL,
    steps = c(.025,.975),                      # p-value truncation points
    X = NULL,                                  # mean parameter design matrix
    U = NULL                                   # variance component design matrix
) {
  
  params <- parse_beta_params(
    theta = theta,
    yi = yi,
    sei = sei,
    pi = pi,
    beta = beta,
    gamma = gamma,
    zeta = zeta,
    alpha = steps,
    X = X,
    U = U,
    calc_Ai = TRUE
  )

  # calculate log likelihood (Eq. 9-10)
  # likelihood contributions
  log_lik_i <- log(params$weight_vec) - (yi - params$mu)^2 / (2 * params$eta) - log(params$eta) / 2 - log(params$Ai)
  
  # weighted log likelihood (Eq. 51)
  log_lik <- if (is.null(ai)) sum(log_lik_i) else sum(ai * log_lik_i)
    
    
  return(log_lik)

}

#-------------------------------------------------------------------------------
# Score function

beta_score <- function(
    theta,                                     # full parameter vector
    yi,                                        # outcome vector
    sei,                                       # sampling standard errors
    pi = pnorm(yi / sei, lower.tail = FALSE),  # one-sided p-values
    ai = NULL,                                 # analytic weight
    beta = NULL,                               # mean parameter coefficients
    gamma = NULL,                              # variance component coefficients
    zeta = NULL,                              # selection model coefficients
    steps = c(.025,.975),                      # p-value truncation points
    X = NULL,                                  # mean parameter design matrix
    U = NULL,                                  # variance component design matrix
    contributions = FALSE                      # whether to return matrix of score contributions
) {

  params <- parse_beta_params(
    theta = theta,
    yi = yi,
    sei = sei,
    pi = pi,
    beta = beta,
    gamma = gamma,
    zeta = zeta,
    alpha = steps,
    X = X,
    U = U,
    calc_Ai = TRUE,
    calc_Ai_deriv = TRUE
  )

  # calculate score components (Eq. 15-17)
  
  zeta_i <- (yi - params$mu) / params$eta

  # score contribution for mu (Eq. 15) should be x_dim X k
  S_beta_ij <- zeta_i  - params$dA_dmu / params$Ai
  if (!is.null(X)) S_beta_ij <- X * S_beta_ij
  
  # score contribution for gamma (Eq. 16) should be u_dim X k
  S_gamma_ij <- params$tausq * (zeta_i^2 / 2 - 1 / (2 * params$eta)  - params$dA_deta / params$Ai)
  if (!is.null(U)) S_gamma_ij <- U * S_gamma_ij
  
  # score contribution for lambda_1
  # is this correct ?? 
  S_lambda_ij_1 <- params$g_dot[1] * (params$dw_dlambda1  / params$weight_vec - params$dA_dlambda1 / params$Ai)
  S_lambda_ij_2 <- params$g_dot[2] * (params$dw_dlambda2  / params$weight_vec - params$dA_dlambda2 / params$Ai)
  
  score_contributions <- get_score_contributions(
    S_beta_ij = S_beta_ij,
    S_gamma_ij = S_gamma_ij,
    S_zeta_ij = cbind(S_lambda_ij_1, S_lambda_ij_2),
    ai = ai,
    contributions = contributions
  )

  return(score_contributions)

}

#-------------------------------------------------------------------------------
# Hessian function

beta_hessian <- function(
    theta,                                     # full parameter vector
    yi,                                        # outcome vector
    sei,                                       # sampling standard errors
    pi = pnorm(yi / sei, lower.tail = FALSE),  # one-sided p-values
    ai = NULL,                                 # analytic weight
    beta = NULL,                               # mean parameter coefficients
    gamma = NULL,                              # variance component coefficients
    zeta = NULL,                              # selection model coefficients
    steps = c(.025,.975),                      # p-value truncation points
    X = NULL,                                  # mean parameter design matrix
    U = NULL                                   # variance component design matrix
) {
  
  params <- parse_beta_params(
    theta = theta,
    yi = yi,
    sei = sei,
    pi = pi,
    beta = beta,
    gamma = gamma,
    zeta = zeta,
    alpha = steps,
    X = X,
    U = U,
    calc_Ai = TRUE,
    calc_Ai_deriv = TRUE
  )

  k <- params$k

  # analytic weights
  if (is.null(ai)) ai <- rep(1, k)
  
  # calculate second derivatives of the weights (Eq. 38-40)
  dw_dlambda1_dlambda1 <- params$weight_vec * (log(params$pi_tilde))^2
  dw_dlambda2_dlambda2 <- params$weight_vec * (log(1 - params$pi_tilde))^2
  dw_dlambda1_dlambda2 <- params$weight_vec * log(params$pi_tilde) * log(1 - params$pi_tilde)
  
  # calculate second derivatives of A_ij w.r.t. beta, eta, lambda (Eq. 41-50)
  
  EY_A_mu_mu <-  E_Y_f_vec(
    f_exp = "((Y - mu) ^ 2 / eta - 1)", 
    sei = sei, mu = params$mu, eta = params$eta,
    lambda = params$lambda, alpha = params$alpha
  )

  d_1ij <- dnorm(params$c_1ij)
  d_2ij <- dnorm(params$c_2ij)
  
  dA_dmu_dmu <- (params$alpha_lambda[1] * params$c_1ij * d_1ij
                 + EY_A_mu_mu
                 - params$alpha_lambda[2] * params$c_2ij * d_2ij) / params$eta
  
  EY_A_mu_eta <- E_Y_f_vec(
    f_exp = " ((Y - mu) * ((Y - mu)^2 - 3 * eta) / eta^(3/2))", 
    sei = sei, mu = params$mu, eta = params$eta,
    lambda = params$lambda, alpha = params$alpha
  )
  
  dA_dmu_deta <- (params$alpha_lambda[1] * (params$c_1ij^2 - 1) * d_1ij
                  + EY_A_mu_eta
                  - params$alpha_lambda[2] * (params$c_2ij^2 - 1) * d_2ij) / (2 * params$eta^1.5)
  
  EY_A_mu_lambda1 <- E_Y_f_vec(
    f_exp = "((Y - mu) / sqrt(eta)) * pnorm(-Y/sei, log.p = TRUE)", 
    sei = sei, mu = params$mu, eta = params$eta,
    lambda = params$lambda, alpha = params$alpha
  )
  
  log_alpha_1 <- log(params$alpha[1])
  log_alpha_2 <- log(params$alpha[2])
  log_c_alpha_1 <- log(1 - params$alpha[1])
  log_c_alpha_2 <- log(1 - params$alpha[2])
  
  dA_dmu_dlambda1 <- (log_alpha_1 * params$alpha_lambda[1] * d_1ij
                      + EY_A_mu_lambda1
                      - log_alpha_2 * params$alpha_lambda[2] * d_2ij) / sqrt(params$eta)
  
  EY_A_mu_lambda2 <- E_Y_f_vec(
    f_exp = "((Y - mu) / sqrt(eta)) * pnorm(Y/sei, log.p = TRUE)", 
    sei = sei, mu = params$mu, eta = params$eta,
    lambda = params$lambda, alpha = params$alpha
  )

  dA_dmu_dlambda2 <- (log_c_alpha_1 * params$alpha_lambda[1] * d_1ij
                      + EY_A_mu_lambda2
                      - log_c_alpha_2 * params$alpha_lambda[2] * d_2ij) / sqrt(params$eta)
  
  EY_A_eta_eta <- E_Y_f_vec(
    f_exp = "( (Y - mu)^4 / eta^2 - 6 * (Y - mu)^2 / eta + 3)", 
    sei = sei, mu = params$mu, eta = params$eta,
    lambda = params$lambda, alpha = params$alpha
  )
 
  dA_deta_deta <- (params$alpha_lambda[1] * (params$c_1ij^3 - 3 * params$c_1ij) * d_1ij
                   + EY_A_eta_eta
                   - params$alpha_lambda[2] * (params$c_2ij^3 - 3 * params$c_2ij) * d_2ij) / (4 * params$eta^2)
  
  EY_A_eta_lambda1 <- 
    E_Y_f_vec(
      f_exp = "pnorm(-Y / sei, log.p = TRUE) * ((Y - mu)^2 / eta - 1)", 
      sei = sei, mu = params$mu, eta = params$eta,
      lambda = params$lambda, alpha = params$alpha
    )

  dA_deta_dlambda1 <- (log_alpha_1 * params$alpha_lambda[1] * params$c_1ij * d_1ij
                       + EY_A_eta_lambda1
                       - log_alpha_2 * params$alpha_lambda[2] * params$c_2ij * d_2ij) / (2 * params$eta)
  
  EY_A_eta_lambda2 <- E_Y_f_vec(
    f_exp = "pnorm(Y / sei, log.p = TRUE) * ((Y - mu)^2 / eta - 1)", 
    sei = sei, mu = params$mu, eta = params$eta,
    lambda = params$lambda, alpha = params$alpha
  )

  dA_deta_dlambda2 <- (log_c_alpha_1 * params$alpha_lambda[1] * params$c_1ij * d_1ij
                       + EY_A_eta_lambda2
                       - log_c_alpha_2 * params$alpha_lambda[2] * params$c_2ij * d_2ij) / (2 * params$eta)
  
  EY_A_lambda1_lambda1 <- E_Y_f_vec(
    f_exp = "pnorm(-Y / sei, log.p = TRUE)^2", 
    sei = sei, mu = params$mu, eta = params$eta,
    lambda = params$lambda, alpha = params$alpha
  )
  
  dA_dlambda1_dlambda1 <- 
    log_alpha_1^2 * params$alpha_lambda[1] * params$B_0ij + 
    EY_A_lambda1_lambda1 + 
    log_alpha_2^2 * params$alpha_lambda[2] * params$B_2ij
  
  EY_A_lambda2_lambda2 <- E_Y_f_vec(
    f_exp = "pnorm(Y / sei, log.p = TRUE)^2", 
    sei = sei, mu = params$mu, eta = params$eta,
    lambda = params$lambda, alpha = params$alpha
  )
  
  dA_dlambda2_dlambda2 <- 
    log_c_alpha_1^2 * params$alpha_lambda[1] * params$B_0ij + 
    EY_A_lambda2_lambda2 + 
    log_c_alpha_2^2 * params$alpha_lambda[2] * params$B_2ij
  
  EY_A_lambda1_lambda2 <- E_Y_f_vec(
    f_exp = "pnorm(-Y / sei, log.p = TRUE) * pnorm(Y / sei, log.p = TRUE)", 
    sei = sei, mu = params$mu, eta = params$eta,
    lambda = params$lambda, alpha = params$alpha
  )

  dA_dlambda1_dlambda2 <- 
    log_alpha_1 * log_c_alpha_1 * params$alpha_lambda[1] * params$B_0ij + 
    EY_A_lambda1_lambda2 + 
    log_alpha_2 * log_c_alpha_2 * params$alpha_lambda[2] * params$B_2ij
  
  
  # calculate Hessian components (p. 10)

  #----------------------------------------------------------
  
  # same stuff from step 
  
  zeta_i <- (yi - params$mu) / params$eta
  
  # H_beta_beta is x_dim * x_dim
  H_beta_beta_right <- params$dA_dmu^2 / params$Ai^2 - dA_dmu_dmu / params$Ai - 1 / params$eta 
  H_beta_beta <- matrix_diag_crossprod(d = ai * H_beta_beta_right, A = X, B = X)
  
  # H_beta_gamma is x_dim * u_dim
  H_beta_gamma_right <- params$tausq * (params$dA_dmu * params$dA_deta / params$Ai^2 - dA_dmu_deta / params$Ai - zeta_i / params$eta)
  H_beta_gamma <- matrix_diag_crossprod(d = ai * H_beta_gamma_right, A = X, B = U)
  
  # H_gamma_gamma is u_dim * u_dim
  H_gamma_gamma_right <- params$tausq * (zeta_i^2 / 2 - 1 / (2 * params$eta) - params$dA_deta / params$Ai) + 
    params$tausq^2 * (params$dA_deta^2 / params$Ai^2 - dA_deta_deta / params$Ai -  zeta_i^2 / params$eta + 1 / (2 * params$eta^2))
  H_gamma_gamma <- matrix_diag_crossprod(d = ai * H_gamma_gamma_right, A = U, B = U)
  
  
  H_beta_zeta_1_right <-  as.matrix(params$g_dot[1] * (params$dA_dmu * params$dA_dlambda1 / params$Ai^2 - dA_dmu_dlambda1/ params$Ai))
  H_beta_zeta_1 <- matrix_diag_crossprod(A = X, B = ai * H_beta_zeta_1_right)
  
  H_beta_zeta_2_right <-  as.matrix(params$g_dot[2] * (params$dA_dmu * params$dA_dlambda2 / params$Ai^2 - dA_dmu_dlambda2/ params$Ai))
  H_beta_zeta_2 <- matrix_diag_crossprod(A = X, B = ai * H_beta_zeta_2_right)
  
  H_gamma_zeta_1_right <- as.matrix(params$tausq * params$g_dot[1] * (params$dA_deta * params$dA_dlambda1 / params$Ai^2 - dA_deta_dlambda1 / params$Ai))
  H_gamma_zeta_1 <- matrix_diag_crossprod(A = U, B = ai * H_gamma_zeta_1_right)
  
  H_gamma_zeta_2_right <- as.matrix(params$tausq * params$g_dot[2] * (params$dA_deta * params$dA_dlambda2 / params$Ai^2 - dA_deta_dlambda2 / params$Ai))
  H_gamma_zeta_2 <- matrix_diag_crossprod(A = U, B = ai * H_gamma_zeta_2_right)
  
  
  # H zeta zeta -----------------------------------------------------------

  H_zeta_1_zeta_1_right <- 
    params$g_dot[1] * params$g_dot[1]  * 
    (dw_dlambda1_dlambda1 / params$weight_vec - params$dw_dlambda1^2 / params$weight_vec^2
     + params$dA_dlambda1^2 / params$Ai^2 - dA_dlambda1_dlambda1 / params$Ai) + 
    params$g_dot[1] * (params$dw_dlambda1 / params$weight_vec - params$dA_dlambda1 / params$Ai)
  H_zeta_1_zeta_1 <- sum(ai * H_zeta_1_zeta_1_right)
  
  H_zeta_2_zeta_2_right <- 
    params$g_dot[2] * params$g_dot[2]  * 
    (dw_dlambda2_dlambda2 / params$weight_vec - params$dw_dlambda2^2 / params$weight_vec^2
     + params$dA_dlambda2^2 / params$Ai^2 - dA_dlambda2_dlambda2 / params$Ai) + 
    params$g_dot[2] * (params$dw_dlambda2 / params$weight_vec - params$dA_dlambda2 / params$Ai)
  H_zeta_2_zeta_2 <- sum(ai * H_zeta_2_zeta_2_right)
  
  H_zeta_1_zeta_2_right <- 
    params$g_dot[1] * params$g_dot[2] * 
    (dw_dlambda1_dlambda2 / params$weight_vec - params$dw_dlambda1 * params$dw_dlambda2 / params$weight_vec^2
     + params$dA_dlambda1 * params$dA_dlambda2 / params$Ai^2 - dA_dlambda1_dlambda2 / params$Ai)
  H_zeta_1_zeta_2 <- sum(ai * H_zeta_1_zeta_2_right)
  
  # assemble the full Hessian matrix
  
  H_matrix <- rbind(
    cbind(H_beta_beta, H_beta_gamma, H_beta_zeta_1, H_beta_zeta_2),
    cbind(t(H_beta_gamma), H_gamma_gamma, H_gamma_zeta_1, H_gamma_zeta_2),
    c(H_beta_zeta_1, H_gamma_zeta_1, H_zeta_1_zeta_1, H_zeta_1_zeta_2),
    c(H_beta_zeta_2, H_gamma_zeta_2, H_zeta_1_zeta_2, H_zeta_2_zeta_2)
  )
  
  colnames(H_matrix) <- rownames(H_matrix) <- params$H_names

  
  return(H_matrix)

}

