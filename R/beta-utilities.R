# Parse full parameter vector of beta function model and set up components

parse_beta_params <- function(
    theta,                                     # full parameter vector
    yi,                                        # outcome vector
    sei,                                       # sampling standard errors
    pi = NULL,                                 # one-sided p-values
    beta = NULL,                               # mean parameter coefficients
    gamma = NULL,                              # variance component coefficients
    zeta = NULL,                               # selection model coefficients
    alpha = c(.025,.975),                      # p-value truncation points
    X = NULL,                                  # mean parameter design matrix
    U = NULL,                                  # variance component design matrix
    calc_Ai = FALSE,                           # whether to calculate Ai 
    calc_Ai_deriv = FALSE                      # whether to calculate first derivatives of Ai
) {
  
  # number of observations
  k <- length(sei)
  
  # parameter dimensions and names
  if (is.null(X)) {
    x_name <- "beta"
    x_dim <- 1L
  } else {
    x_name <- paste("beta", colnames(X), sep = "_")
    x_dim <- ncol(X)
  }
  if (is.null(U)) {
    u_name <- "gamma"
    u_dim <- 1L
  } else {
    u_name <- paste("gamma", colnames(U), sep = "_")
    u_dim <- ncol(U)
  }  
  
  z_name <- c("zeta1", "zeta2")
  
  H_names <- c(x_name, u_name, z_name)

  # create specific parameters
  if (is.null(beta)) beta <- theta[1:x_dim]
  if (is.null(gamma)) gamma <- theta[x_dim + 1:u_dim]
  if (is.null(zeta)) zeta <- theta[x_dim + u_dim + 1:2]

  # observation-specific parameters
  mu <- if (is.null(X)) beta else as.vector(X %*% beta)
  tausq <- if (is.null(U)) exp(gamma) else exp(as.vector(U %*% gamma))
  eta <- tausq + sei^2
  
  lambda <- exp(zeta) - 1
  
  alpha_lambda <- (alpha ^ lambda[1]) * ((1 - alpha) ^ lambda[2])

  # calculate weights 
  if (!missing(yi)) {
    if (is.null(pi)) pi <- pnorm(yi / sei, lower.tail = FALSE)
    pi_tilde <- pmin(alpha[2], pmax(alpha[1], pi))
    weight_vec <- (pi_tilde ^ lambda[1]) * ((1 - pi_tilde) ^ lambda[2])
  } else {
    pi_tilde <- NULL
    weight_vec <- NULL
  }

  params <- list(k = k,
                 x_dim = x_dim,
                 u_dim = u_dim,
                 mu = mu, 
                 tausq = tausq,
                 eta = eta,
                 alpha = alpha,
                 zeta = zeta,
                 lambda = lambda,
                 alpha_lambda = alpha_lambda,
                 pi_tilde = pi_tilde,
                 weight_vec = weight_vec,
                 H_names = H_names)
  
  if (calc_Ai) {
    
    c_mat <- (tcrossprod(sei, -qnorm(alpha)) - mu) / sqrt(eta)
    B_0ij <- pnorm(c_mat[,1], lower.tail = FALSE)
    B_2ij <- pnorm(c_mat[,2])
    
    params$c_1ij <- c_mat[,1]
    params$c_2ij <- c_mat[,2]
    params$B_0ij <- B_0ij
    params$B_2ij <- B_2ij
    
    
    E_Y_1 <- E_Y_f_vec(
      f_exp = "1", 
      sei = sei, mu = mu, eta = eta, 
      lambda = lambda, alpha = alpha
    )
    
    params$Ai <- alpha_lambda[1] * B_0ij + E_Y_1 + alpha_lambda[2] * B_2ij

  }
  
  if (calc_Ai_deriv) {
    
    # calculate derivatives of the weights (Eq. 32-33)
    params$g_dot <- exp(zeta)
    params$dw_dlambda1 <- weight_vec * log(pi_tilde)
    params$dw_dlambda2 <- weight_vec * log(1 - pi_tilde)
    
    # calculate derivatives of A_ij w.r.t. beta, eta, lambda (Eq. 34-37)
    
    EY_A_mu <- E_Y_f_vec(
      f_exp = "(Y - mu) / sqrt(eta)", 
      sei = sei, mu = mu, eta = eta,
      lambda = lambda, alpha = alpha
    )
    
    params$dA_dmu <- (alpha_lambda[1] * dnorm(c_mat[,1]) + EY_A_mu  - alpha_lambda[2] * dnorm(c_mat[,2])) / sqrt(eta)

    EY_A_eta <- E_Y_f_vec(
      f_exp = "((Y - mu)^2 / eta - 1)", 
      sei = sei, mu = mu, eta = eta,
      lambda = lambda, alpha = alpha
    )

    params$dA_deta <- 
      (alpha_lambda[1] * c_mat[,1] * dnorm(c_mat[,1]) 
       + EY_A_eta 
       - alpha_lambda[2] * c_mat[,2] * dnorm(c_mat[,2])) / (2 * eta)

    
    EY_A_lambda_1 <- E_Y_f_vec(
      f_exp = "pnorm(-Y/sei, log.p = TRUE)", 
      sei = sei, mu = mu, eta = eta,
      lambda = lambda, alpha = alpha
    )
    
    params$dA_dlambda1 <- 
      log(alpha[1]) * alpha_lambda[1] * B_0ij + 
      EY_A_lambda_1 + 
      log(alpha[2]) * alpha_lambda[2] * B_2ij

    EY_A_lambda_2 <- E_Y_f_vec(
      f_exp = "pnorm(Y/sei, log.p = TRUE)", 
      sei = sei, mu = mu, eta = eta,
      lambda = lambda, alpha = alpha
    )
    
    params$dA_dlambda2 <- 
      log(1 - alpha[1]) * alpha_lambda[1] * B_0ij + 
      EY_A_lambda_2 + 
      log(1 - alpha[2]) * alpha_lambda[2] * B_2ij

  }
  
  return(params)
  
}


# integration ------------------------------------------------------------------

E_Y_f <- function(f_exp, 
                  sei, 
                  mu, 
                  eta, 
                  lambda, 
                  alpha) {
  
  weightfun <- "(pnorm(-Y/sei)^lambda[1]) * (pnorm(Y/sei)^lambda[2])"
  normdens <- "dnorm((Y - mu) / sqrt(eta)) / sqrt(eta)"
  integrand_exp <- paste(f_exp, weightfun, normdens, sep = " * ")
  integrand <- function(Y) Y
  body(integrand) <- str2lang(integrand_exp)
  
  # added the upper and lower bounds here 
  bounds <- sei * qnorm(alpha, lower.tail = FALSE)
  
  res <- integrate(integrand, lower = bounds[2], upper = bounds[1])
  res$value  
  
}

E_Y_f_vec <- function(f_exp, sei, mu, eta, lambda, alpha) {
  mapply(
    E_Y_f, 
    sei = sei, mu = mu, eta = eta,
    MoreArgs = list(
      f_exp = f_exp, 
      lambda = lambda, alpha = alpha
    )
  )
}
