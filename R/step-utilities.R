#-------------------------------------------------------------------------------
# Parse full parameter vector of step function model and set up components

parse_step_params <- function(
    theta,                                      # full parameter vector
    yi,                                         # outcome vector
    sei,                                        # sampling standard errors
    pi = NULL,                                  # one-sided p-values
    ai = NULL,                                  # analytic weight
    beta = NULL,                                # mean parameter coefficients
    gamma = NULL,                               # variance component coefficients
    zeta0 = NULL,                              # selection model coefficients
    zeta = NULL,                               # selection model coefficients
    steps = .025,                               # steps / cut-points
    X = NULL,                                   # mean parameter design matrix
    U = NULL,                                   # variance component design matrix
    Z0 = NULL,                                  # selection model design matrix for highest step
    Z = NULL,
    calc_Ai = FALSE
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
  if (is.null(Z0)) {
    z0_name <- NULL
    z0_dim <- 0L
  } else {
    z0_name <- paste("zeta0", colnames(Z0), sep = "_")
    z0_dim <- ncol(Z0)
  }
  if (is.null(Z)) {
    z_name <- paste0("zeta", 1:length(steps))
    z_dim <- rep(1L, length(steps))
  } else if (is.matrix(Z)) {
    z_name <- paste("zeta", colnames(Z), sep = "_")
    z_dim <- ncol(Z)
  } else {
    z_name_list <- mapply(\(o, zn) paste(o, zn, sep = "_"), 
                          o = paste0("zeta", 1:length(Z)), 
                          zn = lapply(Z, colnames))
    z_name <- unlist(z_name_list)
    z_dim <- sapply(Z, ncol)
  }
  
  H_names <- c(x_name, u_name, z0_name, z_name)
  H <- length(steps) + 1L
  
  # create specific parameters
  if (is.null(beta)) beta <- theta[1:x_dim]
  if (is.null(gamma)) gamma <- theta[x_dim + 1:u_dim]
  if (is.null(zeta0) && z0_dim > 0) zeta0 <- theta[x_dim + u_dim + 1:z0_dim]
  if (is.null(zeta)) zeta <- theta[x_dim + u_dim + z0_dim + 1:sum(z_dim)] 
  if (is.list(Z)) zeta_h <- split(zeta, rep(1:length(steps), z_dim))
  alpha <- c(0, steps, 1)

  # observation-specific parameters
  tausq <- if (is.null(U)) exp(gamma) else exp(as.vector(U %*% gamma))
  eta <- tausq + sei^2
  lambda0 <- if (is.null(Z0)) rep(1, k) else exp(as.vector(Z0 %*% zeta0))  
  lambda <- if (is.null(Z)) {
    sapply(exp(zeta), \(x) rep(x, k))
  } else if (is.matrix(Z)) {
    exp(as.vector(Z %*% zeta))
  } else {
    lambda_h <- mapply(\(z,o) exp(z %*% o), z = Z, o = zeta_h, SIMPLIFY = FALSE)
    lambda <- do.call(cbind, lambda_h)
  } 
  lambda_full <- cbind(lambda0, lambda)
  
  if (!missing(yi)) {
    if (is.null(pi)) pi <- pnorm(yi / sei, lower.tail = FALSE)
    cats <- cut(pi, breaks = c(0, steps, 1), include.lowest = TRUE)
    cats_index <- cbind(1:k, cats)
    weight_vec <- lambda_full[cats_index]
    
    # solve for beta if missing
    if (all(is.na(beta))) {
      wt <- if (is.null(ai)) 1 / (eta * weight_vec) else ai / (eta * weight_vec)
      if (is.null(X)) {
        beta <- stats::weighted.mean(yi, w = wt)
        names(beta) <- "beta"
      } else {
        beta <- stats::lm.wfit(x = X, y = yi, w = wt)$coefficients
        names(beta) <- colnames(X)
      }
    }
  } else {
    weight_vec <- NULL
    cats <- NULL
  }
  
  mu <- if (is.null(X)) beta else as.vector(X %*% beta)
  
  params <- list(k = k,
                 x_dim = x_dim,
                 u_dim = u_dim,
                 z_dim = z_dim,
                 z0_dim = z0_dim,
                 beta = beta,
                 mu = mu, 
                 tausq = tausq,
                 eta = eta,
                 lambda0 = lambda0,
                 lambda = lambda,
                 lambda_full = lambda_full, 
                 weight_vec = weight_vec,
                 H = H,
                 cats = cats,
                 H_names = H_names)
  
  if (calc_Ai) {
    c_mat <- (tcrossprod(sei, -qnorm(steps)) - mu) / sqrt(eta)
    N_mat <- cbind(rep(1,k), pnorm(c_mat), rep(0,k))
    B_mat <- N_mat[,1:H] - N_mat[,2:(H+1L)]
    params$c_mat <- c_mat
    params$B_mat <- B_mat
    params$Ai <- rowSums(B_mat * lambda_full)
  }
  
  return(params)
  
}

#-------------------------------------------------------------------------------
# Calculate score contributions or score vector from components

get_score_contributions <- function(
    S_beta_ij,
    S_gamma_ij,
    S_zeta_ij,
    S_zeta0_ij = NULL,
    ai, 
    contributions
) {
  
  # weighted score contributions (Eq. 52)
  score_contrib <- cbind(S_beta_ij, S_gamma_ij, S_zeta0_ij, S_zeta_ij)
  
  if (is.null(ai)) {
    if (contributions) {
      return(score_contrib)
    } else {
      return(colSums(score_contrib))
    }
  } else {
    if (contributions) {
      return(ai * score_contrib)
    } else {
      return(colSums(ai * score_contrib))
    }
    
  }
  
  
}

#-------------------------------------------------------------------------------
# Calculate crossproducts of the form t(A) %*% diag(d) %*% B
# When A or B can be null

matrix_diag_crossprod <- function(d, A, B) {
  if (missing(d)) {
    if (is.null(A)) {
      if (is.null(B)) {
        1
      } else {
        matrix(colSums(B), nrow = 1L)
      }
    } else {
      if (is.null(B)) {
        matrix(A, ncol = 1L) 
      } else {
        crossprod(A, B)
      }
    }
  } else {
    if (is.null(A)) {
      if (is.null(B)) {
        sum(d)
      } else {
        matrix(colSums(d * B), nrow = 1L)
      }
    } else {
      if (is.null(B)) {
        matrix(colSums(d * A), ncol = 1L) 
      } else {
        crossprod(A, d * B)
      }
    } 
  }
}


#-------------------------------------------------------------------------------
# Calculate S_lambda_ij


calculate_S_zeta_ij <- function(B_mat,
                                  Ai,
                                  weight_vec,
                                  cats,
                                  Z0,
                                  Z,
                                  lambda0,
                                  lambda)
{
  
  
  dA_dlambda <- B_mat
  dw_dlambda <- model.matrix(~ 0 + cats) 
  
  # score contribution for lambda_h (Eq. 17) should be z_dim[h] X k
  if (!is.null(Z0)) {
    S_zeta0_ij <- lambda0 * (dw_dlambda[,1L]  / weight_vec - dA_dlambda[,1L] / Ai) 
    S_zeta0_ij <- Z0 * S_zeta0_ij
  }
  S_zeta_ij <- lambda * (dw_dlambda[,-1L]  / weight_vec - dA_dlambda[,-1L] / Ai) 
  
  if (!is.null(Z)) {
    if (is.matrix(Z)) {
      S_zeta_ij <- Z * S_zeta_ij
    } else {
      S_zeta_ij_list <- apply(S_zeta_ij, 2, identity, simplify = FALSE)
      S_zeta_ij <- mapply(\(z,s) z * s, z = Z, s = S_zeta_ij_list, SIMPLIFY = FALSE)
      S_zeta_ij <- do.call(cbind, S_zeta_ij)
    }
    
  }
  
  scores <- if (!is.null(Z0)) {
    list(S_zeta_ij = S_zeta_ij, S_zeta0_ij = S_zeta0_ij)
  } else {
    list(S_zeta_ij = S_zeta_ij)
  }
  
  return(scores)
  
}



#-------------------------------------------------------------------------------
# Calculate H_zeta

dB_dmu_eta <- function(k,H,c_mat) {
  d1_mat <- cbind(rep(0,k), dnorm(c_mat), rep(0,k))
  d2_mat <- cbind(rep(0,k), c_mat * d1_mat[,2:H,drop=FALSE], rep(0,k))
  
  dB_dmu <- d1_mat[,2:(H+1L),drop=FALSE] - d1_mat[,1:H,drop=FALSE]
  dB_deta <- d2_mat[,2:(H+1L),drop=FALSE] - d2_mat[,1:H,drop=FALSE]
  
  list(dmu = dB_dmu, deta = dB_deta)
}

calculate_H_zeta <- function(
  B_mat,
  Ai,
  weight_vec,
  cats,
  eta,
  c_mat,
  k,
  H,
  tausq,
  X, 
  U,
  Z0,
  Z,
  ai,
  z0_dim,
  z_dim,
  lambda0,
  lambda,
  lambda_full,
  dB = dB_dmu_eta(k, H, c_mat)
){
  
  eta_sqrt <- sqrt(eta)
  
  # First derivatives of A w/r/t mu, eta, lambda_h (Eq. 21-23)
  dA_dmu <- rowSums(dB$dmu * lambda_full) / eta_sqrt
  dA_deta <- rowSums(dB$deta * lambda_full) / (2 * eta)
  
  dA_dlambda <- if (is.null(Z0)) {
    B_mat[,-1,drop=FALSE]
  } else {
    B_mat
  }
  
  if (is.null(Z0)) {
    dA_dmu_dlambda <- dB$dmu[,-1,drop=FALSE] / eta_sqrt
    dA_deta_dlambda <- dB$deta[,-1,drop=FALSE] / (2 * eta)
    
    # First derivatives of weights w/r/t lambda_h (Eq. 20)
    dw_dlambda <- model.matrix(~ 0 + cats)[,-1,drop=FALSE]
    
  } else {
    dA_dmu_dlambda <- dB$dmu / eta_sqrt
    dA_deta_dlambda <- dB$deta / (2 * eta)
    
    # First derivatives of weights w/r/t lambda_h (Eq. 20)
    dw_dlambda <- model.matrix(~ 0 + cats)
    
  }
  
  if (is.null(Z0)) {
    H_beta_zeta_right <- lambda * (dA_dmu * dA_dlambda / Ai^2 - dA_dmu_dlambda/ Ai)
    H_gamma_zeta_right <- tausq * lambda * (dA_deta * dA_dlambda / Ai^2 - dA_deta_dlambda / Ai)
    H_zeta_zeta_right1 <- lambda * (dA_dlambda / Ai)
    H_zeta_zeta_right2 <- lambda * (dw_dlambda / weight_vec)
    Z_full <- if (is.list(Z)) Z else list(Z)
    z_index <- c(0,cumsum(z_dim))
    
  } else {
    H_beta_zeta_right <- lambda_full * ((dA_dmu / Ai^2) * dA_dlambda - dA_dmu_dlambda/ Ai)
    H_gamma_zeta_right <- tausq * lambda_full * ((dA_deta / Ai^2) * dA_dlambda - dA_deta_dlambda / Ai)
    H_zeta_zeta_right1 <- lambda_full * (dA_dlambda / Ai)
    H_zeta_zeta_right2 <- lambda_full * (dw_dlambda / weight_vec)
    Z_full <- if (is.list(Z)) c(list(Z0),Z) else list(Z0, Z)
    z_index <- c(0,cumsum(c(z0_dim, z_dim)))
  }
  
  if (is.null(Z)) {
    if (is.null(Z0)) {
      H_beta_zeta <- matrix_diag_crossprod(A = X, B = ai * H_beta_zeta_right)
      H_gamma_zeta <- matrix_diag_crossprod(A = U, B = ai * H_gamma_zeta_right)
      
      H_zeta_zeta1 <- matrix_diag_crossprod(A = H_zeta_zeta_right1, B = ai * H_zeta_zeta_right1)
      H_zeta_zeta2 <- matrix_diag_crossprod(A = H_zeta_zeta_right2, B = ai * H_zeta_zeta_right2)
      H_zeta_zeta3 <- diag(colSums(ai * H_zeta_zeta_right2) - colSums(ai * H_zeta_zeta_right1), nrow = sum(z_dim))
      H_zeta_zeta <- H_zeta_zeta1 - H_zeta_zeta2 + H_zeta_zeta3
    } else {
      H_beta_zeta_right <- cbind(H_beta_zeta_right[,1] * Z0, H_beta_zeta_right[,-1])
      H_beta_zeta <- matrix_diag_crossprod(A = X, B = ai * H_beta_zeta_right)
      
      H_gamma_zeta_right <- cbind(H_gamma_zeta_right[,1] * Z0, H_gamma_zeta_right[,-1])
      H_gamma_zeta <- matrix_diag_crossprod(A = U, B = ai * H_gamma_zeta_right)
      
      H_zeta0_zeta03 <- crossprod(Z0, ai * (H_zeta_zeta_right2[,1] - H_zeta_zeta_right1[,1]) * Z0)
      H_zeta_zeta3 <- diag(colSums(ai * H_zeta_zeta_right2[,-1,drop=FALSE]) - colSums(ai * H_zeta_zeta_right1[,-1,drop=FALSE]), nrow = sum(z_dim))
      
      H_zeta_zeta_right1 <- cbind(H_zeta_zeta_right1[,1] * Z0, H_zeta_zeta_right1[,-1])
      H_zeta_zeta1 <- matrix_diag_crossprod(A = H_zeta_zeta_right1, B = ai * H_zeta_zeta_right1)
      
      H_zeta_zeta_right2 <- cbind(H_zeta_zeta_right2[,1] * Z0, H_zeta_zeta_right2[,-1])
      H_zeta_zeta2 <- matrix_diag_crossprod(A = H_zeta_zeta_right2, B = ai * H_zeta_zeta_right2)
      
      H_zeta_zeta <- H_zeta_zeta1 - H_zeta_zeta2
      z0_id <- 1:z0_dim
      H_zeta_zeta[z0_id, z0_id] <- H_zeta_zeta[z0_id, z0_id] + H_zeta0_zeta03
      H_zeta_zeta[-z0_id, -z0_id] <- H_zeta_zeta[-z0_id, -z0_id] + H_zeta_zeta3
    }
    
  } else {
    
    H_beta_zeta_right_h <- apply(ai * H_beta_zeta_right, 2, identity, simplify = FALSE)
    H_beta_zeta_h <- mapply(matrix_diag_crossprod, d = H_beta_zeta_right_h, B = Z_full, MoreArgs = list(A = X), SIMPLIFY = FALSE)
    H_beta_zeta <- do.call(cbind, H_beta_zeta_h)
    
    H_gamma_zeta_right_h <- apply(ai * H_gamma_zeta_right, 2, identity, simplify = FALSE)
    H_gamma_zeta_h <- mapply(matrix_diag_crossprod, d = H_gamma_zeta_right_h, B = Z_full, MoreArgs = list(A = U), SIMPLIFY = FALSE)
    H_gamma_zeta <- do.call(cbind, H_gamma_zeta_h)
    
    H_zeta_zeta_right1_h <- apply(H_zeta_zeta_right1, 2, identity, simplify = FALSE)
    H_zeta_zeta1_h <- mapply(\(x,M) x * M, x = H_zeta_zeta_right1_h, M = Z_full, SIMPLIFY = FALSE)
    H_zeta_zeta1_R <- do.call(cbind, H_zeta_zeta1_h)
    H_zeta_zeta1 <- crossprod(H_zeta_zeta1_R, ai * H_zeta_zeta1_R)
    
    H_zeta_zeta_right2_h <- apply(H_zeta_zeta_right2, 2, identity, simplify = FALSE)
    H_zeta_zeta2_h <- mapply(\(x,M) x * M, x = H_zeta_zeta_right2_h, M = Z_full, SIMPLIFY = FALSE)
    H_zeta_zeta2_R <- do.call(cbind, H_zeta_zeta2_h)
    H_zeta_zeta2 <- crossprod(H_zeta_zeta2_R, ai * H_zeta_zeta2_R)
    
    H_zeta_zeta3_h <- mapply(\(x,y,M) crossprod(M, ai * (y - x) * M), 
                               x = H_zeta_zeta_right1_h, 
                               y = H_zeta_zeta_right2_h, 
                               M = Z_full, SIMPLIFY = FALSE)
    H_zeta_zeta <- H_zeta_zeta1 - H_zeta_zeta2
    for (i in 1:length(H_zeta_zeta3_h)) {
      id <- (z_index[i] + 1):(z_index[i+1])
      H_zeta_zeta[id, id] <- H_zeta_zeta[id, id] + H_zeta_zeta3_h[[i]]
    }
  }
  
  
  res <- list(H_beta_zeta = H_beta_zeta, 
              H_gamma_zeta = H_gamma_zeta, 
              H_zeta_zeta = H_zeta_zeta)
  
  return(res)
  
}

