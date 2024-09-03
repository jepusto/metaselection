# -----------------------#
# Data Generating Model  #
# -----------------------#

# simulate empirical distribution of sample size and number of effect size
# n_ES_empirical() is a functional (a function that returns a function). 
# You give it a dataset with primary study sample sizes and numbers of effect sizes per study. 
# It returns a function that generates random samples from the dataset—random study characteristics.
#' 
#' @export

n_ES_empirical <- function(dat) {
  d <- dat
  function(m) d[sample(NROW(dat), size = m, replace = TRUE), ]
  
}

#' 
#' @export

n_ES_param <- function(mean_N, mean_ES, min_N = 20L) {
  function(m) {
    data.frame(
      n = min_N + rpois(m, lambda = mean_N - min_N),
      n_ES = 1L + rpois(m, lambda = mean_ES - 1L)
    )
  }
}

# Reparameterization

r_corr <- function(m, cor_mu, cor_sd) {
  if (cor_sd > 0 & cor_mu > 0) {
    alpha <- cor_mu * ((cor_mu * (1 - cor_mu) / cor_sd^2 ) - 1)
    beta <-  (1 - cor_mu) * ((cor_mu * (1 - cor_mu) / cor_sd^2 ) - 1)
    rbeta(n = m, shape1 = alpha, shape2 = beta)
  } else {
    rep(cor_mu, m)
  }
  
}

# generate primary data for study j ---------------------------------------
 

# Generate subject and primary study level data
## generates i SMD estimates for j participants each primary study k
## based on different outcome measures taken on the sample
## correlated due to sampling error and variation in true effects
## n is per-group sample size ; n_es is number effects per individual n
## r_j is correlation among outcomes, used for the sigma_matrix

r_study <- function(delta_j, # 
                    r_j, # correlation between outcomes
                    n, # sample size
                    n_ES) { # number of es per individual 
  
  # create treatment and control groups
  # divide by 2 and round for uneven sample sizes
  n_tx <- round((n / 2), 0)
  
  # Create Sigma matrix
  Sigma_mat <- r_j + diag(rep(1 - r_j, n_ES), nrow = n_ES)
  Sigma_mat <- as.matrix(Sigma_mat)
  
  # Generate treatment/control deltas and necessary calculations
  tx_j <- mvtnorm::rmvnorm(n_tx, mean = rep(delta_j, length.out = n_ES), sigma = Sigma_mat)
  ybar_tx <- colMeans(tx_j)
  var_tx <- diag(cov(tx_j))
  
  ctl_j <- mvtnorm::rmvnorm(n - n_tx, mean = rep(0, n_ES), sigma = Sigma_mat)
  ybar_ctl <- colMeans(ctl_j)
  var_ctl <- diag(cov(ctl_j))
  
  # Calculate other necessary stats
  s_pooled <- ((n_tx - 1) * var_tx + (n - n_tx - 1) * var_ctl) / (n - 2)
  J <- 1 - (3 / (4 * (n - 2) - 1))
  d <- (ybar_tx - ybar_ctl) / sqrt(s_pooled) * J
  var_d <- J^2 * 4 / n + d^2 / (2 * (n - 2))
  sd_d <- sqrt(var_d)
  Va <- J^2 * 4 / n
  sda <- sqrt(Va)
  df <- n - 2
  t_i <- d / sqrt(Va)
  p_onesided <- pt(t_i, df = df, lower.tail = FALSE) 
  
  # Keep only necessary values
  data.frame(n, n_ES, d, var_d, sd_d, Va, sda, t_i, p_onesided, esid = 1:n_ES)
}

# censoring functions ---------------------------------------------
#' @export

step_fun <- function(cut_vals = .025, weights = 1) {
  
  if (length(cut_vals) != length(weights)) stop("cut_vals and weights must be the same length, doofus!")
  
  wt_vec <- c(1, weights) / max(c(1, weights))
  
  cut_vals_full <- c(0, cut_vals, 1)
  
  function(pvals) {
    
    pval_buckets <- cut(pvals, cut_vals_full, labels = FALSE, include.lowest = TRUE)
    wt_vec[pval_buckets]
    
  }
}

#' @export

beta_wts_fun <- function(delta_1 = 1, delta_2 = 1,
                         trunc_1 = 0.025, trunc_2 = 0.975) {

  if (delta_1 + delta_2 > 2) {
    max_p <- (delta_1 - 1) / (delta_1 + delta_2 - 2)
    if (is.nan(max_p)) max_p <- 0.5
    max_p <- max(min(max_p, trunc_2), trunc_1)
  } else {
    max_p <- if (delta_1 > delta_2) trunc_2 else trunc_1
  }
  
  max_val <- max_p^(delta_1 - 1) * (1 - max_p)^(delta_2 - 1)
  
  function(pvals) {
    
    # truncating 
    pvals <- ifelse(pvals < trunc_1, trunc_1, pvals)
    pvals <- ifelse(pvals > trunc_2, trunc_2, pvals)
    
    pvals^(delta_1 - 1) * (1 - pvals)^(delta_2 - 1) / max_val
    
  }
  
}


# generate meta analytic data ---------------------------------------------

# Generate each meta-analytic dataset 
# Sample k studies from r_study function
# establish multi-variate normal distribution of n_ES per study
# var/covariance structure created using mean correlation (cor_mu) and sd (cor_sd)
# What was n_ES_gen (number generated?)


#' @title Generate meta-analytic data 
#' 
#' @description Generate meta-analytic correlated or correlated and hierarchical effects data with options to simulate selective outcome reporting
#' 
#' @param mean_smd number indicating the true mean effect size 
#' @param tau number characterizing between-study heterogeneity in effects
#' @param omega number characterizing within-study heterogeneity in effects
#' @param m number of studies in the simulated meta-analysis
#' @param cor_mu number indicating the average correlation between outcomes
#' @param cor_sd number indicating standard deviation of correlation between outcomes
#' @param censor_fun a function used to censor effects. This package provides functions `step_fun()` and `beta_wts_fun()`... ???
#' @param n_ES_sim a function used to simulate the distribution of primary study sample sizes and the number of effect sizes per study
#' @param m_multiplier number indicating a multiplier for buffer for the number of studies 
#' @param id_start integer indicating the starting value for id
#' @param paste_ids logical with \code{TRUE} (the default) indicating that the study id and effect size id should be pasted together
#' @param include_sel_prob logical with \code{TRUE} indicating 
#' 
#' @returns A data frame containing the simulated meta-analytic dataset.
#' 
#' 
#' @examples
#' 
#' example_dat <- r_meta(
#'   mean_smd = 0, 
#'   tau = .1, omega = .01,
#'   m = 50, 
#'   cor_mu = .4, cor_sd = 0.001, 
#'   censor_fun = step_fun(cut_vals = .025, weights = 0.4), 
#'   n_ES_sim = n_ES_param(40, 3)
#' )
#' 
#' @export



r_meta <- function(
    mean_smd, # true mean of effect sizes
    tau, # between-study heterogeneity
    omega, # within-study heterogeneity
    m, # number of studies in each meta analysis
    cor_mu, # average correlation between outcomes
    cor_sd, # sd correlation between outcomes
    censor_fun, # censoring function
    n_ES_sim, # distribution of sample sizes and number of outcomes
    m_multiplier = 2,
    id_start = 0L,
    paste_ids = TRUE,
    include_sel_prob = FALSE
) {
  
  # Calculate buffer for number of studies to generate
  m_star <- round(m * m_multiplier)
  
  # Sample correlation values for matrix
  r_j <- r_corr(m = m_star, cor_mu, cor_sd) 

  # Simulate delta_j for each individual within each study. 
  delta_j <- rnorm(m_star, mean = mean_smd, sd = tau)
  
  # Sample a sample-size and n_ES combo 
  n_select <- as.data.frame(n_ES_sim(m_star))
  
  # generate study design parameters
  study_parms <- data.frame(delta_j = delta_j, r_j = r_j, n = n_select$n, n_ES = n_select$n_ES)
  
  # add within-study heterogeneity 
  delta_ij <- purrr::pmap(study_parms, .f = \(delta_j, n_ES, ...) rnorm(n_ES, mean = delta_j, sd = omega))
  study_parms$delta_j <- delta_ij
  
  # generate the studies
  studies <- purrr::pmap_dfr(study_parms, .f = r_study, .id = "studyid")
  studies$studyid <- id_start + as.integer(studies$studyid)

  # apply censoring process
  K_star <- nrow(studies)
  p_sel <- censor_fun(studies$p_onesided)
  if (include_sel_prob) studies$selection_prob <- p_sel
  observed <- as.logical(stats::rbinom(K_star, size = 1L, prob = p_sel))
  studies <- studies[observed,,drop=FALSE]
  
  # Count how many unique studies are kept
  n_studies <- length(unique(studies$studyid))
  
  # recurse to get desired number of studies per meta-analytic dataset
  if (n_studies < m) {
    
    more_studies <- r_meta(
      mean_smd = mean_smd, 
      tau = tau, omega = omega, 
      m = m - n_studies, 
      cor_mu = cor_mu, cor_sd = cor_sd, 
      censor_fun = censor_fun, 
      n_ES_sim = n_ES_sim, 
      m_multiplier = m_multiplier,
      id_start = id_start + m_star,
      paste_ids = FALSE
    )
    
    if (nrow(studies) > 0L) {
      studies <- rbind(studies, more_studies)
    } else {
      studies <- more_studies
    }
  }
  
  # keep only desired number of studies
  if (paste_ids) studies$esid <- paste(studies$studyid, studies$esid, sep = "-")
  studies$studyid <- factor(studies$studyid)
  studies <- droplevels(studies[studies$studyid %in% (levels(studies$studyid)[1:m]),])
  
  return(studies)
  
}


# Map over multiple parameter values to generate meta-regression

r_meta_categories <- function(
  mean_smd, # true mean of effect sizes
  tau, # between-study heterogeneity
  omega, # within-study heterogeneity
  m, # number of studies in each meta analysis
  cor_mu, # average correlation between outcomes
  cor_sd, # sd correlation between outcomes
  censor_fun, # censoring function 
  n_ES_sim, # distribution of sample sizes and number of outcomes
  m_multiplier = 2,
  id_start = 0L,
  paste_ids = TRUE
) {   
  
  max_param_length <- max(lengths(list(mean_smd, tau, omega, m, cor_mu, cor_sd, censor_fun)))
  
  params <- data.frame(
    mean_smd = rep(mean_smd, length.out = max_param_length),
    tau = rep(tau, length.out = max_param_length),
    omega = rep(omega, length.out = max_param_length),
    m = rep(m, length.out = max_param_length),
    cor_mu = rep(cor_mu, length.out = max_param_length),
    cor_sd = rep(cor_sd, length.out = max_param_length)
  )
  
  if (is.list(censor_fun)) {
    params$censor_fun <- rep(censor_fun, length.out = max_param_length)
  } else {
    params$censor_fun <- rep(list(censor_fun), length.out = max_param_length)
  }
  
  studies <- purrr::pmap_dfr(
    params, 
    r_meta, 
    n_ES_sim = n_ES_sim, 
    m_multiplier = m_multiplier,
    paste_ids = FALSE,
    .id = "X"
  )

  studies$studyid <- factor(paste(studies$X, studies$studyid, sep = "-"))
  if (paste_ids) studies$esid <- paste(studies$studyid, studies$esid, sep = "-")
  
  if (max_param_length > 26L) {
    if (max_param_length > 26L^2) stop("Too many categories!")
    labs <- do.call(paste0, args = rev(expand.grid(LETTERS, LETTERS[1:ceiling(max_param_length / 26)])))[1:max_param_length]
  } else {
    labs <- LETTERS[1:max_param_length] 
  }
  
  studies$X <- factor(studies$X, labels = labs)
  
  return(studies)
  
}