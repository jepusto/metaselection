#' @title Define default prior functions for selection model parameters
#'
#' @description Creates a set of priors for use in estimating selection models.
#'   beta parameters are assigned independent normal priors. log(tau) parameters
#'   are assigned independent log-gamma priors. log(lambda) parameters are
#'   assigned independent log-gamma priors.
#'
#'
#' @param beta_mean numeric vector of prior means for beta (mean regression)
#'   parameters.
#' @param beta_sd numeric vector of prior standard deviations for beta (mean
#'   regression) parameters.
#' @param tau_mode numeric vector of prior modes for tau (heterogeneity SD)
#'   regression parameters.
#' @param tau_alpha numeric vector of prior precisions for tau (heterogeneity SD)
#'   regression parameters.
#' @param lambda_mode numeric vector of prior modes for lambda (selection)
#'   parameters.
#' @param lambda_precision numeric vector of prior precisions for lambda
#'   (selection) parameters.
#'
#' @returns An object of class \code{"selmodel_prior"} containing the following
#'   components:
#' \describe{
#'   \item{\code{log_prior}}{A function with arguments \code{beta},\code{gamma},\code{zeta} that returns the log of the prior density over these parameters.}
#'   \item{\code{score_prior}}{A function with arguments \code{beta},\code{gamma},\code{zeta} that returns the vector of scores for the prior density over these parameters.}
#'   \item{\code{hessian_prior}}{A function with arguments \code{beta},\code{gamma},\code{zeta} that returns the Hessian matrix of the prior density over these parameters.}
#' }
#'
#' @export
#'
#' @examples
#' # set very informative priors on beta and lambda
#' strong_priors <- default_priors(beta_mean = 0.4, beta_sd = 0.1, lambda_mode = 0.2, lambda_sd = 0.1)
#'
#' # set very weak priors
#' weak_priors <- default_priors(beta_sd = 10, tau_sd = 1, lambda_sd = 4)
#' 

default_priors <- function(
    beta_mean = 0,
    beta_precision = 1 / 2,
    beta_L = 2,
    tau_mode = 0.20,
    tau_alpha = 1,
    lambda_mode = 0.8,
    lambda_precision = 1 / 2,
    lambda_L = 2
) {
  
  gamma_lambda <- tau_alpha / tau_mode
  gamma_alpha <- tau_alpha
  
  zeta_mean <- log(lambda_mode)
  
  log_prior <- function(beta, gamma, zeta, include_zeta = TRUE) {
    beta_log_prior <- -beta_precision * abs(beta - beta_mean)^beta_L
    gamma_log_prior <- 0.5 * gamma_alpha * gamma - gamma_lambda * exp(gamma / 2)
    if (include_zeta) {
      zeta_log_prior <- -lambda_precision * (zeta - zeta_mean)^lambda_L
      sum(beta_log_prior) + sum(gamma_log_prior) + sum(zeta_log_prior)
    } else {
      sum(beta_log_prior) + sum(gamma_log_prior)
    }
    
  }
  
  score_prior <- function(beta, gamma, zeta) {
    score_beta <- -1 * beta_precision * beta_L * sign(beta - beta_mean) * abs(beta - beta_mean)^(beta_L - 1)
    score_gamma <- 0.5 * (gamma_alpha - gamma_lambda * exp(gamma / 2))
    score_zeta <- - 1 * lambda_precision * lambda_L * sign(zeta - zeta_mean) * abs(zeta - zeta_mean)^(lambda_L - 1)
    c(score_beta, score_gamma, score_zeta)
  }
  
  hessian_prior <- function(beta, gamma, zeta) {
    beta_hessian <- -1 * beta_precision * beta_L * (beta_L - 1) * abs(beta - beta_mean)^(beta_L - 2)
    gamma_hessian <- -0.25 * gamma_lambda * exp(gamma / 2)
    zeta_hessian <- -1 * lambda_precision * lambda_L * (lambda_L - 1) * abs(zeta - zeta_mean)^(lambda_L - 2)
    diag(c(beta_hessian, gamma_hessian, zeta_hessian))
  }
  
  res <- list(
    log_prior = log_prior,
    score_prior = score_prior,
    hessian_prior = hessian_prior
  )
  
  class(res) <- "selmodel_prior"
  
  return(res)
}
