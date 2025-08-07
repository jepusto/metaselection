#' @title Define prior functions for selection model parameters
#'
#' @description Creates a set of priors for use in estimating selection models
#'
#'
#' @param beta_mean numeric vector of prior means for beta (mean regression) parameters
#' @param beta_sd numeric vector of prior standard deviations for beta (mean regression) parameters
#' @param tau_mode numeric vector of prior modes for gamma (log-heterogeneity regression) parameters
#' @param tau_sd numeric vector of prior standard deviations for gamma (log-heterogeneity regression) parameters
#' @param zeta_mean numeric vector of prior means for zeta (selection) parameters
#' @param zeta_sd numeric vector of prior standard deviations for zeta (selection) parameters
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
#' 

define_priors <- function(
    beta_mean = 0,
    beta_sd = 1,
    tau_mode = 0.15,
    tau_sd = 0.15,
    zeta_mean = 0,
    zeta_sd = 2,
    zeta0_mean = 0,
    zeta0_sd = 2,
) {
  
  tau_lambda <- (tau_mode + sqrt(tau_mode^2 + 4 * tau_sd^2)) / (2 * tau_sd^2)
  tau_alpha <- tau_lambda * tau_mode + 1
  
  log_prior <- function(beta, gamma, zeta) {
    beta_log_prior <- dnorm(beta, mean = beta_mean, sd = beta_sd, log = TRUE)
    gamma_log_prior <- dgamma(exp(gamma), shape = tau_alpha, rate = tau_lambda, log = TRUE)
    zeta_log_prior <- dnorm(zeta, mean = zeta_mean, sd = zeta_sd, log = TRUE)
    sum(beta_log_prior) + sum(gamma_log_prior) + sum(zeta_log_prior)
  }
  
  score_prior <- function(beta, gamma, zeta) {
    score_beta <- (beta - beta_mean) / beta_sd^2
    score_gamma <- ()
    score_zeta <- (zeta - zeta_mean) / zeta_sd^2
    c(score_beta, score_gamma, score_zeta)
  }
  
  hessian_prior   <- function(beta, gamma, zeta) {
    beta_hessian <- rep(1 / beta_sd^2, length.out = length(beta))
    gamma_hessian <- ()
    zeta_hessian <- rep(1 / zeta_sd^2, length.out = length(zeta))
    diag(c(beta_hessian, gamma_hessian, zeta_hessian))
  }
}
