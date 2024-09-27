#' @title Print results from `selmodel()`
#'
#' @description Print relevant results from `selmodel()`
#' 
#' 
#' @param x Fitted model of class \code{"selmodel"}.
#' @param transf_gamma logical with `TRUE` indicating that the heterogeneity parameter estimates (called gamma) should be transformed by exponentiating.
#' @param transf_zeta logical with `TRUE` indicating that the selection parameter estimates (called zeta) should be transformed by exponentiating.
#' @param digits Minimum number of significant digits to be used, with a default of 3.
#' @param ... further arguments passed to \code{print.data.frame()}.
#'
#' @export



print.selmodel <- function(x, transf_gamma = FALSE, transf_zeta = FALSE, digits = 3, ...) {
  
  estimates <- x$est
  
  if (transf_gamma | transf_zeta) {
    transf_variables <- intersect(
      names(estimates),
      c("Est","CI_lo","CI_hi","percentile_lower","percentile_upper","basic_lower","basic_upper","student_lower","student_upper")
    )
  }
  if (transf_gamma) {
    gamma_params <- grepl("^gamma", estimates$param)
    estimates[gamma_params, transf_variables] <- exp(estimates[gamma_params,transf_variables])
    estimates$SE[gamma_params] <- estimates$Est[gamma_params] * estimates$SE[gamma_params]
    estimates$param <- sub("^gamma","tau2", estimates$param)
  }
  
  if (transf_zeta) {
    zeta_params <- grepl("^zeta", estimates$param)
    estimates[zeta_params, transf_variables] <- exp(estimates[zeta_params, transf_variables])
    estimates$SE[zeta_params] <- estimates$Est[zeta_params] * estimates$SE[zeta_params]
    estimates$param <- sub("^zeta","lambda_", estimates$param)
  }
  
  col_names <- setdiff(names(estimates), c("estimator","bootstrap","bootstraps"))
  print(estimates[,col_names], digits = digits, ...)
  
  return(estimates)

}