#' @title Print results from a `selmodel` object
#'
#' @description Print relevant results from a fitted `selmodel` object.
#' 
#' 
#' @param x Fitted model of class \code{"selmodel"}.
#' @param transf_gamma logical with `TRUE` (the default) indicating that the heterogeneity parameter estimates (called gamma) should be transformed by exponentiating.
#' @param transf_zeta logical with `TRUE` (the default) indicating that the selection parameter estimates (called zeta) should be transformed by exponentiating.
#' @param digits Minimum number of significant digits to be used, with a default of 3.
#' @param ... further arguments passed to \code{print.data.frame()}.
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
#'   estimator = "ML",
#'   bootstrap = "none"
#' )
#' 
#' print(res_ML)
#' print(res_ML, transf_gamma = FALSE, transf_zeta = FALSE)
#'  



print.selmodel <- function(x, transf_gamma = TRUE, transf_zeta = TRUE, digits = 3, ...) {
  
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
    estimates$param <- sub("^zeta","lambda", estimates$param)
  }
  
  col_names <- setdiff(names(estimates), c("estimator","bootstrap","bootstraps"))
  print(estimates[,col_names], digits = digits, ..., row.names = FALSE)
  
}