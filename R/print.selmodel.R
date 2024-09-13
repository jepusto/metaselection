#' @title Print results from `selmodel()`
#'
#' @description Print relevant results from `selmodel()`
#' 
#' 
#' @param x Fitted model of class \code{"selmodel"}.
#' @param exponentiate logical with `TRUE` indicating that the gamma and zeta estimates should be exponentiated
#' @param ... Not used.
#'
#' @export



print.selmodel <- function(x, exponentiate = TRUE, ...) {
  
  model <- if ("step.selmodel" %in% class(x)) "Step Function Model with Robust Variance Estimation"
  
  estimates <- x$est
  
  if (exponentiate) {
    
    estimates$Est[estimates$param == "gamma"] <- exp(estimates$Est)[estimates$param == "gamma"] 
    
  }
  
  call <- x$cl
  steps <- x$steps
  
  cat(model, "\n")
  # cat("\nCall:", call ,"\n\n")
  print(estimates)
  cat("---\n")
  cat("Signif. codes: < .01 *** < .05 ** < .10 *\n")
  cat("---\n")
}