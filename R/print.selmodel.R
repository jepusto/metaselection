#' @title Print results from `selmodel()`
#'
#' @description Print relevant results from `selmodel()`
#' 
#' 
#' @param x Object from selmodel class.
#' @param exponentiate logical with `TRUE` indicating that the gamma and zeta estimates should be exponentiated
#'
#'
#' @export



print.selmodel <- function(mod, expoentiate = TRUE){
  
  model <- if("step.selmodel" %in% class(mod)) "Step Function Model with Robust Variance Estimation"
  
  estimates <- mod$est
  
  if(expoentiate){
    
    estimates$Est[estimates$param == "gamma"] <- exp(estimates$Est)[estimates$param == "gamma"] 
    
  }
  
  call <- mod$cl
  steps <- mod$steps
  
  cat(model, "\n")
  # cat("\nCall:", call ,"\n\n")
  print(estimates)
  cat("---\n")
  cat("Signif. codes: < .01 *** < .05 ** < .10 *\n")
  cat("---\n")
}