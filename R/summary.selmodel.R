#' @title Summarize results from `selmodel()`
#'
#' @description Summarize relevant results from `selmodel()`
#' 
#' 
#' @param x Fitted model of class \code{"selmodel"}.
#' @param transf_gamma logical with `TRUE` indicating that the heterogeneity parameter estimates (called gamma) should be transformed by exponentiating.
#' @param transf_zeta logical with `TRUE` indicating that the selection parameter estimates (called zeta) should be transformed by exponentiating.
#' @param ... Not used.
#'
#' @export



summary.selmodel <- function(x, digits = 3, transf_gamma = FALSE, transf_zeta = FALSE, ...) {
  
  # title -------------------------------------------------------------------
  
  # need to edit the language for the title 
  model <- if ("step.selmodel" %in% class(x)) "Step Function Model" else if("beta.selmodel" %in% class(x)) "Beta Density Model"
  boot_model <- if(grepl("^boot", class(x)[1])) "With Cluster Bootstrapping" else NULL
  if(!is.null(boot_model)) model <- paste(model, boot_model)

  # info  -------------------------------------------------------------------
  n_clusters <- x$n_clusters
  n_effects <- x$n_effects
  steps <- paste(x$steps, collapse = ", ")
  call <- x$cl
  
  estimates <- x$est
  estimator <- estimates$estimator[1]
  estimator <- ifelse(estimator == "ML", "Maximum likelihood", "Hybrid estimating equations")
  

  # bootstrap information  --------------------------------------------------
                 
  if(grepl("^boot", class(x)[1])){
    
    boot_type <- estimates$bootstrap[1]
    R <- estimates$bootstraps[2]
    
    boot_ci_type_names <- names(estimates)[endsWith(names(estimates), 'lower')]
    boot_ci_type <- sapply(strsplit(boot_ci_type_names, "_"), `[`, 1)
    boot_ci_type <- paste(boot_ci_type, collapse = ", ")

  }
  

  # betas model results -----------------------------------------------------

  vars_display <- intersect(
    names(estimates),
    c("param", "Est", "SE", "CI_lo","CI_hi","percentile_lower","percentile_upper","basic_lower","basic_upper","student_lower","student_upper")
  )

  beta_params <- grepl("^beta", estimates$param)
  beta_estimates <- estimates[beta_params, vars_display]
  beta_estimates <- format(beta_estimates, digits = digits)

  # heterogeneity -----------------------------------------------------------

  # do we need to have a transf_gamma argument or just display tau2 and the se 
  
  transf_variables <- intersect(
    names(estimates),
      c("Est","CI_lo","CI_hi","percentile_lower","percentile_upper","basic_lower","basic_upper","student_lower","student_upper")
    )

  gamma_params <- grepl("^gamma", estimates$param)
  
  if(transf_gamma){
    
    estimates[gamma_params, transf_variables] <- exp(estimates[gamma_params,transf_variables])
    estimates$SE[gamma_params] <- estimates$Est[gamma_params] * estimates$SE[gamma_params]
    estimates$param <- sub("^gamma","tau2", estimates$param)
    gamma_estimates <- estimates[gamma_params, vars_display]
  
  } else {
    
    gamma_estimates <- estimates[gamma_params, vars_display]
    
  }
  
  gamma_estimates <- format(gamma_estimates, digits = digits)
  

  # zetas  ------------------------------------------------------------------

  zeta_params <- grepl("^zeta", estimates$param)
  
  if(transf_zeta){
    
    estimates[zeta_params, transf_variables] <- exp(estimates[zeta_params, transf_variables])
    estimates$SE[zeta_params] <- estimates$Est[zeta_params] * estimates$SE[zeta_params]
    estimates$param <- sub("^zeta","lambda_", estimates$param)
    zeta_estimates <- estimates[zeta_params, vars_display]

  } else{
    
    zeta_estimates <- estimates[zeta_params, vars_display]
      
  }
  
  zeta_estimates <- format(zeta_estimates, digits = digits)
  

  # output ------------------------------------------------------------------

  cat(model, "\n", "\n")
  cat("Call:", "\n")
  print(call)
  cat("\n")
  if(!is.null(n_clusters)) cat(paste0("Number of clusters = ", n_clusters, ";", " Number of effects = ", n_effects), "\n", "\n")
  if(is.null(n_clusters)) cat(paste0("Number of effects = ", n_effects), "\n", "\n") 
  cat(paste("Steps: ", steps), "\n")
  cat(paste("Estimator: ", estimator), "\n")
  if(grepl("^boot", class(x)[1]))  cat(paste0("Bootstrap type: ", boot_type, "; Number of replications: ", R, "; CI type: ", boot_ci_type), "\n")
  cat("\n")
  cat("Model results:", "\n")
  print(beta_estimates, row.names = FALSE)
  cat("\n", "\n")
  cat("Heterogeneity estimates:", "\n")
  print(gamma_estimates, row.names = FALSE)
  cat("\n", "\n")
  cat("Selection model results:", "\n")
  print(zeta_estimates, row.names = FALSE)
  
  
}
