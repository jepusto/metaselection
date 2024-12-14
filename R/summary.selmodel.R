#' @title Summarize results from a `selmodel` object
#'
#' @description Summarize relevant results from a `selmodel` object.
#' 
#' @param object Fitted model of class \code{"selmodel"}.
#' @inheritParams print.selmodel
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
#'   estimator = "CML",
#'   bootstrap = "none"
#' )
#' 
#' summary(res_ML)
#' summary(res_ML, transf_gamma = FALSE, transf_zeta = FALSE)



summary.selmodel <- function(object, transf_gamma = TRUE, transf_zeta = TRUE, digits = 3, ...) {
  
  # title -------------------------------------------------------------------
  
  # need to edit the language for the title 
  model <- if (inherits(object, "step.selmodel")) "Step Function Model" else if (inherits(object, "beta.selmodel")) "Beta Density Model"
  
  # info  -------------------------------------------------------------------
  n_clusters <- object$n_clusters
  n_effects <- object$n_effects
  steps <- paste(object$steps, collapse = ", ")
  call <- object$cl
 
  
  estimates <- object$est
  estimator <- estimates$estimator[1]
  estimator <- ifelse(estimator %in% c("ML","CML"), "composite marginal likelihood", "augmented and reweighted Gaussian likelihood")
  vcov_type <- object$vcov_type
  clog_lik <- object$ll
  wt_partial_log_lik <- object$wpll
  
  # for the p table
  ptable <- object$ptable
  selection_mods <- !(is.null(object$selmods))
  steps_original <- object$steps
  

  # bootstrap information  --------------------------------------------------
                 
  if (inherits(object, "boot.selmodel")) {
    model <- paste(model, "with Cluster Bootstrapping")
    
    boot_type <- object$bootstrap_type[1]
    R <- unique(estimates$bootstraps)
  }
  

  all_vars <- c("param", "Est", "SE", "p_value", "CI_lo", "CI_hi", "percentile_lower", "percentile_upper", "basic_lower", "basic_upper", "student_lower", "student_upper")
  vars_display <- intersect(names(estimates), all_vars)
  
  # mean model results -----------------------------------------------------

  beta_params <- grepl("^beta", estimates$param)
  beta_estimates <- estimates[beta_params, vars_display]

  # heterogeneity -----------------------------------------------------------

  # do we need to have a transf_gamma argument or just display tau2 and the se 
  
  transf_variables <- intersect(names(estimates), setdiff(all_vars, c("param","SE", "p_value")))

  gamma_params <- grepl("^gamma", estimates$param)
  
  if (transf_gamma) {
    
    estimates[gamma_params, transf_variables] <- exp(estimates[gamma_params,transf_variables])
    estimates$SE[gamma_params] <- estimates$Est[gamma_params] * estimates$SE[gamma_params]
    estimates$param <- sub("^gamma","tau2", estimates$param)
    gamma_estimates <- estimates[gamma_params, vars_display]
  
  } else {
    
    gamma_estimates <- estimates[gamma_params, vars_display]
    
  }
  

  # zetas  ------------------------------------------------------------------

  zeta_params <- grepl("^zeta", estimates$param)
  

  
  if (transf_zeta) {
    
    estimates[zeta_params, transf_variables] <- exp(estimates[zeta_params, transf_variables])
    estimates$SE[zeta_params] <- estimates$Est[zeta_params] * estimates$SE[zeta_params]
    estimates$param <- sub("^zeta","lambda", estimates$param)
    
  }
  
  zeta_estimates <- estimates[zeta_params, vars_display]
  
  if (inherits(object, "step.selmodel")) {
    
    zeta_estimates <- clean_zetas(
      zeta_estimates_dat = zeta_estimates,
      transform_zeta = transf_zeta,
      selmods = selection_mods,
      pval_table = ptable,
      steps = steps_original,
      clusters = n_clusters
    )
      
  }


  # output ------------------------------------------------------------------

  cat(model, "\n", "\n")
  cat("Call:", "\n")
  print(call)
  cat("\n")
  if (!is.null(n_clusters)) {
    cat("Number of clusters = ", n_clusters, "; ", "Number of effects = ", n_effects, "\n", "\n", sep = "")
  } else {
    cat("Number of effects = ", n_effects, "\n", "\n", sep = "") 
  }
  
  cat("Steps:", steps, "\n")
  cat("Estimator:", estimator, "\n")
  cat("Variance estimator:", vcov_type, "\n")
  if (inherits(object, "boot.selmodel")) {
    cat("Bootstrap type:", boot_type,"\n")
    cat("Number of bootstrap replications:", R, "\n")
  }  
  
  cat("\n", "Log composite likelihood of selection model: ", clog_lik, "\n", sep = "")
  if (inherits(object, "step.selmodel")) {
    cat("Inverse selection weighted partial log likelihood:", wt_partial_log_lik, "\n")
  }

  cat("\n", "Mean effect estimates:", sep = "")
  print_with_header(beta_estimates, digits = digits, ...)
  
  cat("\n", "Heterogeneity estimates:", sep = "")
  print_with_header(gamma_estimates, digits = digits, ...)
  
  cat("\n", "Selection process estimates:", sep = "")
  print_with_header(zeta_estimates, digits = digits, ...)
}

print_with_header <- function(x, digits, ...) UseMethod("print_with_header")

#' @exportS3Method

print_with_header.data.frame <- function(x, digits, ...) {
 
  all_vars <- data.frame(
    param = c(" ", "Coef."),
    Est   = c(" ", "Estimate"),
    SE    = c(" ", "Std. Error"),
    p_value = c(" ", "p-value"),
    CI_lo = c("Large", "Lower"),
    CI_hi = c("Sample", "Upper"),
    percentile_lower = c("Percentile","Lower"),
    percentile_upper = c("Bootstrap", "Upper"),
    basic_lower = c("Basic","Lower"),
    basic_upper = c("Bootstrap", "Upper"),
    student_lower = c("Studentized","Lower"),
    student_upper = c("Bootstrap", "Upper")
  )
  
  
  vars_display <- intersect(names(all_vars), names(x))
  
  x_format <- format(x[,vars_display], digits = digits, ...)
  
  
  x_small <- x[,vars_display]
  x_format[is.na(x_small)] <- "---"

  x_print <- rbind(all_vars[,vars_display], x_format)
  
  print(unname(x_print), row.names = FALSE, na.print = "")
}

#' @exportS3Method

print_with_header.list <- function(x, digits, ...) {
  for (h in seq_along(x)) {
    cat("\n", names(x)[[h]])
    print_with_header.data.frame(x[[h]], digits = digits, ...)
  }
}



clean_zetas <- function(zeta_estimates_dat,
                        transform_zeta,
                        selmods,
                        pval_table,
                        steps,
                        clusters){
  
  if (transform_zeta) {
    
    param_name <- "lambda"
    
  } else{
    
    param_name <- "zeta"
    
  }
  
  
  first_step <- zeta_estimates_dat[1, ]
  first_step[1, ] <- NA
  first_step$Est <- as.numeric(transform_zeta)
  first_step$param <- paste0(param_name, "0")
  
  zeta_estimates <- rbind(first_step, zeta_estimates_dat)
  
  if (selmods == FALSE) {
    
    zeta_estimates <- cbind(pval_table, zeta_estimates)
    
  } else if (selmods == TRUE) {
    
    zeta_estimates$zeta <- sapply(strsplit(zeta_estimates$param, "_"),"[[",1)
    z_name <- paste0(param_name, 0:length(steps))
    pval_table$zeta <- z_name
    
    zeta_estimates <- merge(pval_table, zeta_estimates, by = "zeta")
    
    
  }
  
  if (!is.null(clusters)) {
    
    zeta_estimates$label <- with(zeta_estimates, paste0("Step: ", step, "; Studies: ", m, "; Effects: ", k))
  
  } else if (is.null(clusters)) {
    
    zeta_estimates$label <- with(zeta_estimates, paste0("Step: ", step, "; Effects: ", k))
    
    
  }
  
  
  zeta_estimates <- split(zeta_estimates, zeta_estimates$label)
  
  
  return(zeta_estimates)
  
}

