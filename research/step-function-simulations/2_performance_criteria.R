#------------------------------------------------------
# Calculate performance measures
# (For some simulations, it may make more sense
# to do this as part of the simulation driver.)
#------------------------------------------------------


calc_performance <- function(results, winz = Inf, B_target = 1999) {
  
  require(simhelpers)
  
  results_template <- 
    results %>%
    group_by(model, estimator, param) %>%
    summarize(.groups = "drop")
  
  if (all(c("ARGL","CML") %in% unique(results$estimator))) {
    results_fallback <- 
      results %>%
      group_by(rep, model, param) %>%
      filter(estimator %in% c("ARGL","CML")) %>%
      summarize(
        across(
          c(Est, SE, p_value, CI_lo, CI_hi, R_conv, max_method, true_param),
          \(x) if_else(is.na(x[estimator=="CML"]), x[estimator=="ARGL"], x[estimator=="CML"])
        ), 
        .groups = "drop"
      ) %>%
      mutate(
        estimator = "CML-fallback"
      )
    
    results <- bind_rows(results, results_fallback)
  }
  
    
  results_NA <- 
    results %>% 
    mutate(var_est = SE ^ 2) %>%
    group_by(model, estimator, param) 
  
  
  if (nrow(results) == 0) {
    
    res <- 
      results_template %>%
      mutate(
        K_absolute = 0,
        K_coverage = 0,
        K_relvar = 0,
      )
    
  } else {

    res <- 
      results_NA %>%
      summarize(
        calc_absolute(estimates = Est, true_param = true_param, criteria = c("bias", "variance", "rmse"), winz = winz),
        calc_relative_var(estimates = Est, var_estimates = var_est, criteria = c("relative bias","relative rmse"), winz = winz),
        calc_coverage(lower_bound = CI_lo, upper_bound = CI_hi, true_param = true_param, criteria = c("coverage","width"), winz = winz),
        .groups = "drop"
      )
    
    if ("boot_CIs" %in% names(results)) {
      
      boot_res <- 
        results_NA %>%
        filter(sapply(boot_CIs, is.data.frame)) %>%
        summarize(
          extrapolate_coverage(CI_subsamples = boot_CIs, true_param = true_param, B_target = B_target, winz = winz, cover_na_val = 0, width_na_val = Inf, nested = TRUE),
          .groups = "drop"
        )
      
      res <- left_join(res, boot_res, by = c("model","param","estimator"))
      
    }
  }

  return(res)
 
}