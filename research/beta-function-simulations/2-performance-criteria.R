#------------------------------------------------------
# Calculate performance measures
# (For some simulations, it may make more sense
# to do this as part of the simulation driver.)
#------------------------------------------------------

possibly_extrapolate_coverage <- purrr::possibly(simhelpers::extrapolate_coverage, otherwise = tibble::tibble(.rows = 1L))

calc_performance <- function(results, winz = Inf, B_target = 1999) {
  
  require(simhelpers)
  
  results_template <- 
    results %>%
    group_by(model, estimator, param) %>%
    summarize(.groups = "drop")
  
  results_all <- 
    results %>% 
    filter(param %in% c("gamma", "zeta1", "zeta2")) %>%
    mutate(
      across(c(Est, CI_lo, CI_hi, true_param), ~ exp(.x)),
      SE = SE * Est,
      param = case_match(param, "gamma" ~ "tau2", "zeta1" ~ "lambda1", "zeta2" ~ "lambda2")
    ) %>%
    bind_rows(results) %>%
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
      results_all %>%
      summarize(
        calc_absolute(estimates = Est, true_param = true_param, criteria = c("bias", "variance", "rmse"), winz = winz),
        calc_relative_var(estimates = Est, var_estimates = var_est, criteria = c("relative bias","relative rmse"), winz = winz),
        calc_coverage(lower_bound = CI_lo, upper_bound = CI_hi, true_param = true_param, criteria = c("coverage","width"), winz = winz),
        .groups = "drop"
      )  
    
    if ("boot_CIs" %in% names(results)) {
      
      boot_res <- 
        results_all %>%
        filter(sapply(boot_CIs, is.data.frame)) %>%
        mutate(max_boots = sapply(boot_CIs, \(x) max(x$bootstraps))) %>%
        filter(max_boots == max(max_boots)) %>%
        summarize(
          possibly_extrapolate_coverage(CI_subsamples = boot_CIs, true_param = true_param, B_target = B_target, winz = winz, cover_na_val = 0, width_na_val = Inf, nested = TRUE),
          .groups = "drop"
        )
      
      res <- left_join(res, boot_res, by = c("model","param","estimator"))
      
    }
  }
  
  return(res)
  
}
