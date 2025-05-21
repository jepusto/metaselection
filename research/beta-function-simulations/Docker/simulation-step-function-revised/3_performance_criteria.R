#------------------------------------------------------
# Calculate performance measures
# (For some simulations, it may make more sense
# to do this as part of the simulation driver.)
#------------------------------------------------------


calc_performance <- function(results, filter_n_study = FALSE) {
  
  require(simhelpers)
  
  if (filter_n_study) results <- filter(results, n_study_filter)
  
  results_NA <- 
    results %>% 
    mutate(var = SE ^ 2) %>%
    group_by(estimator, param, sd) 
  

  abs <- 
    results_NA %>%
    group_modify(~ calc_absolute(.x, estimates = Est, 
                                 true_param = true_param,
                                 criteria = c("bias", "variance", "rmse"))) %>%
    ungroup()
  
  ci_cov <- 
    results_NA %>%
    group_modify(~ calc_coverage(.x, lower_bound = CI_lo, 
                                 upper_bound = CI_hi, 
                                 true_param = true_param,
                                 criteria = c("coverage","width"))) %>%
    ungroup()
  
  
  rel_var <- 
    results_NA %>%
    group_modify(~ calc_relative_var(.x, 
                                     estimates = Est,
                                     var_estimates = var)) %>%
    ungroup()
 
  best_methods <- 
    results %>%
    filter(!is.na(max_method)) %>%
    group_by(estimator, param, sd, max_method) %>%
    summarize(
      n = n(), .groups = "drop_last"
    ) %>%
    mutate(
      pr = n / sum(n)
    ) %>%
    select(-n) %>%
    pivot_wider(names_from = max_method, values_from = pr) %>%
    nest(.key = "best_methods")
  
  
  abs %>%
    left_join(ci_cov, by = c("estimator", "param","sd")) %>%
    left_join(rel_var, by = c("estimator", "param","sd")) %>%
    left_join(best_methods, by = c("estimator","param","sd"))
  
}