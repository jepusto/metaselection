#-----------------------------------------------------------
# Simulation Driver
#-----------------------------------------------------------

run_sim <- function(iterations, 
                    mean_smd,
                    tau,
                    omega,
                    m, 
                    cor_mu,
                    cor_sd, 
                    delta_1, 
                    delta_2,
                    censor_fun = beta_fun,
                    n_ES_dat = wwc_es,
                    n_multiplier = 1,
                    steps = c(.025, .975),
                    mean_mods = NULL,
                    var_mods = NULL,
                    sel_mods = NULL,
                    sel_zero_mods = NULL,
                    conf_level = .95,
                    step_models = c("3PSM","4PSM"),
                    comparison_methods = "All", 
                    rho = cor_mu,
                    seed = NULL,
                    summarize_performance = TRUE,
                    winz = 2.5,
                    bootstrap = "two-stage",
                    boot_CI = c("large-sample","basic","percentile","student","bias-corrected"),
                    R_beta = c(49,99,199,299,399),
                    retry_bootstrap = 3L,
                    R_step = R_beta,
                    ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  
  total_var <- tau^2 + omega^2
  true_params <- data.frame(
    param = c("beta", "gamma", "zeta1", "zeta2"),
    true_param = c(mean_smd, log(total_var), log(delta_1), log(delta_2))
  )
  
  censor_fun_param <- censor_fun(delta_1 = delta_1, delta_2 = delta_2)
  
  n_ES_dat$n <- round(n_multiplier * n_ES_dat$n)
  n_ES_sim <- n_ES_empirical(n_ES_dat)
  
  
  results <-
    
    map_dfr(1:iterations, ~ {
      
      dat <- r_meta(mean_smd = mean_smd,
                    tau = tau,
                    omega = omega,
                    m = m,
                    cor_mu = cor_mu,
                    cor_sd = cor_sd,
                    censor_fun =  censor_fun_param,  
                    n_ES_sim = n_ES_sim)
      
      if (is.null(beta_methods) & is.null(step_models) & is.null(comparison_methods)) return(dat)
      
      res_selection <- estimate_step_models(
        dat = dat,
        selection_type = "beta",
        steps = steps,
        sd_char = "sd_d",
        mean_mods = mean_mods,
        var_mods = var_mods,
        sel_mods = sel_mods,
        sel_zero_mods = sel_zero_mods,
        conf_level = conf_level,
        bootstrap = bootstrap,
        boot_CI = boot_CI,
        R = R_beta,
        retry_bootstrap = retry_bootstrap
      ) %>%
        mutate(model = "beta")
      
      if ("3PSM" %in% step_models) {
        res_3PSM <- estimate_step_models(
          dat = dat,
          selection_type = "step",
          steps = .025,
          sd_char = "sd_d",
          mean_mods = mean_mods,
          var_mods = var_mods,
          sel_mods = sel_mods,
          sel_zero_mods = sel_zero_mods,
          conf_level = conf_level,
          bootstrap = bootstrap,
          boot_CI = boot_CI,
          R = R_step,
          retry_bootstrap = retry_bootstrap
        ) %>%
          mutate(model = "3PSM")
        
        res_selection <- bind_rows(
          res_selection, 
          res_3PSM
        )
      }
      
      if ("4PSM" %in% step_models) {
        res_4PSM <- estimate_step_models(
          dat = dat,
          selection_type = "step",
          steps = c(.025,.50),
          sd_char = "sd_d",
          mean_mods = mean_mods,
          var_mods = var_mods,
          sel_mods = sel_mods,
          sel_zero_mods = sel_zero_mods,
          conf_level = conf_level,
          bootstrap = bootstrap,
          boot_CI = boot_CI,
          R = R_step,
          retry_bootstrap = retry_bootstrap
        ) %>%
          mutate(model = "4PSM")
        
        res_selection <- bind_rows(
          res_selection, 
          res_4PSM
        )
      }
      res_comp <- estimate_comparison_methods(
        dat = dat, 
        rho = rho, 
        methods = comparison_methods
      ) %>%
        mutate(model = "Comparison")
      
      bind_rows(res_selection, res_comp)
      
    }, .id = "rep")
  
  if (!summarize_performance || (is.null(beta_methods) & is.null(comparison_methods))) {
    
    return(results)
    
  } else {
    
    summary_res <- 
      results %>%
      merge(true_params) %>%
      calc_performance(winz = winz)
    
    return(summary_res)
    
  }

  
}

