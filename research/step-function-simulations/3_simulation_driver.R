#-----------------------------------------------------------
# Simulation Driver
#-----------------------------------------------------------

run_sim <- function(
    iterations, 
    mean_smd,
    tau,
    omega,
    m, 
    cor_mu,
    cor_sd, 
    censor_fun = step_fun,
    steps,
    n_ES_dat = wwc_es,
    n_multiplier = 1,
    estimate_4psm = FALSE, 
    mean_mods = NULL,
    var_mods = NULL,
    sel_mods = NULL,
    sel_zero_mods = NULL,
    priors = "None",
    conf_level = .95,
    stepfun_methods = c("CML","ARGL"), 
    comparison_methods = "All", 
    rho = cor_mu,
    seed = NULL,
    summarize_performance = TRUE,
    filter_n_study = 2L,
    CML_optimizer = "BFGS",
    CML_optimizer_control = list(),
    ARGL_optimizer_control = list(),
    bootstrap = "multinomial",
    CI_type = c("large-sample","basic","percentile","student","bias-corrected","BCa"),
    R = c(49,99,199,299),
    retry_bootstrap = 3L,
    winz = Inf,
    B_target = 1999,
    ...
) {
  
  require(metaselection)
  
  prior_spec <- switch(
    priors, 
    "None" = NULL, 
    "Default" = default_priors(
      beta_mean = 0,
      beta_precision = 1 / 2,
      beta_L = 2,
      tau_mode = 0.2,
      tau_alpha = 1,
      lambda_mode = 0.8,
      lambda_precision = 1 / 2,
      lambda_L = 2
    ),
    "Weaker" = default_priors(
      beta_mean = 0,
      beta_precision = 2, 
      tau_mode = 0.2, 
      tau_alpha = 0.5,
      lambda_mode = 0.8,
      lambda_precision = 1 / 27,
      lambda_L = 4
    )
  )
  
  if (!is.null(seed)) set.seed(seed)
  
  if (length(steps) == 1) {
    true_params <- data.frame(param = c("beta", "gamma", "zeta1"))
  } else {
    true_params <- data.frame(param = c("beta", "gamma", paste0("zeta", 1:length(steps))))
  }
  
  censor_fun_param <- censor_fun(...)
  total_var <- tau^2 + omega^2
  true_params$true_param <- c(mean_smd, log(total_var), log(censor_fun_param(steps + 1e-8)))
  
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
      
      if (is.null(stepfun_methods) & is.null(comparison_methods)) return(dat)
      
      res_selection <- estimate_step_models(
        dat = dat,
        estimators = stepfun_methods,
        steps = 0.025,
        mean_mods = mean_mods,
        var_mods = var_mods,
        sel_mods = sel_mods,
        sel_zero_mods = sel_zero_mods,
        priors = prior_spec,
        CML_optimizer = CML_optimizer,
        CML_optimizer_control = CML_optimizer_control,
        ARGL_optimizer_control = ARGL_optimizer_control,
        CI_type = CI_type,
        bootstrap = bootstrap,
        conf_level = conf_level,
        R = R,
        retry_bootstrap = retry_bootstrap
      ) %>%
        mutate(model = "3PSM")
      
      if (estimate_4psm == TRUE) {
      
        res_selection_4psm <- estimate_step_models(
          dat = dat,
          estimators = stepfun_methods,
          steps = c(0.025, .5),
          mean_mods = mean_mods,
          var_mods = var_mods,
          sel_mods = sel_mods,
          sel_zero_mods = sel_zero_mods,
          priors = prior_spec,
          CML_optimizer = CML_optimizer,
          CML_optimizer_control = CML_optimizer_control,
          ARGL_optimizer_control = ARGL_optimizer_control,
          CI_type = CI_type,
          bootstrap = bootstrap,
          conf_level = conf_level,
          R = R,
          retry_bootstrap = retry_bootstrap
        ) %>%
          mutate(model = "4PSM")
        
        res_selection <- bind_rows(res_selection, res_selection_4psm)
      
      }
      
      res_comp <- estimate_comparison_methods(
        dat = dat, 
        rho = rho, 
        methods = comparison_methods
      ) %>%
        mutate(model = "Comparison")
      
      bind_rows(res_selection, res_comp)

    }, .id = "rep")
  
  if (!summarize_performance || (is.null(stepfun_methods) & is.null(comparison_methods))) {
    
    return(results)
    
  } else {
    
    summary_res <- 
      results %>%
      merge(true_params) %>%
      calc_performance(winz = winz, B_target = B_target)
    
    return(summary_res)
    
  }

  
}

