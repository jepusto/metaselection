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
                    censor_fun = step_fun,
                    n_ES_sim = n_ES_empirical(wwc_es),
                    steps,
                    mean_mods = NULL,
                    var_mods = NULL,
                    sel_mods = NULL,
                    sel_zero_mods = NULL,
                    make_sandwich = TRUE,
                    conf_level = .95,
                    stepfun_methods = c("step-MLE","step-hybrid-Broyden"), 
                    comparison_methods = "All", 
                    rho = cor_mu,
                    seed = NULL,
                    summarize_performance = TRUE,
                    filter_n_study = 0L,
                    ML_optimizer = "BFGS",
                    ML_optimizer_control = list(),
                    Broyden_optimizer_control = list(),
                    Newton_optimizer_control = list(),
                    rootSolve_optimizer_control = list(),
                    ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if (length(steps) == 1) {
    true_params <- data.frame(param = c("beta", "gamma", "omega1"))
  } else {
    true_params <- data.frame(param = c("beta", "gamma", paste0("omega", 1:length(steps))))
  }
  
  censor_fun_param <- censor_fun(...)
  total_var <- tau^2 + omega^2
  true_params$true_param <- c(mean_smd, log(total_var), log(censor_fun_param(steps + 1e-8)))
  
  
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
      
      p_counts <- pval_distribution(dat, steps = steps)
      
      res_selection <- estimate_step_models(
        dat = dat,
        estimators = stepfun_methods,
        steps = steps,
        sd_char = "sd_d",
        mean_mods = mean_mods,
        var_mods = var_mods,
        sel_mods = sel_mods,
        sel_zero_mods = sel_zero_mods,
        make_sandwich = make_sandwich,
        conf_level = conf_level,
        ML_optimizer = ML_optimizer,
        ML_optimizer_control = ML_optimizer_control,
        Broyden_optimizer_control = Broyden_optimizer_control,
        Newton_optimizer_control = Newton_optimizer_control,
        rootSolve_optimizer_control = rootSolve_optimizer_control
      )
      
      res_comp <- estimate_comparison_methods(
        dat = dat, 
        rho = rho, 
        methods = comparison_methods
      )
      
      bind_rows(res_selection, res_comp) %>%
        mutate(
          n_study_filter = all(p_counts$n_Study >= filter_n_study)
        )

    }, .id = "rep")
  
  if (!summarize_performance || (is.null(stepfun_methods) & is.null(comparison_methods))) {
    return(results)
  } else {
    summary_res <- 
      results %>%
      merge(true_params) %>%
      calc_performance()
    
    if (filter_n_study > 0L) {
      summary_res_filter <- 
        results %>%
        merge(true_params) %>%
        calc_performance(filter_n_study = TRUE) %>%
        mutate(filtered = TRUE)
      
      summary_res <- 
        summary_res %>%
        mutate(filtered = FALSE) %>%
        bind_rows(summary_res_filter)
      
    }

    return(summary_res)
  }

  
}

