#-------------------------------------------------------------------------------
# Model-fitting/estimation/testing functions

pval_distribution <- function(dat, steps) {
  p_buckets <- cut(dat$p_onesided, breaks = c(0, steps, 1), include.lowest = TRUE, right = FALSE)
  es_tab <- table(p_buckets)
  pr_es <- prop.table(es_tab)
  study_tab <- table(dat$studyid, p_buckets) > 0L
  n_study <- colSums(study_tab)
  pr_study <- colMeans(study_tab)
  data.frame(
    range = names(pr_es),
    n_ES = as.numeric(es_tab),
    Pr_ES = as.numeric(pr_es),
    n_Study = as.numeric(n_study),
    pr_Study = as.numeric(pr_study)
  )
}
                                        
estimate_step_models <- function(
    dat,
    estimators = c("step-MLE", "step-hybrid-Broyden", "step-hybrid-Newton","step-hybrid-rootSolve"),
    steps,
    sd_char = "sd_d",
    mean_mods = NULL,
    var_mods = NULL,
    sel_mods = NULL,
    sel_zero_mods = NULL,
    make_sandwich = TRUE,
    conf_level = .95,
    ML_optimizer = "BFGS",
    ML_optimizer_control = list(),
    Broyden_optimizer_control = list(),
    Newton_optimizer_control = list(),
    rootSolve_optimizer_control = list()
) {
  
  error_res <- list(
    est = data.frame(
            param = c("beta","gamma",paste0("omega", seq_along(steps))),
            Est = rep(NA_real_, 2 + length(steps))
          ),
    info = list(convcode = NA_integer_, termcd = NA_integer_)
  )
  
  sd <- as.symbol(sd_char)
  
  res <- list()
  
  if ("step-MLE" %in% estimators) {
    res_MLE <- tryCatch(
      eval(substitute(
        selection_model(
          data = dat,
          yi = d,
          sei = sd,
          pi = p_onesided,
          cluster = studyid,
          steps = steps,
          mean_mods = mean_mods,
          var_mods = var_mods,
          sel_mods = sel_mods,
          sel_zero_mods = sel_zero_mods,
          make_sandwich = make_sandwich,
          conf_level = conf_level,
          estimator = "ML",
          optimizer = ML_optimizer, 
          optimizer_control = ML_optimizer_control
        ),
        list(sd = sd)
      )), error = function(e) error_res
    )
    res_MLE$est$R_conv <- res_MLE$info$convcode
    res_MLE$est$max_method = res_MLE$method
    res_MLE$est$estimator <- "ML"
    if (!is.na(res_MLE$info$convcode)) {
      res_MLE$est$score <- step_score(
        theta = res_MLE$est$Est, 
        yi = dat$d,
        sei = dat$sd_d,
        pi = dat$p_onesided,
        steps = steps
      )
      res_MLE$est$hybrid <- step_hybrid_score(
        theta = res_MLE$est$Est, 
        yi = dat$d,
        sei = dat$sd_d,
        pi = dat$p_onesided,
        steps = steps
      )
    }
    res$MLE <- res_MLE$est
  }
  
  if ("step-hybrid-Broyden" %in% estimators) {
    
    res_hybrid <- tryCatch(
      eval(substitute(
        selection_model(
          data = dat,
          yi = d,
          sei = sd,
          pi = p_onesided,
          cluster = studyid,
          steps = steps,
          mean_mods = mean_mods,
          var_mods = var_mods,
          sel_mods = sel_mods,
          sel_zero_mods = sel_zero_mods,
          make_sandwich = make_sandwich,
          conf_level = conf_level,
          estimator = "hybrid",
          optimizer = "nleqslv", 
          optimizer_control = c(method = "Broyden", Broyden_optimizer_control)
        ), 
        list(sd = sd)
      )), error = function(e) error_res
    )
    
    res_hybrid$est$estimator <- "hybrid-Broyden"
    if ("nleqslv" %in% names(res_hybrid$info)) {
      res_hybrid$est$R_conv <- res_hybrid$info$nleqslv$termcd
      res_hybrid$est$score <- step_score(
        theta = res_hybrid$est$Est, 
        yi = dat$d,
        sei = dat$sd_d,
        pi = dat$p_onesided,
        steps = steps
      )
      res_hybrid$est$hybrid <- step_hybrid_score(
        theta = res_hybrid$est$Est, 
        yi = dat$d,
        sei = dat$sd_d,
        pi = dat$p_onesided,
        steps = steps
      )
    }
    res$hybrid_Broyden <- res_hybrid$est
  }
  
  if ("step-hybrid-Newton" %in% estimators) {
    
    res_hybrid <- tryCatch(
      eval(substitute(
        selection_model(
          data = dat,
          yi = d,
          sei = sd,
          pi = p_onesided,
          cluster = studyid,
          steps = steps,
          mean_mods = mean_mods,
          var_mods = var_mods,
          sel_mods = sel_mods,
          sel_zero_mods = sel_zero_mods,
          make_sandwich = make_sandwich,
          conf_level = conf_level,
          estimator = "hybrid",
          optimizer = "nleqslv", 
          optimizer_control = c(method = "Newton", Newton_optimizer_control)
        ), 
        list(sd = sd)
      )), error = function(e) error_res
    )
    
    res_hybrid$est$estimator <- "hybrid-Newton"
    if ("nleqslv" %in% names(res_hybrid$info)) {
      res_hybrid$est$R_conv <- res_hybrid$info$nleqslv$termcd
      res_hybrid$est$score <- step_score(
        theta = res_hybrid$est$Est, 
        yi = dat$d,
        sei = dat$sd_d,
        pi = dat$p_onesided,
        steps = steps
      )
      res_hybrid$est$hybrid <- step_hybrid_score(
        theta = res_hybrid$est$Est, 
        yi = dat$d,
        sei = dat$sd_d,
        pi = dat$p_onesided,
        steps = steps
      )
    }
    res$hybrid_Newton <- res_hybrid$est
  }
  
  if ("step-hybrid-rootSolve" %in% estimators) {
    
    res_hybrid <- tryCatch(
      eval(substitute(
        selection_model(
          data = dat,
          yi = d,
          sei = sd_d,
          pi = p_onesided,
          cluster = studyid,
          steps = steps,
          mean_mods = mean_mods,
          var_mods = var_mods,
          sel_mods = sel_mods,
          sel_zero_mods = sel_zero_mods,
          make_sandwich = make_sandwich,
          conf_level = conf_level,
          estimator = "hybrid",
          optimizer = "rootSolve", 
          optimizer_control = rootSolve_optimizer_control
        ),
        list(sd = sd)
      )), error = function(e) error_res
    )
    
    res_hybrid$est$estimator <- "hybrid-rootSolve"
    if ("rootsolve" %in% names(res_hybrid$info)) {
      res_hybrid$est$R_conv <- as.integer(res_hybrid$info$rootsolve$f_norm > 1e-10)
      res_hybrid$est$score <- step_score(
        theta = res_hybrid$est$Est, 
        yi = dat$d,
        sei = dat$sd_d,
        pi = dat$p_onesided,
        steps = steps
      )
      res_hybrid$est$hybrid <- step_hybrid_score(
        theta = res_hybrid$est$Est, 
        yi = dat$d,
        sei = dat$sd_d,
        pi = dat$p_onesided,
        steps = steps
      )
    }
    res$hybrid_rootSolve <- res_hybrid$est
  }
  
  bind_rows(res) %>%
    mutate(sd = sd_char)
}


# che, fec, pet, and peese -----------------------------------------------------


estimate_meta  <- function(dat, 
                           rho = .8,
                           rand = NULL,
                           moderator = ~ 1,  # change to sd or var for pet peese
                           V_mat = NULL,
                           method) { 
  require(metafor)

  if (is.null(V_mat)) {
    V_mat <- vcalc(var_avg, 
                   cluster = studyid, obs = esid, 
                   data = dat, rho = rho)  
  }
  
  optimizers <- c("nlminb","nloptr","Rvmmin","BFGS")
  mod <- "Non-converged"
  i <- 1L
  
  while (!inherits(mod, "rma.mv") & i <= 4L) {
    mod <- tryCatch(
      rma.mv(
        yi = d,
        V = V_mat,
        random = rand,
        mods = moderator,
        data = dat,
        sparse = TRUE,
        control = list(optimizer=optimizers[i])
      ),
      error = function(e) "Non-converged"
    )
    i <- i + 1L
  }
  
  if (inherits(mod, "rma.mv")) {
    mod <- robust(mod, cluster = studyid, clubSandwich = TRUE)
    
    res <- data.frame(
      param = "beta",
      Est = as.numeric(mod$beta[1,]),
      SE = mod$se[1],
      df = mod$dfs[1],
      CI_lo = mod$ci.lb[1],
      CI_hi = mod$ci.ub[1],
      p = mod$pval[1],
      R_conv = NA,
      estimator = method
    )
  } else {
    res <- NULL
  }

  return(res)
  
}

estimate_comparison_methods <- function(dat, rho = .8, methods = "All") {
  
  res <- list()
  
  if (is.null(methods)) return(res)
  
  if (identical(methods, "All")) methods <- c("CHE","FEC","PET","PEESE")
  
  require(metafor)
  
  dat$var_avg <- unsplit(tapply(dat$Va, dat$studyid, mean), dat$studyid)
  V_mat <- vcalc(var_avg, 
                 cluster = studyid, obs = esid, 
                 data = dat, rho = rho)
  
  if ("CHE" %in% methods) {
    if (sum(table(dat$studyid) > 1L) >= 3) {
      res$che <- estimate_meta(
        dat = dat,
        rand = ~ 1 | studyid / esid,
        V_mat = V_mat,
        method = "CHE"
      )
    } else {
      res$che <- estimate_meta(
        dat = dat,
        rand = ~ 1 | studyid,
        V_mat = V_mat,
        method = "CHE"
      )
    }
  }
  
  if ("FEC" %in% methods) {
    res$fec <- estimate_meta(
      dat = dat,
      V_mat = V_mat,
      method = "FEC"
    )
  }
  
  
  if ("PET" %in% methods) {
    res$pet <- estimate_meta(
      dat = dat, 
      moderator = ~ sd_d,
      V_mat = V_mat,
      method = "PET"
    )
  }
  
  if ("PEESE" %in% methods) {
    res$peese <- estimate_meta(
      dat = dat, 
      moderator = ~ var_d,
      V_mat = V_mat,
      method = "PEESE"
    )
  }
  
  
  if (all(c("PET","PEESE") %in% methods)) {
    res$pet_peese <- if (res$pet$p < 0.1 & res$pet$est > 0) res$peese else res$pet
    res$pet_peese$estimator <- "PET/PEESE"
  }
  
  bind_rows(res)
}
