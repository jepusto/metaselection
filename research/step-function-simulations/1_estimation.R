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
    estimators = c("CML","ARGL"),
    steps,
    mean_mods = NULL,
    var_mods = NULL,
    sel_mods = NULL,
    sel_zero_mods = NULL,
    priors = NULL,
    vcov_type = "robust",
    CI_type = c("large-sample","basic","percentile","student","bias-corrected","BCa"),
    bootstrap = "multinomial",
    conf_level = .95,
    CML_optimizer = NULL,
    CML_optimizer_control = list(),
    ARGL_optimizer_control = list(),
    R = 399,
    retry_bootstrap = 3L,
    use_jac = FALSE
) {
  

  if (bootstrap == "none") CI_type <- "large-sample"
  
  error_res <- list(
    est = data.frame(
            param = c("beta","gamma",paste0("zeta", seq_along(steps))),
            Est = rep(NA_real_, 2 + length(steps))
          ),
    info = list(convcode = NA_integer_, termcd = NA_integer_)
  )
  
  res <- list()
  
  if ("ARGL" %in% estimators) {
    
    res_hybrid <- tryCatch(
      selection_model(
        data = dat,
        yi = d,
        sei = sd_d,
        pi = p_onesided,
        cluster = studyid,
        selection_type = "step",
        steps = steps,
        mean_mods = mean_mods,
        var_mods = var_mods,
        sel_mods = sel_mods,
        sel_zero_mods = sel_zero_mods,
        priors = priors,
        estimator = "ARGL",
        vcov_type = vcov_type,
        CI_type = CI_type,
        conf_level = conf_level,
        optimizer = "nleqslv", 
        optimizer_control = c(method = "Broyden", ARGL_optimizer_control),
        bootstrap = bootstrap,
        R = R,
        retry_bootstrap = retry_bootstrap
      ), error = function(e) error_res
    )
    
    res_hybrid$est$estimator <- "ARGL"
    
    res$hybrid <- res_hybrid$est
  }
  
  if ("CML" %in% estimators) {
    theta <- if ("ARGL" %in% estimators) res$hybrid$est$Est else NULL
    res_MLE <- tryCatch(
      selection_model(
        data = dat,
        yi = d,
        sei = sd_d,
        pi = p_onesided,
        cluster = studyid,
        selection_type = "step",
        steps = steps,
        mean_mods = mean_mods,
        var_mods = var_mods,
        sel_mods = sel_mods,
        sel_zero_mods = sel_zero_mods,
        priors = priors,
        estimator = "CML",
        vcov_type = vcov_type,
        CI_type = CI_type,
        conf_level = conf_level,
        optimizer = CML_optimizer, 
        optimizer_control = CML_optimizer_control,
        bootstrap = bootstrap,
        R = R,
        retry_bootstrap = retry_bootstrap,
        theta = theta,
        use_jac = use_jac
      ), error = function(e) error_res
    )
    res_MLE$est$R_conv <- res_MLE$info$convcode
    res_MLE$est$max_method = res_MLE$method
    res_MLE$est$estimator <- "CML"
    
    res$MLE <- res_MLE$est
  }
  
  bind_rows(res)
}


# che, fec, pet, and peese -----------------------------------------------------


estimate_meta  <- function(dat, 
                           rho = .8,
                           rand = NULL,
                           sigma2 = NA,
                           moderator = ~ 1,  # change to sd or var for pet peese
                           V_mat = NULL,
                           W_mat = NULL) { 
  require(metafor)
  require(clubSandwich)

  if (is.null(V_mat)) {
    dat$var_avg <- unsplit(tapply(dat$Va, dat$studyid, mean), dat$studyid)
    
    V_mat <- vcalc(var_avg, 
                   cluster = studyid, obs = esid, 
                   data = dat, rho = rho, sparse = TRUE)  
  }
  
  optimizers <- c("nlminb","nloptr","Rvmmin","BFGS")
  mod <- "Non-converged"
  i <- 1L
  
  while (!inherits(mod, "rma.mv") & i <= 4L) {
    mod <- tryCatch(
      if (is.null(W_mat)) {
        rma.mv(
          yi = d,
          V = V_mat,
          random = rand,
          mods = moderator,
          data = dat,
          sigma2 = sigma2,
          sparse = TRUE,
          control = list(optimizer=optimizers[i])
        ) 
      } else {
        rma.mv(
          yi = d,
          V = V_mat,
          W = W_mat,
          random = rand,
          mods = moderator,
          data = dat,
          sigma2 = sigma2,
          sparse = TRUE,
          control = list(optimizer=optimizers[i])
        )
      },
      error = function(e) "Non-converged"
    )
    i <- i + 1L
  }
 
  return(mod)
}

tidy_meta_results <- function(mod, dat,method) {
  if (inherits(mod, "rma.mv")) {
    
    est <- tryCatch(
      conf_int(mod, cluster = dat$studyid, vcov = "CR2", coefs = 1L, p_values = TRUE),
      error = function(e) "NPD"
    )
    
    if (inherits(est, "conf_int_clubSandwich")) {
      res <- data.frame(
        param = c("beta","gamma"),
        Est = c(est$beta, log(sum(mod$sigma2))),
        SE = c(est$SE, NA_real_),
        df = c(est$df, NA_real_),
        CI_lo = c(est$CI_L, NA_real_),
        CI_hi = c(est$CI_U, NA_real_),
        p_value = c(est$p_val, NA_real_),
        R_conv = NA,
        estimator = method
      )
      return(res)
    }
    
  }
    
  return(NULL)
}

estimate_comparison_methods <- function(dat, rho = .8, methods = "All") {
  
  res <- list()
  
  if (is.null(methods) || methods == "None") return(data.frame())
  
  if (identical(methods, "All")) methods <- c("CHE","CHE-ISCW","PET","PEESE")
  if ("CHE-ISCW" %in% methods) methods <- union(methods, "CHE")
  
  require(metafor)
  
  dat$var_avg <- unsplit(tapply(dat$Va, dat$studyid, mean), dat$studyid)
  V_mat <- vcalc(
    var_avg, cluster = studyid, obs = esid, 
    data = dat, rho = rho, sparse = TRUE
  )
  
  if ("CHE" %in% methods) {
    
    if (sum(table(dat$studyid) > 1L) >= 3) {
      che_mod <- estimate_meta(
        dat = dat,
        rand = ~ 1 | studyid / esid,
        V_mat = V_mat
      )
      
    } else {
      che_mod <- estimate_meta(
        dat = dat,
        rand = ~ 1 | studyid,
        V_mat = V_mat
      )
    }
    
    res$che <- tidy_meta_results(che_mod, dat = dat, method = "CHE")
  }
  
  if ("CHE-ISCW" %in% methods) {
    
    sigma2 <- che_mod$sigma2
    W_mat <- solve(V_mat)
    
    if (sum(table(dat$studyid) > 1L) >= 3) {
      iscw_mod <- estimate_meta(
        dat = dat,
        rand = ~ 1 | studyid / esid,
        sigma2 = sigma2,
        V_mat = V_mat,
        W_mat = W_mat
      )
      
    } else {
      iscw_mod <- estimate_meta(
        dat = dat,
        rand = ~ 1 | studyid,
        sigma2 = sigma2,
        V_mat = V_mat,
        W_mat = W_mat
      )
    }
    
    res$fec <- tidy_meta_results(iscw_mod, dat = dat, method =  "CHE-ISCW")
  }
  
  if ("PET" %in% methods) {
    res$pet <- estimate_meta(
      dat = dat, 
      moderator = ~ sd_d,
      V_mat = V_mat
    ) |>
      tidy_meta_results(dat = dat, method = "PET")
  }
  
  if ("PEESE" %in% methods) {
    res$peese <- estimate_meta(
      dat = dat, 
      moderator = ~ var_d,
      V_mat = V_mat
    ) |>
      tidy_meta_results(dat = dat, method = "PEESE")
  }
  
  
  if (all(c("PET","PEESE") %in% methods) && !is.null(res$pet) && !is.null(res$peese)) {
    res$pet_peese <- if (res$pet$p_value[1] < 0.1 & res$pet$est[1] > 0) res$peese else res$pet
    res$pet_peese$estimator <- "PET/PEESE"
  }
  
  bind_rows(res)
}
