# che, fec, pet, and peese -----------------------------------------------------


estimate_meta  <- function(dat, 
                           rho = .8,
                           rand = NULL,
                           moderator = ~ 1,  # change to sd or var for pet peese
                           V_mat = NULL,
                           method) { 
  require(metafor)
  
  if (is.null(V_mat)) {
    V_mat <- vcalc(var_d, 
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

estimate_comparison_methods <- function(dat, 
                                        yi,
                                        vi, 
                                        study_id, 
                                        es_id, 
                                        rho = .8,
                                        methods = "All") {
  
  
  dat <- dat %>%
    mutate(d = {{yi}},
           var_d = {{vi}}, 
           studyid = {{study_id}}, 
           esid = {{es_id}}) %>%
    mutate(sd_d = sqrt(var_d))
  
  res <- list()
  
  if (is.null(methods)) return(res)
  
  if (identical(methods, "All")) methods <- c("CHE","FEC","PET","PEESE")
  
  require(metafor)
  
  # why did we do this?
  #dat$var_avg <- unsplit(tapply(dat$var_d, dat$studyid, mean), dat$studyid)
  
  V_mat <- vcalc(var_d, 
                 cluster = studyid, 
                 obs = esid, 
                 data = dat, 
                 rho = rho)
  
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

estimate_che <- function(dat, 
                         equation){

  V_mat <- vcalc(
    vi = vi, 
    cluster = study_id,
    obs = es_id, 
    data = dat,
    rho = 0.8,
    sparse = TRUE
  )
  
  
  
  res <- rma.mv(yi = as.formula(equation),
                V = V_mat,
                random = ~ 1 | study_id / es_id, 
                data = dat, 
                sparse = TRUE, 
                digits = 3) |> 
    robust(cluster = study_id, clubSandwich = TRUE)


return(res)

}

