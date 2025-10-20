skip_if_not_installed("metadat")
skip_if_not_installed("metafor", minimum_version = "4.8-0")
library(metadat)
suppressPackageStartupMessages(library(metafor))

dat <- dat.lehmann2018
dat$sei <- sqrt(dat$vi)

test_that("step_loglik() and step_score() agree with metafor::selmodel().", {
  
  # Fixed effect selection model, no predictors
  FE1 <- rma.uni(yi = yi, sei = sei, data = dat, method = "FE")
  check_against_metafor_selmodel(FE1)
  
  # Fixed effect selection model, two steps
  check_against_metafor_selmodel(FE1, steps = c(.025, .50), tol_score = 2e-4)
  
  # Fixed effect selection model with meta-regression
  FE2 <- rma.uni(yi = yi, sei = sei, mods = ~ Color_Match + Gender + sqrt(Control_N), 
                 data = dat, method = "FE")
  check_against_metafor_selmodel(FE2, tol_score = 2e-4)
  
  # Fixed effect selection model with meta-regression, two steps
  check_against_metafor_selmodel(FE2, steps = c(.025, .50), tol_score = 2e-3)
  
  # Random effects 3-parameter selection model, no predictors
  RE1 <- rma.uni(yi = yi, sei = sei, data = dat, method = "ML")
  check_against_metafor_selmodel(RE1)
  
  # Random effects 4-parameter selection model, no predictors
  check_against_metafor_selmodel(
    RE1, steps = c(.025, .50), 
    tol_score = 1e-4, tol_param = 5e-4,
  )
  
  # Random effects 3-parameter selection model with meta-regression
  RE2 <- rma.uni(yi = yi, sei = sei, mods = ~ Color_Match + Gender, 
                 data = dat, method = "ML")
  check_against_metafor_selmodel(RE2)
  
  # Random effects 4-parameter selection model with meta-regression
  check_against_metafor_selmodel(RE2, steps = c(.025, .50), tol_score = 1e-4)
  
  # Random effects 3-parameter selection model with meta-regression
  RE3 <- rma.uni(yi = yi, sei = sei, mods = ~ Gender + sqrt(Control_N), 
                 data = dat, method = "ML")
  check_against_metafor_selmodel(RE3, tol_score = 5e-4, tol_param = 5e-4)
  
  # Random effects 4-parameter selection model with meta-regression
  check_against_metafor_selmodel(RE3, steps = c(.025, .50), tol_score = 5e-4)
  
})

test_that("selection_model() agrees with composites of metafor::selmodel().", {
  
  # fit 3PSM to subsets of data with metafor 
  dat_split <- by(dat, dat$Gender, identity)
  metafor_3PSM <- lapply(dat_split, \(dat) selmodel(rma.uni(yi = yi, sei = sei, data = dat, method = "ML"), 
                                                    type = "stepfun", steps = c(.025), control=list(optimizer = "nlminb", rel.tol = 1e-10)))
  metafor_3PSM_params <- 
    sapply(metafor_3PSM, get_selmodel_params) |>
    t() |>
    as.vector()
  
  mod_3PSM_fit <- 
    selection_model(
      data = dat, 
      yi = yi, sei = sei,
      steps = c(.025),
      mean_mods = ~ 0 + Gender,
      var_mods = ~ 0 + Gender,
      sel_mods = ~ 0 + Gender,
      priors = NULL
    )
  
  expect_equal(mod_3PSM_fit$est$Est, metafor_3PSM_params, tolerance = 1e-4)
  
  # fit moderated 4PSM to subsets of data with metafor 
  metafor_4PSM <- lapply(dat_split, \(dat) selmodel(rma.uni(yi = yi, sei = sei, mods = ~ Preregistered, data = dat, method = "ML"), 
                                                    type = "stepfun", steps = c(.025, .500), control=list(optimizer = "nlminb", rel.tol = 1e-10)))
  metafor_4PSM_params <- 
    sapply(metafor_4PSM, get_selmodel_params) |>
    t() |>
    as.vector()
  
  mod_4PSM_fit <- 
    selection_model(
      data = dat, 
      yi = yi, sei = sei,
      steps = c(.025, .500),
      mean_mods = ~ 0 + Gender + Gender:Preregistered,
      var_mods = ~ 0 + Gender,
      sel_mods = ~ 0 + Gender,
      priors = NULL
    )
  
  expect_equal(mod_4PSM_fit$est$Est, metafor_4PSM_params, tolerance = 1e-4)
  
})


test_that("beta_loglik() and beta_score() agree with metafor::selmodel().", {
  
  # Fixed effect selection model, no predictors
  FE1 <- rma.uni(yi = yi, sei = sei, data = dat, method = "FE")
  check_against_metafor_selmodel(
    FE1, type = "beta", steps = c(.01, .99), 
    tol_LRT = 1e-4, tol_score = 2e-3
  )
  
  # Fixed effect selection model with meta-regression
  FE2 <- rma.uni(yi = yi, sei = sei, mods = ~ Gender + sqrt(Control_N), 
                 data = dat, method = "FE")
  check_against_metafor_selmodel(
    FE2, type = "beta", steps = c(.025, .975), 
    tol_LRT = 1e-4, tol_score = 1e-2
  )
  
  # Random effects selection model, no predictors
  RE1 <- rma.uni(yi = yi, sei = sei, data = dat, method = "ML")
  check_against_metafor_selmodel(
    RE1, type = "beta", steps = NULL, 
    tol_LRT = 3e-3, tol_score = 5e-3, tol_SE = 5e-1
  )
  
  # Random effects selection model with meta-regression
  RE2 <- rma.uni(yi = yi, sei = sei, mods = ~ Color_Match + Gender, 
                 data = dat, method = "ML")
  check_against_metafor_selmodel(
    RE2, type = "beta", steps = c(.0001, .9999), 
    tol_LRT = 2e-3, tol_score = 2e-2, tol_SE = Inf
  )
  
})
