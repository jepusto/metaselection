library(tidyverse)

devtools::load_all()
source("simulation-step-bootstrap/2_estimation.R")
load("simulations-bootstrap/test/example_dat.RData")

estimate_step_models(dat = example_dat, 
                     estimators = "step-MLE",
                     steps = 0.025,
                     bootstrap = NULL)

estimate_step_models(dat = example_dat, 
                     estimators = "step-MLE",
                     steps = c(0.025, .5),
                     bootstrap = NULL)


estimate_step_models(dat = example_dat, 
                     estimators = "step-MLE",
                     steps = 0.025,
                     R = 9)

estimate_step_models(dat = example_dat, 
                     estimators = "step-hybrid-Broyden",
                     steps = 0.025,
                     R = 399)

estimate_step_models(dat = example_dat, 
                     estimators = "step-hybrid-Newton",
                     steps = 0.025,
                     R = 99)

estimate_step_models(dat = example_dat, 
                     estimators = "step-hybrid-rootSolve",
                     steps = 0.025,
                     R = 99)



dat <- example_dat
estimators = c("step-MLE")
steps <- 0.025
sd_char = "sd_d"
mean_mods = NULL
var_mods = NULL
sel_mods = NULL
sel_zero_mods = NULL
make_sandwich = TRUE
conf_level = .95
ML_optimizer = "BFGS"
ML_optimizer_control = list()
Broyden_optimizer_control = list()
Newton_optimizer_control = list()
rootSolve_optimizer_control = list()
bootstrap = "exp"
boot_CI = "percentile"
R = 9
sd <- as.symbol(sd_char)

estimate_step_models(dat = dat, 
                     estimators = "step-MLE",
                     steps = 0.025,
                     R = 9)

estimate_step_models(dat = dat, 
                     estimators = "step-hybrid-Broyden",
                     steps = 0.025,
                     R = 9)

estimate_step_models(dat = dat, 
                     estimators = "step-hybrid-Newton",
                     steps = 0.025,
                     R = 9)

estimate_step_models(dat = dat, 
                     estimators = "step-hybrid-rootSolve",
                     steps = 0.025,
                     R = 9)

# this runs
res_MLE <- 
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
      estimator = "ML",
      optimizer = ML_optimizer, 
      optimizer_control = ML_optimizer_control,
      bootstrap = bootstrap,
      boot_CI = boot_CI,
      R = R
    )


# this doesn't 
res_MLE <- 
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
    estimator = "ML",
    optimizer = ML_optimizer, 
    optimizer_control = ML_optimizer_control,
    bootstrap = "multinomial",
    boot_CI = boot_CI,
    R = R
  )




