library(tidyverse)

devtools::load_all()
load("simulations-bootstrap/test/example_dat.RData")
source("simulation-beta-function/2_estimation.R")


pkg_res <-
  selection_model(
    data = example_dat,
    yi = d,
    sei = sd_d,
    pi = p_onesided,
    cluster = studyid,
    steps = c(.025, .975),
    selection_type = "beta",
    estimator = "ML"
  )
pkg_res$est



estimate_step_models(dat = example_dat, R = 19)

pkg_boot <- 
  selection_model(
    data = example_dat,
    yi = d,
    sei = sd_d,
    pi = p_onesided,
    cluster = studyid,
    steps = c(.025, .975),
    selection_type = "beta",
    estimator = "ML",
    bootstrap = "exp",
    boot_CI = c("large-sample","percentile", "boot-t"),
    R = 49
  )

pkg_boot$est
